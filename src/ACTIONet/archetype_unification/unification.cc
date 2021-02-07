#include "ACTIONet.h"
#include "dagConstruct.h"

vector<double> Corrector::vals;

template <class Function>
inline void ParallelFor(size_t start, size_t end, size_t thread_no,
                        Function fn) {
  if (thread_no <= 0) {
    thread_no = std::thread::hardware_concurrency();
  }

  if (thread_no == 1) {
    for (size_t id = start; id < end; id++) {
      fn(id, 0);
    }
  } else {
    std::vector<std::thread> threads;
    std::atomic<size_t> current(start);

    // keep track of exceptions in threads
    // https://stackoverflow.com/a/32428427/1713196
    std::exception_ptr lastException = nullptr;
    std::mutex lastExceptMutex;

    for (size_t threadId = 0; threadId < thread_no; ++threadId) {
      threads.push_back(std::thread([&, threadId] {
        while (true) {
          size_t id = current.fetch_add(1);

          if ((id >= end)) {
            break;
          }

          try {
            fn(id, threadId);
          } catch (...) {
            std::unique_lock<std::mutex> lastExcepLock(lastExceptMutex);
            lastException = std::current_exception();
            /*
             * This will work even when current is the largest value that
             * size_t can fit, because fetch_add returns the previous value
             * before the increment (what will result in overflow
             * and produce 0 instead of current + 1).
             */
            current = end;
            break;
          }
        }
      }));
    }
    for (auto &thread : threads) {
      thread.join();
    }
    if (lastException) {
      std::rethrow_exception(lastException);
    }
  }
}

namespace ACTIONet {
mat unsigned_cluster_batch(sp_mat A, vec resolutions,
                           uvec initial_clusters = uvec(), int seed = 0);
mat signed_cluster_batch(sp_mat A, vec resolutions,
                         uvec initial_clusters = uvec(), int seed = 0);

double Kappa(double p, double q) {
  double a = 0.0, b = 0.0;
  if ((1e-300 < p) & (1e-300 < q)) {
    a = p * log(p / q);
  }
  if ((p < (1 - 1e-300)) & (q < (1 - 1e-300))) {
    b = (1 - p) * log((1 - p) / (1 - q));
  }

  double k = a + b;
  return (k);
}

double log_HGT_tail(int population_size, int success_count, int sample_size,
                    int observed_success) {
  if (observed_success == 0) return (0);

  double success_rate = (double)success_count / population_size;
  double expected_success = sample_size * success_rate;
  double delta = (observed_success / expected_success) - 1.0;
  if (delta < 0) {
    return (0);
  }

  double log_tail_bound =
      sample_size * Kappa((1.0 + delta) * success_rate, success_rate);

  return (log_tail_bound);
}

double assess_overlap(uvec i1, uvec i2, int population_size) {
  int success_count = i1.n_elem;
  int sample_size = i2.n_elem;

  uvec shared = intersect(i1, i2);
  int observed_success = shared.n_elem;

  double log_pval = log_HGT_tail(population_size, success_count, sample_size,
                                 observed_success);

  return (log_pval);
}

mat compute_overlap_matrix(mat C) {
  int N = C.n_cols;
  mat O = zeros(N, N);

  vector<uvec> indices(N);
  for (int i = 0; i < N; i++) {
    uvec idx = find(C.col(i) > 0);
    indices[i] = idx;
  }

  for (int i = 0; i < N; i++) {
    uvec i1 = indices[i];
    for (int j = i + 1; j < N; j++) {
      uvec i2 = indices[j];

      O(i, j) = O(j, i) = assess_overlap(i1, i2, N);
    }
  }

  return (O);
}

field<vec> run_HDBSCAN(mat &X, int minPoints = 5, int minClusterSize = 5) {
  Hdbscan hdbscan(X);
  hdbscan.execute(minPoints, minClusterSize, "Euclidean");

  vec labels(X.n_rows);
  vec membershipProbabilities(X.n_rows);
  vec outlierScores(X.n_rows);

  for (int i = 0; i < X.n_rows; i++) {
    labels[i] = hdbscan.labels_[i];
    membershipProbabilities[i] = hdbscan.membershipProbabilities_[i];
    outlierScores[hdbscan.outlierScores_[i].id] =
        hdbscan.outlierScores_[i].score;
  }

  field<vec> out(3);
  out(0) = labels;
  out(1) = membershipProbabilities;
  out(2) = outlierScores;

  return (out);
}

mat NetEnh(mat A) {
  A.diag().zeros();
  mat P = normalise(A, 1, 1);

  mat D = diagmat(sqrt(sum(P) + 1e-16));
  mat W = P * D;
  mat P2 = W * trans(W);
  // P2.diag().zeros();

  return (P2);
}

field<mat> nndsvd(mat &A, int dim = 100, int max_iter = 5) {
  dim = std::min(dim, (int)A.n_cols);
  field<mat> SVD_res = HalkoSVD(A, dim, max_iter, 0);

  mat U = SVD_res(0);
  vec s = SVD_res(1);
  mat V = SVD_res(2);

  mat W(size(U));
  mat H(size(V));

  W.col(0) = abs(U.col(0)) * sqrt(s(0));
  H.col(0) = abs(V.col(0)) * sqrt(s(0));
  uvec mask_idx;
  for (int i = 1; i < U.n_cols; i++) {
    vec u = U.col(i);
    vec v = V.col(i);

    vec up = u;
    mask_idx = find(u < 0);
    if (mask_idx.n_elem > 0) up(mask_idx).zeros();

    vec un = -u;
    mask_idx = find(u > 0);
    if (mask_idx.n_elem > 0) un(mask_idx).zeros();

    vec vp = v;
    mask_idx = find(v < 0);
    if (mask_idx.n_elem > 0) vp(mask_idx).zeros();

    vec vn = -v;
    mask_idx = find(v > 0);
    if (mask_idx.n_elem > 0) vn(mask_idx).zeros();

    double n_up = norm(up);
    double n_un = norm(un);
    double n_vp = norm(vp);
    double n_vn = norm(vn);

    double termp = n_up * n_vp;
    double termn = n_un * n_vn;
    if (termp >= termn) {
      W.col(i) = sqrt(s(i) * termp) * up / n_up;
      H.col(i) = sqrt(s(i) * termp) * vp / n_vp;
    } else {
      W.col(i) = sqrt(s(i) * termn) * un / n_un;
      H.col(i) = sqrt(s(i) * termn) * vn / n_vn;
    }
  }

  field<mat> out(3);
  out(0) = W;
  out(1) = H;
  out(2) = s;

  return (out);
}

unification_results unify_archetypes(sp_mat &G, mat &S_r, mat &C_stacked,
                                     double alpha = 0.99,
                                     double sensitivity = 0,
                                     int thread_no = 0) {
  if (thread_no <= 0) {
    thread_no = SYS_THREADS_DEF;
  }
  stdout_printf("Unifying %d archetypes (%d threads):\n", C_stacked.n_cols,
                thread_no);

  stdout_printf("\tParameters: alpha = %0.2f, sensitivity = %0.2f\n", alpha,
                sensitivity);
  FLUSH;

  unification_results output;

  // Smooth archetypes using ACTIONet
  printf("Running network diffusion\n");
  mat X0 = join_rows(ones(C_stacked.n_rows), C_stacked);
  X0 = normalise(X0, 1, 0);
  sp_mat X0_sp = sp_mat(X0);
  mat PR = compute_network_diffusion_fast(G, X0_sp, thread_no, alpha, 5);

  mat C_imputed = PR.cols(1, PR.n_cols - 1);

  mat LOR = C_imputed;
  for (int i = 0; i < LOR.n_cols; i++) {
    LOR.col(i) = LOR.col(i) / PR.col(0);
  }
  LOR = log2(LOR);
  uvec zero_idx = find(C_imputed == 0);
  LOR(zero_idx).zeros();
  LOR.transform([](double val) { return (val < 0 ? 0 : val); });
  mat LOR_norm = normalise(LOR, 1, 0);

  field<mat> nndsvd_out = nndsvd(LOR_norm, min((int)LOR_norm.n_cols, 100));
  mat Wpos = nndsvd_out(0);
  mat Hpos = nndsvd_out(1);

  SPA_results res = run_SPA(Wpos, Wpos.n_cols);
  uvec selected_columns = res.selected_columns;
  vec x = res.column_norms;

  // Select number of unified archetypes
  int state_no = -1;
  vec coverage = cumsum(x) / sum(x);
  if (sensitivity == 0) {
    double sx = sum(x);
    double sx_sq = sum(square(x));
    state_no = (int)round(sx * sx / sx_sq);
  } else if (sensitivity <= 1) {
    state_no = min(find(sensitivity <= coverage));
  } else {
    state_no = min((int)selected_columns.n_elem, (int)sensitivity);
  }
  printf("%d unified states (sensitivity = %.2f)\n", state_no,
         coverage(state_no - 1));

  /*
  mat subH = Hpos.cols(selected_columns(span(0, state_no-1)));
  mat weights = normalise(subH, 1, 0);
  */
  mat A = normalise(LOR, 1, 0);
  mat B = normalise(Wpos.cols(selected_columns(span(0, state_no - 1))), 1, 0);
  mat weights = run_simplex_regression(A, B, false);

  // Compute unified archetypes
  mat C_unified = C_imputed * weights;
  mat W_r_unified = S_r * C_unified;

  mat H_unified = run_simplex_regression(W_r_unified, S_r, false);
  uvec assigned_archetypes = trans(index_max(H_unified, 0));

  /*
    // Construct Ontology!
    graph_undirected inputNetwork(Sim);
    DAGraph ontology;
    ontology.setTerminalName("archetype");
    nodeDistanceObject nodeDistances;

    double threshold = 0.05;
    double density = 0.5;

    dagConstruct::constructDAG(inputNetwork, ontology, nodeDistances, threshold,
                               density);

    vector<vector<int>> dag_nodes;
    vector<int> dag_nodes_type;

    for (int k = 0; k < Sim.n_rows; k++) {
      vector<int> v;
      v.push_back(k);
      dag_nodes.push_back(v);
      dag_nodes_type.push_back(0);
    }

    for (map<pair<unsigned, unsigned>, string>::iterator edgesIt =
             ontology.edgesBegin();
         edgesIt != ontology.edgesEnd(); ++edgesIt) {
      unsigned ii = edgesIt->first.first;
      unsigned jj = edgesIt->first.second;

      vector<int> v;
      int tt;
      if (edgesIt->second == "archetype") {
        v.push_back(jj);
        tt = 0;
      } else {  // Internal node
        for (int kk = 0; kk < dag_nodes[jj].size(); kk++) {
          v.push_back(dag_nodes[jj][kk]);
        }
        tt = 1;
      }

      if (ii >= dag_nodes.size()) {
        dag_nodes.push_back(v);
        dag_nodes_type.push_back(tt);
      } else {  // merge
        for (int kk = 0; kk < v.size(); kk++) {
          dag_nodes[ii].push_back(v[kk]);
        }

        if (edgesIt->second != "archetype") dag_nodes_type[ii] = 1;
      }

      // cout << ontology.getName(edgesIt->first.first) << "\t" <<
      // ontology.getName(edgesIt->first.second) << "\t" << edgesIt->second <<
      // "\t"
      // << ontology.getWeight(edgesIt->first.first) << endl;
    }

    // Get internal adjacency matrix of DAGs
    int dag_node_counts = dag_nodes.size();
    mat dag_adj = zeros(dag_node_counts, dag_node_counts);
    vec dag_node_annotations = zeros(dag_node_counts);
    for (map<pair<unsigned, unsigned>, string>::iterator edgesIt =
             ontology.edgesBegin();
         edgesIt != ontology.edgesEnd(); ++edgesIt) {
      unsigned ii = edgesIt->first.first;
      unsigned jj = edgesIt->first.second;
      double w = ontology.getWeight(edgesIt->first.first);

      if (edgesIt->second == "archetype") {
        dag_node_annotations(ii) = 1;
      } else {
        dag_node_annotations(ii) = 2;
      }

      dag_adj(ii, jj) = 1;
    }

    output.dag_adj = dag_adj;
    output.dag_node_annotations = dag_node_annotations;
  */

  output.C_unified = C_unified;
  output.H_unified = H_unified;
  output.assigned_archetypes = assigned_archetypes;
  // output.selected_archetypes = selected_archetypes(idx);

  FLUSH;
  return (output);
}

}  // namespace ACTIONet
