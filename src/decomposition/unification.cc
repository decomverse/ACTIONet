#include "ACTIONet.h"
//#include "dagConstruct.h"

// vector<double> Corrector::vals;

namespace ACTIONet
{
  mat unsigned_cluster_batch(sp_mat A, vec resolutions,
                             uvec initial_clusters = uvec(), int seed = 0);
  mat signed_cluster_batch(sp_mat A, vec resolutions,
                           uvec initial_clusters = uvec(), int seed = 0);

  double Kappa(double p, double q)
  {
    double a = 0.0, b = 0.0;
    if ((1e-300 < p) & (1e-300 < q))
    {
      a = p * log(p / q);
    }
    if ((p < (1 - 1e-300)) & (q < (1 - 1e-300)))
    {
      b = (1 - p) * log((1 - p) / (1 - q));
    }

    double k = a + b;
    return (k);
  }

  double log_HGT_tail(int population_size, int success_count, int sample_size,
                      int observed_success)
  {
    if (observed_success == 0)
      return (0);

    double success_rate = (double)success_count / population_size;
    double expected_success = sample_size * success_rate;
    double delta = (observed_success / expected_success) - 1.0;
    if (delta < 0)
    {
      return (0);
    }

    double log_tail_bound =
        sample_size * Kappa((1.0 + delta) * success_rate, success_rate);

    return (log_tail_bound);
  }

  double assess_overlap(uvec i1, uvec i2, int population_size)
  {
    int success_count = i1.n_elem;
    int sample_size = i2.n_elem;

    uvec shared = intersect(i1, i2);
    int observed_success = shared.n_elem;

    double log_pval = log_HGT_tail(population_size, success_count, sample_size,
                                   observed_success);

    return (log_pval);
  }

  mat compute_overlap_matrix(mat C)
  {
    int N = C.n_cols;
    mat O = zeros(N, N);

    vector<uvec> indices(N);
    for (int i = 0; i < N; i++)
    {
      uvec idx = find(C.col(i) > 0);
      indices[i] = idx;
    }

    for (int i = 0; i < N; i++)
    {
      uvec i1 = indices[i];
      for (int j = i + 1; j < N; j++)
      {
        uvec i2 = indices[j];

        O(i, j) = O(j, i) = assess_overlap(i1, i2, N);
      }
    }

    return (O);
  }

  field<vec> run_HDBSCAN(mat &X, int minPoints = 5, int minClusterSize = 5)
  {
    Hdbscan hdbscan(X);
    hdbscan.execute(minPoints, minClusterSize, "Euclidean");

    vec labels(X.n_rows);
    vec membershipProbabilities(X.n_rows);
    vec outlierScores(X.n_rows);

    for (int i = 0; i < X.n_rows; i++)
    {
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

  mat NetEnh(mat A)
  {
    A.diag().zeros();
    mat P = normalise(A, 1, 1);

    mat D = diagmat(1.0 / (sqrt(sum(P) + 1e-16)));
    mat W = P * D;
    mat P2 = W * trans(W);

    return (P2);
  }

  field<mat> nndsvd(mat &A, int dim = 100, int max_iter = 5)
  {
    dim = std::min(dim, (int)A.n_cols);
    field<mat> SVD_res = HalkoSVD(A, dim, max_iter, 0, 0);

    mat U = SVD_res(0);
    vec s = SVD_res(1);
    mat V = SVD_res(2);

    mat W(size(U));
    mat H(size(V));

    uvec mask_idx;
    for (int i = 0; i < U.n_cols; i++)
    {
      vec u = U.col(i);
      vec v = V.col(i);

      vec up = u;
      mask_idx = find(u < 0);
      if (mask_idx.n_elem > 0)
        up(mask_idx).zeros();

      vec un = -u;
      mask_idx = find(u > 0);
      if (mask_idx.n_elem > 0)
        un(mask_idx).zeros();

      vec vp = v;
      mask_idx = find(v < 0);
      if (mask_idx.n_elem > 0)
        vp(mask_idx).zeros();

      vec vn = -v;
      mask_idx = find(v > 0);
      if (mask_idx.n_elem > 0)
        vn(mask_idx).zeros();

      double n_up = norm(up);
      double n_un = norm(un);
      double n_vp = norm(vp);
      double n_vn = norm(vn);

      double termp = n_up * n_vp;
      double termn = n_un * n_vn;
      if (termp >= termn)
      {
        W.col(i) = sqrt(s(i) * termp) * up / n_up;
        H.col(i) = sqrt(s(i) * termp) * vp / n_vp;
      }
      else
      {
        W.col(i) = sqrt(s(i) * termn) * un / n_un;
        H.col(i) = sqrt(s(i) * termn) * vn / n_vn;
      }
    }

    field<mat> out(5);
    out(0) = W;
    out(1) = H;
    out(2) = s;
    out(3) = U;
    out(4) = V;

    return (out);
  }

  field<mat> orient_SVD(field<mat> &SVD_out)
  {
    mat U = SVD_out(0);
    vec s = SVD_out(1);
    mat V = SVD_out(2);

    int dim = U.n_cols;

    mat Up = U;
    mat Vp = V;

    uvec mask_idx;
    for (int i = 0; i < U.n_cols; i++)
    {
      vec u = U.col(i);
      vec v = V.col(i);

      vec up = u;
      mask_idx = find(u < 0);
      if (mask_idx.n_elem > 0)
        up(mask_idx).zeros();

      vec un = -u;
      mask_idx = find(u > 0);
      if (mask_idx.n_elem > 0)
        un(mask_idx).zeros();

      vec vp = v;
      mask_idx = find(v < 0);
      if (mask_idx.n_elem > 0)
        vp(mask_idx).zeros();

      vec vn = -v;
      mask_idx = find(v > 0);
      if (mask_idx.n_elem > 0)
        vn(mask_idx).zeros();

      double n_up = norm(up);
      double n_un = norm(un);
      double n_vp = norm(vp);
      double n_vn = norm(vn);

      double termp = n_up * n_vp;
      double termn = n_un * n_vn;
      if (termp < termn)
      {
        Up.col(i) *= -1;
        Vp.col(i) *= -1;
      }
    }

    field<mat> out(3);
    out(0) = Up;
    out(1) = s;
    out(2) = Vp;

    return (out);
  }

  field<mat> convexSVD(mat &A, int dim = 100, int max_iter = 5)
  {
    field<mat> out(4);

    dim = std::min(dim, (int)A.n_cols);
    field<mat> SVD_res = HalkoSVD(A, dim, max_iter, 0, 0);
    SVD_res = orient_SVD(SVD_res);

    mat U = SVD_res(0);
    vec s = SVD_res(1);
    mat V = SVD_res(2);

    out(0) = U;
    out(1) = s;
    out(2) = V;

    return (out);
  }

  field<mat> recursiveNMU(mat M, int dim = 100, int max_SVD_iter = 5,
                          int max_iter_inner = 40)
  {
    dim = std::min(dim, (int)M.n_cols);

    mat W(M.n_rows, dim);
    mat H(M.n_cols, dim);
    vec obj(dim);
    vec factor_weights(dim);
    vec ss(dim);

    // sparse M_sp = sparse(M);

    double denom = sum(sum(square(M)));
    for (int k = 0; k < dim; k++)
    {
      /*
      mat U, V;
      vec s;
      svds( U, s, V, M_sp, 1, tol );
      */

      field<mat> SVD_res = HalkoSVD(M, 1, max_SVD_iter, 0, 0);
      mat U = SVD_res(0);
      vec s = SVD_res(1);
      mat V = SVD_res(2);

      vec x = abs(U.col(0)) * sqrt(s(0));
      vec y = abs(V.col(0)) * sqrt(s(0));

      W.col(k) = x;
      H.col(k) = y;

      mat R = M - x * trans(y);
      mat lambda = -R;
      lambda.transform([](double val)
                       { return (val < 0 ? 0 : val); });

      for (int j = 0; j < max_iter_inner; j++)
      {
        mat A = M - lambda;

        x = A * y;
        x.transform([](double val)
                    { return (val < 0 ? 0 : val); });
        x /= (max(x) + 1e-16);

        y = trans(trans(x) * A);
        y.transform([](double val)
                    { return (val < 0 ? 0 : val); });
        y /= dot(x, x);

        if ((sum(x) != 0) && (sum(y) != 0))
        {
          W.col(k) = x;
          H.col(k) = y;
          R = M - x * trans(y);
          lambda = lambda - R / ((double)j + 1);
          lambda.transform([](double val)
                           { return (val < 0 ? 0 : val); });
        }
        else
        {
          lambda /= 2.0;
          x = W.col(k);
          y = H.col(k);
        }
      }

      mat oldM = M;
      M -= x * trans(y);
      M.transform([](double val)
                  { return (val < 0 ? 0 : val); });

      obj(k) = (sum(sum(square(M))) / denom);
      ss(k) = s(0);

      double w_norm1 = sum(abs(W.col(k))); // abs() is reducndant
      W.col(k) /= w_norm1;
      double h_norm1 = sum(abs(H.col(k))); // abs() is reducndant
      H.col(k) /= h_norm1;

      factor_weights(k) = w_norm1 * h_norm1;
    }

    field<mat> out(5);
    out(0) = W;
    out(1) = H;
    out(2) = factor_weights;
    out(3) = obj;
    out(4) = ss;

    return (out);
  }

  field<mat> recursiveNMU_mine(mat M, int dim = 100, int max_SVD_iter = 1000,
                               int max_iter_inner = 100)
  {
    dim = std::min(dim, (int)M.n_cols);

    mat W(M.n_rows, dim);
    mat H(M.n_cols, dim);
    vec factor_weights(dim);

    mat M0 = M;

    vec s;
    mat U, V;
    double denom = sum(sum(square(M)));
    for (int k = 0; k < dim; k++)
    {
      field<mat> SVD_res = IRLB_SVD(M, 1, max_SVD_iter, 0, 0);
      U = SVD_res(0);
      s = SVD_res(1);
      V = SVD_res(2);

      vec w = trans(trans(U.col(0)) * M);
      int selected_columns = index_max(w);

      vec u = M0.col(selected_columns);
      W.col(k) = u;
      H.col(k).zeros();
      H(selected_columns, k) = 1;

      vec v = M.col(selected_columns);
      M -= v * (trans(v) * M) / dot(v, v);
      M.transform([](double val)
                  { return (val < 0 ? 0 : val); });

      factor_weights(k) = (sum(sum(abs(M))) / M.n_cols);
    }

    field<mat> out(3);
    out(0) = W;
    out(1) = H;
    out(2) = factor_weights;

    return (out);
  }

  vec sweepcut(sp_mat &A, vec s, int min_size, int max_size)
  {
    A.diag().zeros();
    int nV = A.n_rows;
    if (max_size == -1)
    {
      max_size = nV;
    }

    vec w = vec(sum(A, 1));
    double total_vol = sum(w);

    vec conductance = datum::inf * ones(max_size);

    uvec perm = sort_index(s, "descend");
    vec x = zeros(nV);
    x(perm(span(0, min_size - 1))).ones();
    double vol = sum(w(perm(span(0, min_size - 1))));

    double cut_size = vol;
    for (int i = 0; i < min_size; i++)
    {
      for (int j = 0; j < min_size; j++)
      {
        cut_size -= A(i, j);
      }
    }

    for (int i = min_size; i <= max_size - 1; i++)
    {
      int u = perm(i);
      vol += w[u];

      x(u) = 1;

      sp_mat::col_iterator it = A.begin_col(u);
      for (; it != A.end_col(u); it++)
      {
        int v = it.row();
        if (x[v] == 0)
        { // v is in S_prime (not yet selected)
          cut_size += (*it);
        }
        else
        {
          cut_size -= (*it);
        }
      }

      double vol_prime = total_vol - vol;
      conductance(i) = cut_size / min(vol, vol_prime);
    }

    return (conductance);
  }

  unification_results unify_archetypes(mat &S_r, mat &C_stacked, mat &H_stacked,
                                       double backbone_density, double resolution,
                                       int min_cluster_size, int thread_no,
                                       int normalization)
  {
    if (thread_no <= 0)
    {
      thread_no = SYS_THREADS_DEF;
    }

    stdout_printf("Unifying %d archetypes (%d threads):\n", C_stacked.n_cols,
                  thread_no);
    // stdout_printf("\tParameters: backbone density = %0.2f, outlier threshold
    // (on core #) = %d\n", backbone_density, outlier_threshold);
    FLUSH;

    unification_results output;

    H_stacked = normalise(H_stacked, 1, 0);
    sp_mat H_stacked_sp = sp_mat(H_stacked);
    mat H_arch = spmat_mat_product_parallel(H_stacked_sp, C_stacked, thread_no);
    H_arch.replace(datum::nan, 0); // replace each NaN with 0

    SPA_results SPA_out = run_SPA(H_arch, H_arch.n_cols);
    uvec candidates = SPA_out.selected_columns;
    vec scores = SPA_out.column_norms;
    double x1 = sum(scores);
    double x2 = sum(square(scores));
    double arch_no = round((x1 * x1) / x2);
    candidates = candidates(span(0, arch_no - 1));

    stdout_printf("# unified archetypes: %d\n", (int)arch_no);
    output.selected_archetypes = candidates;
    /*
      double M = 16, ef_construction = 200, ef = 200;
      sp_mat backbone =
          ACTIONet::buildNetwork(H_arch, "k*nn", "jsd", backbone_density, 1, M,
                                 ef_construction, ef, true, 10);
      output.dag_adj = backbone;

      vec initial_clusters = zeros(backbone.n_rows);
      for (int i = backbone.n_rows - 1; i >= 0; i--) {
        int v = candidates(i);
        uvec N = find(vec(backbone.col(v)) > 0);

        initial_clusters(v) = i;
        initial_clusters(N) = i * ones(N.n_elem);
      }

      uvec initial_clusters_uvec = conv_to<uvec>::from(initial_clusters);
      vec arch_clusters = ACTIONet::unsigned_cluster(backbone, resolution,
                                                     initial_clusters_uvec, 0);
      int idx = 1;
      for (; idx <= max(arch_clusters); idx++) {
        uvec selected_archs = find(arch_clusters == idx);
        int cc = selected_archs.n_elem;
        if (cc < min_cluster_size) {
          break;
        }
      }
      int arch_no = idx - 1;

      stdout_printf("# unified archetypes: %d\n", arch_no);
      output.selected_archetypes = conv_to<uvec>::from(arch_clusters - 1);

      // Compute unified archetypes
      vec c;
      mat C_unified = zeros(C_stacked.n_rows, arch_no);
      for (idx = 1; idx <= arch_no; idx++) {
        uvec cur_archs = find(arch_clusters == idx);
        if (cur_archs.n_elem == 1) {
          c = C_stacked.col(cur_archs(0));
        } else {
          c = mean(C_stacked.cols(cur_archs), 1);
        }
        c = c / sum(c);
        C_unified.col(idx - 1) = c;
      }
    */
    mat C_unified = C_stacked.cols(candidates);

    mat X_r = normalize_mat(S_r, normalization, 0);

    mat W_r_unified = X_r * C_unified;

    mat H_unified = run_simplex_regression(W_r_unified, X_r, false);

    uvec assigned_archetypes = trans(index_max(H_unified, 0));

    /*
        // Construct Ontology!
        mat Sim = 1+cor(W_r_unified);

        graph_undirected inputNetwork(Sim);
        DAGraph ontology;
        ontology.setTerminalName("archetype");
        nodeDistanceObject nodeDistances;

        double threshold = 0.05;
        double density = 0.5;

        dagConstruct::constructDAG(inputNetwork, ontology, nodeDistances,
       threshold, density);

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
    // output.arch_membership_weights = H_NMU.cols(selected_cols);

    FLUSH;
    return (output);
  }

} // namespace ACTIONet
