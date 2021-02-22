#include <ACTIONet.h>

#include <hnswlib.h>
#include <atomic>
#include <thread>

template <class Function>
inline void ParallelFor(size_t start, size_t end, size_t thread_no,
                        Function fn) {
  if (thread_no <= 0) {
    thread_no = SYS_THREADS_DEF;
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
double Sim(const double *pVect1, const double *pVect2, const double *log_vec,
           int N) {
  double half = 0.5;

  double sum1 = 0, sum2 = 0;
  for (size_t i = 0; i < N; i++) {
    double p = pVect1[i];
    double q = pVect2[i];
    double m = (p + q) * half;

    int p_idx = (int)floor(p * 1000000.0);
    int q_idx = (int)floor(q * 1000000.0);
    int m_idx = (int)floor(m * 1000000.0);

    double lg_p = log_vec[p_idx];
    double lg_q = log_vec[q_idx];
    double lg_m = log_vec[m_idx];

    sum1 += (p * lg_p) + (q * lg_q);
    sum2 += m * lg_m;
  }

  double JS = std::max(half * sum1 - sum2, 0.0);

  return (double)(1.0 - sqrt(JS));
}

mat computeFullSim(mat &H, int thread_no) {
  double log_vec[1000001];
  for (int i = 0; i <= 1000000; i++) {
    log_vec[i] = (double)log2((double)i / 1000000.0);
  }
  log_vec[0] = 0;

  H = clamp(H, 0, 1);
  H = normalise(H, 1, 0);  // make the norm (sum) of each column 1

  int sample_no = H.n_cols;
  int dim = H.n_rows;
  // stdout_printf("sample # = %d, dim = %d\n", sample_no, dim);

  mat G = zeros(sample_no, sample_no);
  ParallelFor(0, sample_no, thread_no, [&](size_t i, size_t threadId) {
    for (int j = 0; j < sample_no; j++) {
      G(i, j) = Sim(H.colptr(i), H.colptr(j), log_vec, dim);
    }
  });

  G = clamp(G, 0.0, 1.0);

  return (G);
}

// k^{*}-Nearest Neighbors: From Global to Local (NIPS 2016)
sp_mat build_ACTIONet_JS_KstarNN(mat H_stacked, double density = 1.0,
                                 int thread_no = 0, double M = 16,
                                 double ef_construction = 200, double ef = 10,
                                 bool mutual_edges_only = true) {
  double LC = 1.0 / density;

  if (thread_no <= 0) {
    thread_no = SYS_THREADS_DEF;  // std::thread::hardware_concurrency();
  }

  stdout_printf("Building adaptive network (density = %.2f)\n", density);

  H_stacked = clamp(H_stacked, 0, 1);
  H_stacked =
      normalise(H_stacked, 1, 0);  // make the norm (sum) of each column 1

  double kappa = 5.0;
  int sample_no = H_stacked.n_cols;
  int kNN = min(
      sample_no - 1,
      (int)(kappa *
            round(sqrt(sample_no))));  // start with uniform k=sqrt(N) ["Pattern
                                       // Classification" book by Duda et al.]

  int dim = H_stacked.n_rows;
  int max_elements = H_stacked.n_cols;
  hnswlib::JSDSpace *space = new hnswlib::JSDSpace(dim);
  hnswlib::HierarchicalNSW<float> *appr_alg =
      new hnswlib::HierarchicalNSW<float>(space, max_elements, M,
                                           ef_construction);
  appr_alg->setEf(ef);

  stdout_printf("\tBuilding index ... ");
  fmat X = conv_to<fmat>::from(H_stacked);    
  ParallelFor(0, max_elements, thread_no, [&](size_t j, size_t threadId) {
    appr_alg->addPoint(X.colptr(j), static_cast<size_t>(j));
  });
  stdout_printf("done\n");

  stdout_printf("\tIdentifying nearest neighbors ... ");
  mat idx = zeros(sample_no, kNN + 1);
  mat dist = zeros(sample_no, kNN + 1);
  //		for(int i = 0; i < sample_no; i++) {
  ParallelFor(0, sample_no, thread_no, [&](size_t i, size_t threadId) {
    std::priority_queue<std::pair<float, hnswlib::labeltype>> result =
        appr_alg->searchKnn(X.colptr(i), kNN + 1);

    if (result.size() != (kNN + 1)) {
      stdout_printf(
          "Unable to find %d results. Probably ef (%f) or M (%f) is too "
          "small\n",
          kNN, ef, M);
    }

    for (size_t j = 0; j <= kNN; j++) {
      auto &result_tuple = result.top();
      dist(i, kNN - j) = result_tuple.first;
      idx(i, kNN - j) = result_tuple.second;

      result.pop();
    }
  });
  stdout_printf("done\n");

  delete (appr_alg);
  delete (space);

  dist = clamp(dist, 0.0, 1.0);
  idx = clamp(idx, 0, sample_no - 1);

  stdout_printf("\tConstructing adaptive-nearest neighbor graph ... ");
  mat Delta;
  mat beta = LC * dist;
  vec beta_sum = zeros(sample_no);
  vec beta_sq_sum = zeros(sample_no);

  mat lambda = zeros(size(beta));
  // lambda.col(0) = datum::inf*ones(sample_no);
  // lambda.col(1) = beta.col(1) + 1;
  int k;
  for (k = 1; k <= kNN; k++) {
    beta_sum += beta.col(k);
    beta_sq_sum += square(beta.col(k));

    lambda.col(k) = (1.0 / (double)k) *
                    (beta_sum + sqrt(k + square(beta_sum) - k * beta_sq_sum));
  }
  lambda.replace(datum::nan, 0);

  lambda = trans(lambda);
  vec node_lambda = zeros(sample_no);
  beta = trans(beta);

  Delta = lambda - beta;
  Delta.shed_row(0);

  sp_mat G(sample_no, sample_no);
  // for(int v = 0; v < sample_no; v++) {
  ParallelFor(0, sample_no, 1, [&](size_t v, size_t threadId) {
    vec delta = Delta.col(v);

    // uvec rows = find(delta > 0, 1, "last");
    uvec rows = find(delta < 0, 1, "first");
    int neighbor_no = rows.n_elem == 0 ? kNN : (rows(0));

    int dst = v;
    rowvec v_dist = dist.row(v);
    rowvec v_idx = idx.row(v);
    for (int i = 1; i < neighbor_no; i++) {
      int src = v_idx(i);
      G(src, dst) = 1.0 - v_dist(i);
    }
  });
  stdout_printf("done\n");

  stdout_printf("\tFinalizing network ... ");
  G.replace(datum::nan, 0);  // replace each NaN with 0

  sp_mat Gt = trans(G);

  sp_mat G_sym;
  if (mutual_edges_only == false) {
    G_sym = (G + Gt);
    G_sym.for_each([](sp_mat::elem_type &val) { val /= 2.0; });
  } else {  // Default to MNN
    stdout_printf("\n\t\tKeeping mutual nearest-neighbors only ... ");
    G_sym = sqrt(G % Gt);
  }
  stdout_printf("done\n");

  G_sym.diag().zeros();

  return (G_sym);
}

sp_mat build_ACTIONet_JS_KstarNN_v2(mat H_stacked, double density = 1.0,
                                    int thread_no = 0, double M = 16,
                                    double ef_construction = 200,
                                    double ef = 10,
                                    bool mutual_edges_only = true) {
  double LC = 1.0 / density;

  if (thread_no <= 0) {
    thread_no = SYS_THREADS_DEF;  // std::thread::hardware_concurrency();
  }

  stdout_printf("Building adaptive network (%d threads):\n", thread_no);
  stdout_printf("\tParameters: density = %0.2f, mutual_edges_only = %s\n",
                density, mutual_edges_only ? "TRUE" : "FALSE");
  FLUSH;

  H_stacked = clamp(H_stacked, 0, 1);
  H_stacked =
      normalise(H_stacked, 1, 0);  // make the norm (sum) of each column 1

  double kappa = 5.0;
  int sample_no = H_stacked.n_cols;
  // int kNN = min(sample_no-1, (int)(kappa*round(sqrt(sample_no)))); // start
  // with uniform k=sqrt(N) ["Pattern Classification" book by Duda et al.]

  int dim = H_stacked.n_rows;
  int max_elements = H_stacked.n_cols;
  hnswlib::JSDSpace *space = new hnswlib::JSDSpace(dim);
  hnswlib::HierarchicalNSW<float> *appr_alg =
      new hnswlib::HierarchicalNSW<float>(space, max_elements, M,
                                           ef_construction);
  appr_alg->setEf(ef);

  stdout_printf("\tBuilding index ... ");
  fmat X = conv_to<fmat>::from(H_stacked); 
   
  ParallelFor(0, max_elements, thread_no, [&](size_t j, size_t threadId) {
    appr_alg->addPoint(X.colptr(j), static_cast<size_t>(j));
  });
  stdout_printf("done\n");

  stdout_printf("\tConstructing k*-NN ... ");

  vector<vector<int>> ii(thread_no);
  vector<vector<int>> jj(thread_no);
  vector<vector<float>> vv(thread_no);

  //		for(int i = 0; i < sample_no; i++) {
  ParallelFor(0, sample_no, thread_no, [&](size_t i, size_t threadId) {
    std::priority_queue<std::pair<float, hnswlib::labeltype>> results =
        appr_alg->searchKStarnn(X.colptr(i), LC);

    while (!results.empty()) {
      auto &res = results.top();
      int j = res.second;
      float v = 1.0 - res.first;
      ii[threadId].push_back(i);
      jj[threadId].push_back(j);
      vv[threadId].push_back(v);

      results.pop();
    }
  });

  vec values;
  umat locations;
  for (int threadId = 0; threadId < thread_no; threadId++) {
    if (threadId == 0) {
      values = conv_to<vec>::from(vv[threadId]);

      uvec iv = conv_to<uvec>::from(ii[threadId]);
      uvec jv = conv_to<uvec>::from(jj[threadId]);
      locations = trans(join_rows(iv, jv));
    } else {
      values = join_vert(values, conv_to<vec>::from(vv[threadId]));

      uvec iv = conv_to<uvec>::from(ii[threadId]);
      uvec jv = conv_to<uvec>::from(jj[threadId]);
      locations = join_rows(locations, trans(join_rows(iv, jv)));
    }
  }
  sp_mat G(locations, values, sample_no, sample_no);

  stdout_printf("done\n");

  delete (appr_alg);
  delete (space);

  stdout_printf("\tFinalizing network ... ");
  G.replace(datum::nan, 0);  // replace each NaN with 0

  sp_mat Gt = trans(G);

  sp_mat G_sym;
  if (mutual_edges_only == false) {
    G_sym = (G + Gt);
    G_sym.for_each([](sp_mat::elem_type &val) { val /= 2.0; });
  } else {  // Default to MNN
    // stdout_printf("\n\t\tKeeping mutual nearest-neighbors only ... ");
    G_sym = sqrt(G % Gt);
  }
  stdout_printf("done\n");

  G_sym.diag().zeros();

  return (G_sym);
}

sp_mat build_ACTIONet_JS_KNN(mat H_stacked, double k, int thread_no = 0,
                             double M = 16, double ef_construction = 200,
                             double ef = 10, bool mutual_edges_only = true) {
  stdout_printf("Building fixed-degree network (k = %d)\n", (int)k);
  if (thread_no <= 0) {
    thread_no = SYS_THREADS_DEF;  // std::thread::hardware_concurrency();
  }

  H_stacked = clamp(H_stacked, 0, 1);
  // H_stacked = normalise(H_stacked, 1, 0); // make the norm (sum) of each
  // column 1

  double kappa = 5.0;
  int sample_no = H_stacked.n_cols;
  int kNN = min(
      sample_no - 1,
      (int)(kappa *
            round(sqrt(sample_no))));  // start with uniform k=sqrt(N) ["Pattern
                                       // Classification" book by Duda et al.]

  int dim = H_stacked.n_rows;
  int max_elements = H_stacked.n_cols;
  hnswlib::JSDSpace *space = new hnswlib::JSDSpace(dim);
  hnswlib::HierarchicalNSW<float> *appr_alg =
      new hnswlib::HierarchicalNSW<float>(space, max_elements, M,
                                           ef_construction);
  appr_alg->setEf(ef);

  stdout_printf("\tBuilding index ... ");
  fmat X = conv_to<fmat>::from(H_stacked);    
  ParallelFor(0, max_elements, thread_no, [&](size_t j, size_t threadId) {
    appr_alg->addPoint(X.colptr(j), static_cast<size_t>(j));
  });
  stdout_printf("done\n");

  sp_mat G(sample_no, sample_no);

  stdout_printf("\tConstructing k*-NN ... ");
  mat idx = zeros(sample_no, kNN + 1);
  mat dist = zeros(sample_no, kNN + 1);
  //		for(int i = 0; i < sample_no; i++) {
  ParallelFor(0, sample_no, thread_no, [&](size_t i, size_t threadId) {
    std::priority_queue<std::pair<float, hnswlib::labeltype>> result =
        appr_alg->searchKnn(X.colptr(i), k);

    if (result.size() != (kNN + 1)) {
      stderr_printf(
          "Unable to find %d results. Probably ef (%f) or M (%f) is too "
          "small\n",
          kNN, ef, M);
    }

    for (size_t j = 0; j < result.size(); j++) {
      auto &result_tuple = result.top();
      G(i, idx(i, kNN - j)) = 1.0 - dist(i, kNN - j);
      result.pop();
    }
  });
  stdout_printf("done\n");

  delete (appr_alg);
  delete (space);

  stdout_printf("\tFinalizing network ... ");
  G.replace(datum::nan, 0);  // replace each NaN with 0

  sp_mat Gt = trans(G);

  sp_mat G_sym;
  if (mutual_edges_only == false) {
    G_sym = (G + Gt);
    G_sym.for_each([](sp_mat::elem_type &val) { val /= 2.0; });
  } else {  // Default to MNN
    stdout_printf("\n\t\tKeeping mutual nearest-neighbors only ... ");
    G_sym = sqrt(G % Gt);
  }
  stdout_printf("done\n");

  G_sym.diag().zeros();

  return (G_sym);
}

sp_mat build_ACTIONet(mat H_stacked, double density = 1.0, int thread_no = 0,
                      double M = 16, double ef_construction = 200,
                      double ef = 10, bool mutual_edges_only = true) {
  if (thread_no <= 0) {
    thread_no = SYS_THREADS_DEF;  // std::thread::hardware_concurrency();
  }

  sp_mat G = build_ACTIONet_JS_KstarNN_v2(
      H_stacked, density, thread_no, M, ef_construction, ef, mutual_edges_only);

  return (G);
}
}  // namespace ACTIONet
