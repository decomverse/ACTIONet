#include "ACTIONet.h"

namespace ACTIONet {

template <class Function>
inline void ParallelFor(size_t start, size_t end, size_t numThreads,
                        Function fn) {
  if (numThreads <= 0) {
    numThreads = SYS_THREADS_DEF;
  }

  if (numThreads == 1) {
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

    for (size_t threadId = 0; threadId < numThreads; ++threadId) {
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

// Solves the weighted Archetypal Analysis (AA) problem
field<mat> run_weighted_AA(mat &A, mat &W0, vec w, int max_it = 50,
                           double min_delta = 0.01) {
  int N = A.n_cols;
  field<mat> decomposition(2);

  if (N != w.n_elem) {
    fprintf(stderr,
            "Number of elements in the weight vector should match the total "
            "number of samples (columns in A)\n");
    return (decomposition);
  }

  w = clamp(w, 0, 1);
  mat A_scaled = A;
  for (int i = 0; i < N; i++) {
    A_scaled.col(i) *= w[i];
  }
  decomposition = run_AA(A_scaled, W0, max_it, min_delta);

  /*
  for(int i = 0; i < N; i++) {
          decomposition(0).row(i) /= w[i];
          decomposition(1).col(i) *= w[i];
  }
  */
  mat C = decomposition(0);
  mat weighted_archs = A_scaled * C;
  mat H = run_simplex_regression(weighted_archs, A, false);
  decomposition(1) = H;

  return (decomposition);
}

ACTION_results run_weighted_ACTION(mat &S_r, vec w, int k_min, int k_max,
                                   int thread_no, int max_it = 50,
                                   double min_delta = 1e-16) {
  int feature_no = S_r.n_rows;

  printf("Running ACTION\n");

  if (k_max == -1) k_max = (int)S_r.n_cols;

  k_min = std::max(k_min, 2);
  k_max = std::min(k_max, (int)S_r.n_cols);

  ACTION_results trace;
  /*
  trace.H.resize(k_max + 1);
  trace.C.resize(k_max + 1);
  trace.selected_cols.resize(k_max + 1);
  */

  trace.H = field<mat>(k_max + 1);
  trace.C = field<mat>(k_max + 1);
  trace.selected_cols = field<uvec>(k_max + 1);

  int N = S_r.n_cols;
  if (N != w.n_elem) {
    fprintf(stderr,
            "Number of elements in the weight vector should match the total "
            "number of samples (columns in S_r)\n");
    return (trace);
  }

  w = clamp(w, 0, 1);
  mat X_r = normalise(S_r, 1);  // ATTENTION!

  mat X_r_scaled = X_r;
  for (int i = 0; i < N; i++) {
    X_r_scaled.col(i) *= w[i];
  }

  int current_k = 0;
  int total = k_min - 1;
  printf("Iterating from k=%d ... %d\n", k_min, k_max);
  ParallelFor(k_min, k_max + 1, thread_no, [&](size_t kk, size_t threadId) {
    total++;
    printf("\tk = %d\n", total);
    SPA_results SPA_res = run_SPA(X_r_scaled, kk);
    trace.selected_cols[kk] = SPA_res.selected_columns;

    mat W = X_r_scaled.cols(trace.selected_cols[kk]);

    field<mat> AA_res = run_AA(X_r_scaled, W, max_it, min_delta);

    mat C = AA_res(0);
    mat weighted_archs = X_r_scaled * C;
    mat H = run_simplex_regression(weighted_archs, X_r, false);
    AA_res(1) = H;

    trace.C[kk] = AA_res(0);
    trace.H[kk] = AA_res(1);
  });

  return trace;
}

ACTION_results run_ACTION_dev(mat &S_r, int k_min, int k_max, int thread_no,
                              bool auto_stop = true, int max_it = 30,
                              double min_delta = 0.01) {
  int feature_no = S_r.n_rows;

  printf("Running ACTION (developmental version)\n");

  if (k_max == -1) k_max = (int)S_r.n_cols;

  k_min = std::max(k_min, 2);
  k_max = std::min(k_max, (int)S_r.n_cols);

  ACTION_results trace;
  /*
  trace.H.resize(k_max + 1);
  trace.C.resize(k_max + 1);
  trace.selected_cols.resize(k_max + 1);
  */

  trace.H = field<mat>(k_max + 1);
  trace.C = field<mat>(k_max + 1);
  trace.selected_cols = field<uvec>(k_max + 1);

  mat X_r = normalise(S_r, 1);  // ATTENTION!

  int current_k = 0;
  printf("Iterating from k=%d ... %d (auto stop = %d)\n", k_min, k_max,
         auto_stop);
  for (int kk = k_min; kk <= k_max; kk++) {
    printf("\tk = %d\n", kk);
    SPA_results SPA_res = run_SPA(X_r, kk);
    trace.selected_cols[kk] = SPA_res.selected_columns;

    mat W = X_r.cols(trace.selected_cols[kk]);
    if (kk > k_min) {
      W.cols(span(0, kk - 2)) = X_r * trace.C[kk - 1];
    }
    // W.print("W");

    field<mat> AA_res = run_AA(X_r, W, max_it, min_delta);

    if (auto_stop) {
      mat C = AA_res(0);
      bool has_trivial_arch = false;
      for (int c = 0; c < C.n_cols; c++) {
        int r = 0, nnz_counts = 0;
        while (r < C.n_rows) {
          if (C(r, c) > 0) nnz_counts++;

          if (nnz_counts > 1) break;

          r++;
        }
        if (nnz_counts <= 1) {
          has_trivial_arch = true;
          break;
        }
      }
      if (has_trivial_arch > 0) {
        printf("\t\tFound trivial archetypes at k = %d\n", kk);

        break;
      }
    }

    current_k = std::max(current_k, kk);

    trace.C[kk] = AA_res(0);
    trace.H[kk] = AA_res(1);
  }
  /*
  trace.H.resize(current_k+1);
  trace.C.resize(current_k+1);
  trace.selected_cols.resize(current_k+1);
  */

  return trace;
}

field<mat> Online_update_AA(mat &Xt, mat &D, mat &A, mat &B) {
  // Compute archetype coefficients using the last learned dictionary
  mat Ct = run_simplex_regression(D, Xt, false);

  // Just in case!
  Ct = clamp(Ct, 0, 1);
  Ct = normalise(Ct, 1);
  mat Ct_T = trans(Ct);

  // Update sufficient statistics
  mat At = A + Ct * Ct_T;
  mat Bt = B + Xt * Ct_T;

  // Update the dictionary using block-coordinate-descent (BCD)
  mat Dt(size(D));
  for (int j = 0; j < D.n_cols; j++) {
    vec u = D.col(j) + (1.0 / At(j, j)) * (Bt.col(j) - D * At.col(j));
    Dt.col(j) = u / std::max(norm(u, 2), 1.0);
  }

  field<mat> decomposition(4, 1);

  decomposition(0) = At;
  decomposition(1) = Bt;
  decomposition(2) = Ct;
  decomposition(3) = Dt;

  return (decomposition);
}

field<mat> run_online_AA(mat &X, mat &D0, field<uvec> samples) {
  int m = X.n_rows;
  int n = X.n_cols;
  int k = D0.n_cols;

  mat At = zeros(k, k);
  mat Bt = zeros(m, k);
  mat Ct(k, n), Ct_T(n, k);

  // Just in case
  X = normalise(X, 2, 0);
  mat Dt = normalise(D0, 2, 0);

  for (int t = 0; t < samples.n_elem; t++) {
    // Extract the next batch
    uvec idx = samples(t);
    mat Xt = X.cols(idx);

    // Compute archetype coefficients using the last learned dictionary
    Ct = run_simplex_regression(Dt, Xt, false);
    Ct_T = trans(Ct);

    // Update sufficient statistics
    At += Ct * Ct_T;
    Bt += Xt * Ct_T;

    // Update the dictionary using block-coordinate-descent (BCD)
    for (int j = 0; j < k; j++) {
      vec u = Dt.col(j) + (1.0 / At(j, j)) * (Bt.col(j) - Dt * At.col(j));
      Dt.col(j) = u / std::max(norm(u, 2), 1.0);
    }
  }

  // Just in case!
  Ct = clamp(Ct, 0, 1);
  Ct = normalise(Ct, 1);

  field<mat> decomposition(4, 1);

  decomposition(0) = At;
  decomposition(1) = Bt;
  decomposition(2) = Ct;
  decomposition(3) = Dt;

  return (decomposition);
}

Online_ACTION_results run_online_ACTION(mat &S_r, field<uvec> samples,
                                        int k_min, int k_max, int thread_no) {
  int feature_no = S_r.n_rows;

  printf("Running Online ACTION\n");

  if (k_max == -1) k_max = (int)S_r.n_cols;

  k_min = std::max(k_min, 2);
  k_max = std::min(k_max, (int)S_r.n_cols);

  Online_ACTION_results trace;

  trace.A = field<mat>(k_max + 1);
  trace.B = field<mat>(k_max + 1);
  trace.C = field<mat>(k_max + 1);
  trace.D = field<mat>(k_max + 1);
  trace.selected_cols = field<uvec>(k_max + 1);

  mat X_r_L1 = normalise(S_r, 1, 0);
  mat X_r_L2 = normalise(S_r, 2, 0);

  int current_k = 0;
  int total = k_min - 1;
  printf("Iterating from k=%d ... %d\n", k_min, k_max);
  ParallelFor(k_min, k_max + 1, thread_no, [&](size_t kk, size_t threadId) {
    total++;
    printf("\tk = %d\n", total);
    SPA_results SPA_res = run_SPA(X_r_L1, kk);
    trace.selected_cols[kk] = SPA_res.selected_columns;

    mat W = X_r_L2.cols(trace.selected_cols[kk]);

    field<mat> AA_res;
    AA_res = run_online_AA(X_r_L2, W, samples);

    trace.A[kk] = AA_res(0);
    trace.B[kk] = AA_res(1);
    trace.C[kk] = AA_res(2);
    trace.D[kk] = AA_res(3);
  });

  return trace;
}

}  // namespace ACTIONet
