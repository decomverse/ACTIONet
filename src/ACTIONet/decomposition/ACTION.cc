#include "ACTIONet.h"
#include "cholmod.h"

namespace ACTIONet {

// Solves the standard Archetypal Analysis (AA) problem
field<mat> run_AA(mat &A, mat &W0, int max_it = 100, double min_delta = 1e-6) {
  int sample_no = A.n_cols;
  int d = A.n_rows;   // input dimension
  int k = W0.n_cols;  // AA components

  mat C = zeros(sample_no, k);
  mat H = zeros(k, sample_no);

  mat W = W0;
  vec c(sample_no);

  double old_RSS = 0;
  // printf("(New) %d- %d\n", k, max_it);

  for (int it = 0; it < max_it; it++) {
    mat C_old = C;
    mat H_old = H;

	double A_norm = norm(A, "fro");
    H = run_simplex_regression(W, A, true);
    mat R = A - W * H;
    mat Ht = trans(H);
    for (int i = 0; i < k; i++) {
      vec w = W.col(i);
      vec h = Ht.col(i);

      double norm_sq = arma::dot(h, h);
      if (norm_sq < double(10e-8)) {
        // singular
        int max_res_idx = index_max(rowvec(sum(square(R), 0)));
        W.col(i) = A.col(max_res_idx);
        c.zeros();
        c(max_res_idx) = 1;
        C.col(i) = c;
      } else {
        // b = (1.0 / norm_sq) *R*ht + w;
        vec b = w;
        cblas_dgemv(CblasColMajor, CblasNoTrans, R.n_rows, R.n_cols,
                    (1.0 / norm_sq), R.memptr(), R.n_rows, h.memptr(), 1, 1,
                    b.memptr(), 1);

        C.col(i) = run_simplex_regression(A, b, false);

        vec w_new = A * C.col(i);
        vec delta = (w - w_new);

        // Rank-1 update: R += delta*h
        cblas_dger(CblasColMajor, R.n_rows, R.n_cols, 1.0, delta.memptr(), 1,
                   h.memptr(), 1, R.memptr(), R.n_rows);

        W.col(i) = w_new;
      }
    }
    double RSS = norm(R, "fro");
	double delta_RSS = abs(RSS - old_RSS)/A_norm;
    old_RSS = RSS;
    /*
    double delta_C = norm(C - C_old, "fro") / norm(C, "fro");
    double delta_H = norm(H - H_old, "fro") / norm(H, "fro");	
    printf("\t<%d, %d>- norm_RSS = %e, delta_RSS = %e, delta_C = %.3e, delta_H = %.3e\n", k, it, RSS/A_norm, delta_RSS, delta_C, delta_H);
    printf("\t<%d, %d>- norm_RSS = %e\n", k, it, delta_RSS);
    */
    
    if(delta_RSS < min_delta)
            break;
  }

  C = clamp(C, 0, 1);
  C = normalise(C, 1);
  H = clamp(H, 0, 1);
  H = normalise(H, 1);

  field<mat> decomposition(2, 1);
  decomposition(0) = C;
  decomposition(1) = H;

  return decomposition;
}

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

// Solves separable NMF problem
SPA_results run_SPA_rows_sparse(sp_mat &A, int k) {
  int m = A.n_rows;
  int n = A.n_cols;

  printf("Computing square values ... ");
  fflush(stdout);
  sp_mat A_sq = square(A);
  printf("done\n");

  cholmod_common chol_c;
  cholmod_start(&chol_c);
  chol_c.final_ll = 1; /* LL' form of simplicial factorization */

  cholmod_sparse_struct *AS = new cholmod_sparse_struct;
  as_cholmod_sparse(AS, A);

  cholmod_sparse_struct *AS_sq = new cholmod_sparse_struct;
  as_cholmod_sparse(AS_sq, A_sq);

  SPA_results res;

  uvec K(k);  // selected columns from A

  printf("Computing norms ... ");
  vec o = ones(n);
  vec normM(m);
  dsdmult('n', m, n, AS_sq, o.memptr(), normM.memptr(), &chol_c);
  vec normM1 = normM;
  printf("done\n");

  mat U(n, k);

  vec norm_trace = zeros(k);
  double eps = 1e-6;
  for (int i = 0; i < k; i++) {
    // Find the column with maximum norm. In case of having more than one column
    // with almost very small diff in norm, pick the one that originally had the
    // largest norm
    double a = max(normM);
    norm_trace(i) = a;

    uvec b = find((a * ones(m, 1) - normM) / a <= eps);

    if (b.n_elem > 1) {
      uword idx = index_max(normM1(b));
      K(i) = b(idx);
    } else {
      K(i) = b(0);
    }

    // Pick row
    U.col(i) = vec(trans(A.row(K(i))));

    // Orthogonalize with respect to current basis
    for (int j = 0; j < i - 1; j++) {
      U.col(i) = U.col(i) - dot(U.col(j), U.col(i)) * U.col(j);
    }
    U.col(i) = U.col(i) / norm(U.col(i), 2);

    // Update column norms
    vec u = U.col(i);
    for (int j = i - 1; 0 <= j; j--) {
      u = u - dot(U.col(j), u) * U.col(j);
    }
    vec r(m);
    dsdmult('n', m, n, AS, u.memptr(), r.memptr(), &chol_c);

    uvec idx = find(U > 0);
    double perc = 100 * idx.n_elem / U.n_elem;
    stdout_printf("\t%d- res_norm = %f, U_density = %.2f%% (%d nnz)\n", i, a, perc,
           idx.n_elem);

    normM = normM - (r % r);
  }

  res.selected_columns = K;
  res.column_norms = norm_trace;

  delete[] AS->x;
  delete[] AS->i;
  delete[] AS->p;
  delete AS;

  delete[] AS_sq->x;
  delete[] AS_sq->i;
  delete[] AS_sq->p;
  delete AS_sq;

  cholmod_finish(&chol_c);

  return res;
}

// Solves separable NMF problem
SPA_results run_SPA(mat &A, int k) {
  SPA_results res;

  int n = A.n_cols;
  uvec K(k);  // selected columns from A
  K.zeros();

  rowvec normM = sum(square(A), 0);
  rowvec normM1 = normM;

  mat U(A.n_rows, k);

  vec norm_trace = zeros(k);
  double eps = 1e-16;

  for (int i = 1; i <= k; i++) {
    // Find the column with maximum norm. In case of having more than one column
    // with almost very small diff in norm, pick the one that originally had the
    // largest norm
    double a = max(normM);
    norm_trace(i - 1) = a;

    uvec b = find((a * ones(1, n) - normM) / a <= eps);
    if (b.n_elem == 0) {
      break;
    } else if (b.n_elem > 1) {
      uword idx = index_max(normM1(b));
      K(i - 1) = b(idx);
    } else {
      K(i - 1) = b(0);
    }

    // Pick column
    U.col(i - 1) = A.col(K(i - 1));

    // Orthogonalize with respect to current basis
    if (i > 1) {
      for (int j = 1; j <= i - 1; j++) {
        U.col(i - 1) =
            U.col(i - 1) - sum(U.col(j - 1) % U.col(i - 1)) * U.col(j - 1);
      }
    }
    double nm = norm(U.col(i - 1), 2);
    if (nm > 0) U.col(i - 1) /= nm;

    // Update column norms
    vec u = U.col(i - 1);
    if (i > 1) {
      for (int j = i - 1; 1 <= j; j--) {
        u = u - sum(U.col(j - 1) % u) * U.col(j - 1);
      }
    }
    normM = normM - square(u.t() * A);
    normM.transform([](double val) { return (val < 0 ? 0 : val); });
  }

  res.selected_columns = K;
  res.column_norms = norm_trace;

  return res;
}

ACTION_results run_ACTION(mat &S_r, int k_min, int k_max, int thread_no,
                          int max_it = 100, double min_delta = 1e-6) {
  if (thread_no <= 0) {
    thread_no = SYS_THREADS_DEF;
  }

  int feature_no = S_r.n_rows;

  stdout_printf("Running ACTION (%d threads):", thread_no);
  FLUSH;

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
  // int total = k_min-1;
  char status_msg[50];

  sprintf(status_msg, "Iterating from k = %d ... %d:", k_min, k_max);
  stderr_printf("\n\t%s %d/%d finished", status_msg, current_k,
                (k_max - k_min + 1));
  FLUSH;

  ParallelFor(k_min, k_max + 1, thread_no, [&](size_t kk, size_t threadId) {
    SPA_results SPA_res = run_SPA(X_r, kk);
    trace.selected_cols[kk] = SPA_res.selected_columns;

    mat W = X_r.cols(trace.selected_cols[kk]);

    field<mat> AA_res;

    AA_res = run_AA(X_r, W, max_it, min_delta);
    // AA_res = run_AA_old(X_r, W);
    trace.C[kk] = AA_res(0);
    trace.H[kk] = AA_res(1);
    current_k++;

    stderr_printf("\r\t%s %d/%d finished", status_msg, current_k,
                  (k_max - k_min + 1));
    FLUSH;
  });
  stdout_printf("\r\t%s %d/%d finished\n", status_msg, current_k,
                (k_max - k_min + 1));

  return trace;
}

ACTION_results run_ACTION_plus(mat &S_r, int k_min, int k_max, int max_it = 100,
                               double min_delta = 1e-6, int max_trial = 3) {
  stdout_printf("Running ACTION++\n");

  int D = std::min((int)S_r.n_rows, (int)S_r.n_cols);
  if (k_max == -1) k_max = D;

  k_min = std::max(k_min, 2);
  k_max = std::min(k_max, D);

  ACTION_results trace;

  trace.H = field<mat>(k_max + 1, 1);
  trace.C = field<mat>(k_max + 1, 1);
  trace.selected_cols = field<uvec>(k_max + 1, 1);

  mat X_r = normalise(S_r, 1);  // ATTENTION!
  SPA_results SPA_res = run_SPA(X_r, D);
  uvec selected_cols = SPA_res.selected_columns;

  mat W = mat(X_r.col(selected_cols(0)));

  field<mat> AA_res;
  int cur_idx = 0, jj, kk;
  stdout_printf("Iterating from k=%d ... %d (max trial = %d)\n", k_min, k_max,
         max_trial);
  for (kk = k_min; kk <= k_max; kk++) {
    stdout_printf("\tk = %d\n", kk);

    for (jj = 0; jj < max_trial; jj++) {
      cur_idx++;
      stdout_printf("\t\tTrial %d: candidate %d = %d ... ", jj + 1, cur_idx + 1,
             selected_cols(cur_idx));
      mat W_tmp = join_rows(W, X_r.col(selected_cols(cur_idx)));

      AA_res = run_AA(X_r, W_tmp, max_it, min_delta);

      vec influential_cells = vec(trans(sum(spones(sp_mat(AA_res(0))), 0)));
      int trivial_counts = (int)sum(influential_cells <= 1);

      if ((trivial_counts == 0)) {
        stdout_printf("success\n");
        selected_cols(kk - 1) = selected_cols(cur_idx);
        break;
      }

      stdout_printf("failed\n");
      if ((cur_idx == (D - 1))) {
        stdout_printf("Reached end of the line!\n");
        break;
      }
    }

    if ((jj == max_trial) || (cur_idx == (D - 1))) {
      break;
    }

    trace.C[kk] = AA_res(0);
    trace.H[kk] = AA_res(1);
    trace.selected_cols(kk) = selected_cols(span(0, kk - 1));

    W = X_r * AA_res(0);
  }

  trace.C = trace.C.rows(0, kk - 1);
  trace.H = trace.H.rows(0, kk - 1);
  trace.selected_cols = trace.selected_cols.rows(0, kk - 1);

  return trace;
}

field<mat> run_AA_with_prior(mat &A, mat &W0, mat &W_prior, int max_it = 50,
                             double min_delta = 1e-16) {
  int sample_no = A.n_cols;
  int d = A.n_rows;   // input dimension
  int k = W0.n_cols;  // AA components

  mat C = zeros(sample_no, k);
  mat H = zeros(k, sample_no);

  mat W = W0;
  vec c(sample_no);

  // printf("(New) %d- %d\n", k, max_it);

  for (int it = 0; it < max_it; it++) {
    mat combined_W = join_rows(W, W_prior);
    mat combined_H = run_simplex_regression(combined_W, A, true);

    H = combined_H.rows(span(0, k - 1));

    // mat C_old = C;
    mat R = A - W * H;
    mat Ht = trans(H);
    for (int i = 0; i < k; i++) {
      vec w = W.col(i);
      vec h = Ht.col(i);

      double norm_sq = arma::dot(h, h);
      if (norm_sq < double(10e-8)) {
        // singular
        int max_res_idx = index_max(rowvec(sum(square(R), 0)));
        W.col(i) = A.col(max_res_idx);
        c.zeros();
        c(max_res_idx) = 1;
        C.col(i) = c;
      } else {
        // b = (1.0 / norm_sq) *R*ht + w;
        vec b = w;
        cblas_dgemv(CblasColMajor, CblasNoTrans, R.n_rows, R.n_cols,
                    (1.0 / norm_sq), R.memptr(), R.n_rows, h.memptr(), 1, 1,
                    b.memptr(), 1);

        C.col(i) = run_simplex_regression(A, b, false);

        vec w_new = A * C.col(i);
        vec delta = (w - w_new);

        // Rank-1 update: R += delta*h
        cblas_dger(CblasColMajor, R.n_rows, R.n_cols, 1.0, delta.memptr(), 1,
                   h.memptr(), 1, R.memptr(), R.n_rows);

        W.col(i) = w_new;
      }
    }
    /*
    double delta = arma::max(rowvec(sum(abs(C - C_old)))) / 2.0;

    //double RSS = norm(R, "fro"); RSS *= RSS;
    //printf("\t<%d, %d>- RSS = %.3e, delta = %.3e\n", l, it, RSS, delta);

    if(delta < min_delta)
            break;
    */
  }

  C = clamp(C, 0, 1);
  C = normalise(C, 1);
  H = clamp(H, 0, 1);
  H = normalise(H, 1);

  field<mat> decomposition(2, 1);
  decomposition(0) = C;
  decomposition(1) = H;

  return decomposition;
}

ACTION_results run_subACTION(mat &S_r, mat &W_parent, mat &H_parent, int kk,
                             int k_min, int k_max, int thread_no,
                             int max_it = 50, double min_delta = 1e-16) {
  int feature_no = S_r.n_rows;

  stdout_printf("Running subACTION (%d threads) for parent archetype %d\n", thread_no, kk + 1);

  if (k_max == -1) k_max = (int)S_r.n_cols;

  k_min = std::max(k_min, 2);
  k_max = std::min(k_max, (int)S_r.n_cols);

  mat X_r = normalise(S_r, 1);  // ATTENTION!

  vec h = vec(trans(H_parent.row(kk)));
  mat W_prior = W_parent;
  W_prior.shed_col(kk);

  mat X_r_scaled = X_r;  // To deflate or not deflate!

  for (int i = 0; i < X_r_scaled.n_cols; i++) {
    X_r_scaled.col(i) *= h[i];
  }

  ACTION_results trace;
  trace.H = field<mat>(k_max + 1);
  trace.C = field<mat>(k_max + 1);
  trace.selected_cols = field<uvec>(k_max + 1);

  int current_k = 0;
  int total = k_min - 1;
  stdout_printf("Iterating from k=%d ... %d\n", k_min, k_max);
  ParallelFor(k_min, k_max + 1, thread_no, [&](size_t kkk, size_t threadId) {
    total++;
    stdout_printf("\tk = %d\n", total);

    SPA_results SPA_res = run_SPA(X_r_scaled, kkk);
    trace.selected_cols[kkk] = SPA_res.selected_columns;

    mat W = X_r.cols(trace.selected_cols[kkk]);

    field<mat> AA_res;

    AA_res = run_AA_with_prior(X_r_scaled, W, W_prior, max_it, min_delta);

    trace.C[kkk] = AA_res(0);
    trace.H[kkk] = AA_res(1);
  });

  return trace;
}

}  // namespace ACTIONet
