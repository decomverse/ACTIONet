#include <ACTIONet.h>

#include <atomic>
#include <thread>

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

namespace ACTIONet {
arma::vec diffusion_solve_FISTA(arma::sp_mat &adj_mat, arma::vec &prob_dist,
                                double alpha, double rho, double epsilon,
                                int max_iter);

mat PR_linsys(sp_mat &G, sp_mat &X, double alpha = 0.85, int thread_no = -1) {
  X = normalise(X, 1, 0);

  /*
  rowvec d = sum(G, 0);
  uvec idx = find(c != 0);
  d[idx] = 1 / d[idx];

  sp_mat D;
  D.diag() = d;

  sp_mat I = speye(size(G));
  */
  sp_mat P = normalise(G, 1, 0);
  sp_mat I = speye(P.n_rows, P.n_cols);
  sp_mat A = I - alpha * P;
  // mat PR = (1-alpha)*spsolve(A, mat(X), "superlu");
  mat PR = (1 - alpha) * spsolve(A, mat(X));

  return (PR);
}

mat compute_network_diffusion(sp_mat &G, sp_mat &X0, int thread_no = 4,
                              double alpha = 0.85, int max_it = 3) {
  thread_no = std::min(thread_no, (int)X0.n_cols);

  int N = G.n_rows;
  vec z = ones(N);
  vec c = vec(trans(sum(G, 0)));
  uvec idx = find(c);
  z(idx) = ones(idx.n_elem) * (1.0 - alpha);
  z = z / N;

  sp_mat P = alpha * normalise(G, 1, 0);
  X0 = normalise(X0, 1, 0);
  mat X = mat(X0);

  X0 *= N;
  rowvec zt = trans(z);

  for (int it = 0; it < max_it; it++) {
    ParallelFor(0, X.n_cols, thread_no, [&](size_t i, size_t threadId) {
      X.col(i) = P * X.col(i) + X0.col(i) * (zt * X.col(i));
    });
  }
  // X = normalise(X, 1)

  return (X);
}

mat compute_network_diffusion_direct(sp_mat &G, sp_mat &X0, int thread_no = 4,
                                     double alpha = 0.85) {
  // A common struct that cholmod always needs
  cholmod_common c;
  cholmod_start(&c);

  int *Ti, *Tj;
  double *Tx;

  // Construct A
  // Construct transition matrix
  vec d = vec(trans(sum(G, 0)));
  uvec zero_idx = find(d == 0);
  d(zero_idx).ones();
  sp_mat P = G;

  for (sp_mat::iterator it = P.begin(); it != P.end(); ++it) {
    (*it) = (*it) / (sqrt(d(it.row()) * d(it.col())));
  }
  P = -alpha * P;
  P.diag().ones();

  printf("Creating A\n");
  fflush(stdout);
  cholmod_triplet *T = cholmod_allocate_triplet(P.n_rows, P.n_cols, P.n_nonzero,
                                                0, CHOLMOD_REAL, &c);
  T->nnz = P.n_nonzero;
  Ti = static_cast<int *>(T->i);
  Tj = static_cast<int *>(T->j);
  Tx = static_cast<double *>(T->x);
  int idx = 0;
  for (sp_mat::const_iterator it = P.begin(); it != P.end(); ++it) {
    Ti[idx] = it.row();
    Tj[idx] = it.col();
    Tx[idx] = (*it);
    idx++;
  }
  cholmod_sparse *A = cholmod_triplet_to_sparse(T, P.n_nonzero, &c);
  cholmod_free_triplet(&T, &c);

  // Construct B
  printf("Creating B\n");
  fflush(stdout);
  vec d_X = vec(trans(sum(X0, 0)));
  zero_idx = find(d_X == 0);
  d_X(zero_idx).ones();
  sp_mat D_X(X0.n_cols, X0.n_cols);
  D_X.diag() = d_X;

  X0 = normalise(X0, 1, 0);

  T = cholmod_allocate_triplet(X0.n_rows, X0.n_cols, X0.n_nonzero, 0,
                               CHOLMOD_REAL, &c);
  T->nnz = X0.n_nonzero;
  Ti = static_cast<int *>(T->i);
  Tj = static_cast<int *>(T->j);
  Tx = static_cast<double *>(T->x);
  idx = 0;
  for (sp_mat::const_iterator it = X0.begin(); it != X0.end(); ++it) {
    Ti[idx] = it.row();
    Tj[idx] = it.col();
    Tx[idx] = (*it);
    idx++;
  }
  cholmod_sparse *B = cholmod_triplet_to_sparse(T, X0.n_nonzero, &c);
  cholmod_free_triplet(&T, &c);

  // Solve linear system
  printf("Chlmod analyze\n");
  fflush(stdout);
  cholmod_factor *L = cholmod_analyze(A, &c);
  printf("Chlmod factor\n");
  fflush(stdout);
  cholmod_factorize(A, L, &c);
  printf("Solve\n");
  fflush(stdout);
  cholmod_sparse *A_inv_B = cholmod_spsolve(CHOLMOD_A, L, B, &c);

  // Export results
  printf("Export\n");
  fflush(stdout);
  T = cholmod_sparse_to_triplet(A_inv_B, &c);
  Ti = (int *)T->i;
  Tj = (int *)T->j;
  Tx = (double *)T->x;
  umat locations(2, T->nnz);
  for (int k = 0; k < T->nnz; k++) {
    locations(0, k) = Ti[k];
    locations(1, k) = Tj[k];
  }
  mat PR = mat(
      sp_mat(locations, (1 - alpha) * vec(Tx, T->nnz), X0.n_rows, X0.n_cols));

  PR = normalise(PR, 1, 0);
  PR = PR * D_X;

  // Free up matrices
  cholmod_free_factor(&L, &c);
  cholmod_free_sparse(&A, &c);
  cholmod_free_sparse(&B, &c);
  cholmod_free_triplet(&T, &c);
  cholmod_finish(&c);

  return (PR);
}

mat compute_network_diffusion_fast(sp_mat &G, sp_mat &X0,
  int thread_no = 4, double alpha = 0.85, int max_it = 5) {
  thread_no = std::min(thread_no, (int)X0.n_cols);

  int n = G.n_rows;

  cholmod_common chol_c;
  cholmod_start(&chol_c);

  int *Ti, *Tj;
  double *Tx;

  sp_mat P = alpha * normalise(G, 1, 0);

  cholmod_triplet *T = cholmod_allocate_triplet(P.n_rows, P.n_cols, P.n_nonzero,
                                                0, CHOLMOD_REAL, &chol_c);
  T->nnz = P.n_nonzero;
  Ti = static_cast<int *>(T->i);
  Tj = static_cast<int *>(T->j);
  Tx = static_cast<double *>(T->x);
  int idx = 0;
  for (sp_mat::const_iterator it = P.begin(); it != P.end(); ++it) {
    Ti[idx] = it.row();
    Tj[idx] = it.col();
    Tx[idx] = (*it);
    idx++;
  }
  cholmod_sparse *AS = cholmod_triplet_to_sparse(T, P.n_nonzero, &chol_c);
  cholmod_free_triplet(&T, &chol_c);

  vec z = ones(n);
  vec cs = vec(trans(sum(G, 0)));
  uvec nnz_idx = find(cs > 0);
  z(nnz_idx) = ones(nnz_idx.n_elem) * (1.0 - alpha);
  z = z / n;

  X0 = normalise(X0, 1, 0);
  mat X = mat(X0);
  X0 *= n;
  rowvec zt = trans(z);

  for (int it = 0; it < max_it; it++) {
    mat Y = X;
    ParallelFor(0, X.n_cols, thread_no, [&](size_t i, size_t threadId) {
      dsdmult('n', n, n, AS, X.colptr(i), Y.colptr(i), &chol_c);
      X.col(i) = Y.col(i) + X0.col(i) * (zt * X.col(i));
    });
  }

  // Free up matrices
  cholmod_free_triplet(&T, &chol_c);
  cholmod_free_sparse(&AS, &chol_c);
  cholmod_finish(&chol_c);

  return (X);
}

}  // namespace ACTIONet
