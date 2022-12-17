#include <ACTIONet.h>

#include <atomic>
#include <thread>

#include <cholmod.h>

sp_mat &as_arma_sparse(cholmod_sparse *chol_A, sp_mat &A,
                       cholmod_common *chol_c);

void dsdmult(char transpose, int m, int n, void *a, double *b, double *c,
             cholmod_common *chol_cp);

cholmod_sparse *as_cholmod_sparse(sp_mat &A, cholmod_sparse *chol_A,
                                  cholmod_common *chol_c);

arma::vec diffusion_solve_FISTA(arma::sp_mat &adj_mat, arma::vec &prob_dist,
                                double alpha, double rho, double epsilon,
                                int max_iter);

namespace ACTIONet
{
  mat PR_linsys(sp_mat &G, sp_mat &X, double alpha = 0.85, int thread_no = -1)
  {
    X = normalise(X, 1, 0);

    sp_mat P = normalise(G, 1, 0);
    sp_mat I = speye(P.n_rows, P.n_cols);
    sp_mat A = I - alpha * P;
    // mat PR = (1-alpha)*spsolve(A, mat(X), "superlu");
    mat PR = (1 - alpha) * spsolve(A, mat(X));

    return (PR);
  }

  mat compute_network_diffusion(sp_mat &G, sp_mat &X0, int thread_no = 4,
                                double alpha = 0.85, int max_it = 3)
  {
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

    for (int it = 0; it < max_it; it++)
    {

      parallelFor(
          0, X.n_cols, [&](size_t i)
          { X.col(i) = P * X.col(i) + X0.col(i) * (zt * X.col(i)); },
          thread_no);
    }

    return (X);
  }

  mat compute_network_diffusion_direct(sp_mat &G, sp_mat &X0, int thread_no = 4,
                                       double alpha = 0.85)
  {
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

    for (sp_mat::iterator it = P.begin(); it != P.end(); ++it)
    {
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
    for (sp_mat::const_iterator it = P.begin(); it != P.end(); ++it)
    {
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
    for (sp_mat::const_iterator it = X0.begin(); it != X0.end(); ++it)
    {
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
    for (int k = 0; k < T->nnz; k++)
    {
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

  mat compute_network_diffusion_fast(sp_mat &G, sp_mat &X0, int thread_no = 4,
                                     double alpha = 0.85, int max_it = 5)
  {
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
    for (sp_mat::const_iterator it = P.begin(); it != P.end(); ++it)
    {
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

    for (int it = 0; it < max_it; it++)
    {
      mat Y = X;
      parallelFor(
          0, X.n_cols, [&](size_t i)
          {
                    dsdmult('n', n, n, AS, X.colptr(i), Y.colptr(i), &chol_c);
                    X.col(i) = Y.col(i) + X0.col(i) * (zt * X.col(i)); },
          thread_no);
    }

    // Free up matrices
    cholmod_free_triplet(&T, &chol_c);
    cholmod_free_sparse(&AS, &chol_c);
    cholmod_finish(&chol_c);

    return (X);
  }

  sp_mat normalize_adj(sp_mat &G, int norm_type)
  {
    vec row_sums = zeros(G.n_rows);
    vec col_sums = zeros(G.n_cols);

    sp_mat::iterator it = G.begin();
    sp_mat::iterator it_end = G.end();
    for (; it != it_end; ++it)
    {
      col_sums[it.col()] += (*it);
      row_sums[it.row()] += (*it);
    }
    uvec idxr = find(row_sums == 0);
    uvec idxc = find(col_sums == 0);

    row_sums.transform([](double val)
                       { return (val == 0 ? 1 : val); });
    col_sums.transform([](double val)
                       { return (val == 0 ? 1 : val); });

    // Update
    sp_mat P = G;
    if (norm_type == 0) // Column-normalize
    {
      for (it = P.begin(); it != P.end(); ++it)
      {
        double w = col_sums[it.col()];
        (*it) /= w;
      }
      for (int k = 0; k < idxc.n_elem; k++)
      {
        int j = idxc(k);
        P(j, j) = 1.0;
      }
    }
    else if (norm_type == 1) // Row-normalize
    {
      for (it = P.begin(); it != P.end(); ++it)
      {
        double w = row_sums[it.row()];
        (*it) /= w;
      }
      for (int k = 0; k < idxr.n_elem; k++)
      {
        int i = idxr(k);
        P(i, i) = 1.0;
      }
    }
    else if (norm_type == 2)
    {
      for (it = P.begin(); it != P.end(); ++it)
      {
        double w = sqrt(row_sums[it.row()] * col_sums[it.col()]);
        (*it) /= w;
      }
    }

    return (P);
  }
  // P is already an stochastic matrix
  mat compute_network_diffusion_Chebyshev(sp_mat &P, mat &X0, int thread_no,
                                          double alpha, int max_it, double res_threshold)
  {
    if (alpha == 1)
    {
      fprintf(stderr, "alpha should be in (0, 1). Value of %.2f was provided.\n", alpha);
      return (X0);
    }
    else if (alpha == 0)
    {
      return (X0);
    }

    alpha = 1 - alpha; // Traditional defitition is to have alpha as weight of prior. Here, alpha is depth of difffusion

    mat mPPreviousScore = X0; // zeros(size(X0));
    mat mPreviousScore = (1 - alpha) * spmat_mat_product_parallel(P, mPPreviousScore, thread_no) + alpha * X0;
    double mu = 0.0, muPPrevious = 1.0, muPrevious = 1 / (1 - alpha);

    if (max_it <= 0)
      return (mPreviousScore);

    mat mScore;
    for (int i = 0; i < max_it; i++)
    {
      double mu = 2.0 / (1.0 - alpha) * muPrevious - muPPrevious;

      mScore = 2 * (muPrevious / mu) * spmat_mat_product_parallel(P, mPreviousScore, thread_no) - (muPPrevious / mu) * mPPreviousScore + (2 * muPrevious) / ((1 - alpha) * mu) * alpha * X0;

      double res = norm(mScore - mPreviousScore);
      if (res < res_threshold)
      {
        break;
      }

      // Change variables
      muPPrevious = muPrevious;
      muPrevious = mu;
      mPPreviousScore = mPreviousScore;
      mPreviousScore = mScore;
    }

    // Temporarty fix. Sometimes diffusion values become small negative numbers
    double m0 = min(min(X0));
    if (0 <= m0)
    {
      mScore = clamp(mScore, 0, max(max(mScore)));
    }

    return (mScore);
  }

} // namespace ACTIONet
