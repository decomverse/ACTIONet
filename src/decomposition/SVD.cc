#include "ACTIONet.h"

#include <cholmod.h>
sp_mat &as_arma_sparse(cholmod_sparse *chol_A, sp_mat &A,
                       cholmod_common *chol_c);

void dsdmult(char transpose, int m, int n, void *a, double *b, double *c,
             cholmod_common *chol_cp);

cholmod_sparse *as_cholmod_sparse(sp_mat &A, cholmod_sparse *chol_A,
                                  cholmod_common *chol_c);

// #define StdNorm(v, n, engine) randN_Marsaglia(v, n, engine)

inline void StdNorm(double *v, int n, std::mt19937_64 engine)
{
  for (int ii = 0; ii < n - 1; ii += 2)
  {
    auto paired = aarand::standard_normal(engine);
    v[ii] = paired.first;
    v[ii + 1] = paired.second;
  }
  auto paired = aarand::standard_normal(engine);
  v[n - 1] = paired.first;
}

namespace ACTIONet
{
  field<mat> orient_SVD(field<mat> SVD_res)
  {
    mat U = SVD_res(0);
    vec s = SVD_res(1);
    mat V = SVD_res(2);

    int dim = s.n_elem;
    uvec mask_idx;

    for (int i = 0; i < dim; i++)
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
        U.col(i) *= -1;
        V.col(i) *= -1;
      }
    }

    field<mat> out(3);
    out(0) = U;
    out(1) = s;
    out(2) = V;

    return (out);
  }

  field<mat> IRLB_SVD(sp_mat &A, int dim, int iters = 1000, int seed = 0,
                      int verbose = 1)
  {
    int m = A.n_rows;
    int n = A.n_cols;

    dim = min(dim, min(m, n) - 1);

    if (verbose)
    {
      stdout_printf("IRLB (sparse) -- A: %d x %d\n", A.n_rows, A.n_cols);
      FLUSH;
    }

    cholmod_common chol_c;
    cholmod_start(&chol_c);
    chol_c.final_ll = 1; /* LL' form of simplicial factorization */

    cholmod_sparse *AS = as_cholmod_sparse(A, AS, &chol_c);

    double eps = 3e-13;
    double tol = 1e-05, svtol = 1e-5;

    int work = dim + 7;
    int lwork = 7 * work * (1 + work);

    double *s = new double[dim];
    double *U = new double[m * work];
    double *V = new double[n * work];

    double *V1 = new double[n * work];
    double *U1 = new double[m * work];
    double *W = new double[m * work];
    double *F = new double[n];
    double *B = new double[work * work];
    double *BU = new double[work * work];
    double *BV = new double[work * work];
    double *BS = new double[work];
    double *BW = new double[lwork];
    double *res = new double[work];
    double *T = new double[lwork];
    double *svratio = new double[work];

    mat tmp(B, work, work, false);
    mat BUmat(BU, work, work, false);
    vec BSvec(BS, work, false);
    mat BVmat(BV, work, work, false);

    double d, S, R, R_F, SS;
    double *x;
    int jj, kk;
    int converged;
    int info, j, k = 0;
    int iter = 0;
    double Smax = 0;

    memset(B, 0, work * work * sizeof(double));
    memset(svratio, 0, work * sizeof(double));

    double alpha = 1, beta = 0;
    int inc = 1;

    vec v, y;

    // Initialize first column of V
    // pcg32 engine(seed);
    std::mt19937_64 engine(seed);

    double ss;
    StdNorm(V, n, engine);

    ss = 0;
    for (int i = 0; i < n; i++)
      ss += V[i];

    /* Main iteration */
    while (iter < iters)
    {
      j = 0;

      /*  Normalize starting vector */
      if (iter == 0)
      {
        d = cblas_dnrm2(n, V, inc);
        d = 1 / d;
        cblas_dscal(n, d, V, inc);
      }
      else
        j = k;

      // Compute Ax
      x = V + j * n;

      dsdmult('n', m, n, AS, x, W + j * m, &chol_c);

      if (iter > 0)
        orthog(W, W + j * m, T, m, j, 1);

      S = cblas_dnrm2(m, W + j * m, inc);
      SS = 1.0 / S;
      cblas_dscal(m, SS, W + j * m, inc);

      /* The Lanczos process */
      while (j < work)
      {
        dsdmult('t', m, n, AS, W + j * m, F, &chol_c);

        SS = -S;
        cblas_daxpy(n, SS, V + j * n, inc, F, inc);
        orthog(V, F, T, n, j + 1, 1);

        if (j + 1 < work)
        {
          R_F = cblas_dnrm2(n, F, inc);
          R = 1.0 / R_F;

          if (R_F < eps)
          { // near invariant subspace

            StdNorm(F, n, engine);
            ss = 0;
            for (int i = 0; i < n; i++)
              ss += F[i];

            orthog(V, F, T, n, j + 1, 1);
            R_F = cblas_dnrm2(n, F, inc);
            R = 1.0 / R_F;
            R_F = 0;
          }

          memmove(V + (j + 1) * n, F, n * sizeof(double));
          cblas_dscal(n, R, V + (j + 1) * n, inc);
          B[j * work + j] = S;
          B[(j + 1) * work + j] = R_F;

          x = V + (j + 1) * n;

          dsdmult('n', m, n, AS, x, W + (j + 1) * m, &chol_c);

          /* One step of classical Gram-Schmidt */
          R = -R_F;
          cblas_daxpy(m, R, W + j * m, inc, W + (j + 1) * m, inc);

          /* full re-orthogonalization of W_{j+1} */
          orthog(W, W + (j + 1) * m, T, m, j + 1, 1);
          S = cblas_dnrm2(m, W + (j + 1) * m, inc);
          SS = 1.0 / S;

          if (S < eps)
          {
            StdNorm(W + (j + 1) * m, m, engine);
            ss = 0;
            for (int i = 0; i < n; i++)
              ss += W[(j + 1) * m + i];

            orthog(W, W + (j + 1) * m, T, m, j + 1, 1);
            S = cblas_dnrm2(m, W + (j + 1) * m, inc);
            SS = 1.0 / S;
            cblas_dscal(m, SS, W + (j + 1) * m, inc);
            S = 0;
          }
          else
            cblas_dscal(m, SS, W + (j + 1) * m, inc);
        }
        else
        {
          B[j * work + j] = S;
        }

        j++;
      }

      arma::svd(BUmat, BSvec, BVmat, tmp, "dc");
      BVmat = trans(BVmat);

      R_F = cblas_dnrm2(n, F, inc);
      R = 1.0 / R_F;
      cblas_dscal(n, R, F, inc);

      /* Force termination after encountering linear dependence */
      if (R_F < eps)
        R_F = 0;

      Smax = 0;
      for (jj = 0; jj < j; ++jj)
      {
        if (BS[jj] > Smax)
          Smax = BS[jj];
        svratio[jj] = fabs(svratio[jj] - BS[jj]) / BS[jj];
      }

      for (kk = 0; kk < j; ++kk)
        res[kk] = R_F * BU[kk * work + (j - 1)];

      /* Update k to be the number of converged singular values. */
      convtests(j, dim, tol, svtol, Smax, svratio, res, &k, &converged, S);
      if (converged == 1)
      {
        iter++;
        break;
      }

      for (jj = 0; jj < j; ++jj)
        svratio[jj] = BS[jj];

      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, k, j, alpha, V, n,
                  BV, work, beta, V1, n);

      memmove(V, V1, n * k * sizeof(double));
      memmove(V + n * k, F, n * sizeof(double));

      memset(B, 0, work * work * sizeof(double));
      for (jj = 0; jj < k; ++jj)
      {
        B[jj * work + jj] = BS[jj];
        B[k * work + jj] = res[jj];
      }

      /*   Update the left approximate singular vectors */
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, k, j, alpha, W, m,
                  BU, work, beta, U1, m);

      memmove(W, U1, m * k * sizeof(double));
      iter++;
    }

    /* Results */
    memmove(s, BS, dim * sizeof(double)); /* Singular values */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, dim, work, alpha, W,
                m, BU, work, beta, U, m);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, dim, work, alpha, V,
                n, BV, work, beta, V1, n);
    memmove(V, V1, n * dim * sizeof(double));

    field<mat> out(3);
    out(0) = mat(U, m, dim);
    out(1) = vec(s, dim);
    out(2) = mat(V, n, dim);

    delete[] s;
    delete[] U;
    delete[] V;
    delete[] V1;
    delete[] U1;
    delete[] W;
    delete[] F;
    delete[] B;
    delete[] BU;
    delete[] BV;
    delete[] BS;
    delete[] BW;
    delete[] res;
    delete[] T;
    delete[] svratio;

    cholmod_free_sparse(&AS, &chol_c);
    cholmod_finish(&chol_c);

    if (converged != 1)
    {
      stderr_printf(
          "IRLB did NOT converge! Try increasing the number of iterations\n");
      FLUSH;
    }

    return (orient_SVD(out));
  }

  field<mat> IRLB_SVD(mat &A, int dim, int iters = 1000, int seed = 0,
                      int verbose = 1)
  {
    double eps = 3e-13;
    double tol = 1e-05, svtol = 1e-5;

    int m = A.n_rows;
    int n = A.n_cols;

    dim = min(dim, min(m, n) - 1);

    int work = dim + 7;
    int lwork = 7 * work * (1 + work);

    if (verbose)
    {
      stdout_printf("IRLB (dense) -- A: %d x %d\n", A.n_rows, A.n_cols);
      FLUSH;
    }

    double *s = new double[dim];
    double *U = new double[m * work];
    double *V = new double[n * work];

    double *V1 = new double[n * work];
    double *U1 = new double[m * work];
    double *W = new double[m * work];
    double *F = new double[n];
    double *B = new double[work * work];
    double *BU = new double[work * work];
    double *BV = new double[work * work];
    double *BS = new double[work];
    double *BW = new double[lwork];
    double *res = new double[work];
    double *T = new double[lwork];
    double *svratio = new double[work];

    mat tmp(B, work, work, false);
    mat BUmat(BU, work, work, false);
    vec BSvec(BS, work, false);
    mat BVmat(BV, work, work, false);

    double d, S, R, R_F, SS;
    double *x;
    int jj, kk;
    int converged;
    int info, j, k = 0;
    int iter = 0;
    double Smax = 0;

    memset(B, 0, work * work * sizeof(double));
    memset(svratio, 0, work * sizeof(double));

    double alpha = 1, beta = 0;
    int inc = 1;

    vec v, y;

    // Initialize first column of V
    // pcg32 engine(seed);
    std::mt19937_64 engine(seed);

    StdNorm(V, n, engine);

    /* Main iteration */
    while (iter < iters)
    {
      j = 0;

      /*  Normalize starting vector */
      if (iter == 0)
      {
        d = cblas_dnrm2(n, V, inc);
        d = 1 / d;
        cblas_dscal(n, d, V, inc);
      }
      else
        j = k;

      // Compute Ax
      x = V + j * n;
      v = vec(x, n, true);
      y = A * v;
      memcpy(W + j * m, y.memptr(), y.n_elem * sizeof(double));

      if (iter > 0)
        orthog(W, W + j * m, T, m, j, 1);

      S = cblas_dnrm2(m, W + j * m, inc);
      SS = 1.0 / S;
      cblas_dscal(m, SS, W + j * m, inc);

      /* The Lanczos process */
      while (j < work)
      {
        v = vec(W + j * m, m, true);
        y = trans(trans(v) * A);
        memcpy(F, y.memptr(), y.n_elem * sizeof(double));

        SS = -S;
        cblas_daxpy(n, SS, V + j * n, inc, F, inc);
        orthog(V, F, T, n, j + 1, 1);

        if (j + 1 < work)
        {
          R_F = cblas_dnrm2(n, F, inc);
          R = 1.0 / R_F;

          if (R_F < eps)
          { // near invariant subspace
            StdNorm(F, n, engine);

            orthog(V, F, T, n, j + 1, 1);
            R_F = cblas_dnrm2(n, F, inc);
            R = 1.0 / R_F;
            R_F = 0;
          }

          memmove(V + (j + 1) * n, F, n * sizeof(double));
          cblas_dscal(n, R, V + (j + 1) * n, inc);
          B[j * work + j] = S;
          B[(j + 1) * work + j] = R_F;

          x = V + (j + 1) * n;
          v = vec(x, n, true);
          y = A * v;
          memcpy(W + (j + 1) * m, y.memptr(), y.n_elem * sizeof(double));

          /* One step of classical Gram-Schmidt */
          R = -R_F;
          cblas_daxpy(m, R, W + j * m, inc, W + (j + 1) * m, inc);

          /* full re-orthogonalization of W_{j+1} */
          orthog(W, W + (j + 1) * m, T, m, j + 1, 1);
          S = cblas_dnrm2(m, W + (j + 1) * m, inc);
          SS = 1.0 / S;

          if (S < eps)
          {
            StdNorm(W + (j + 1) * m, m, engine);

            orthog(W, W + (j + 1) * m, T, m, j + 1, 1);
            S = cblas_dnrm2(m, W + (j + 1) * m, inc);
            SS = 1.0 / S;
            cblas_dscal(m, SS, W + (j + 1) * m, inc);
            S = 0;
          }
          else
            cblas_dscal(m, SS, W + (j + 1) * m, inc);
        }
        else
        {
          B[j * work + j] = S;
        }

        j++;
      }

      arma::svd(BUmat, BSvec, BVmat, tmp, "dc");
      BVmat = trans(BVmat);

      R_F = cblas_dnrm2(n, F, inc);
      R = 1.0 / R_F;
      cblas_dscal(n, R, F, inc);

      /* Force termination after encountering linear dependence */
      if (R_F < eps)
        R_F = 0;

      Smax = 0;
      for (jj = 0; jj < j; ++jj)
      {
        if (BS[jj] > Smax)
          Smax = BS[jj];
        svratio[jj] = fabs(svratio[jj] - BS[jj]) / BS[jj];
      }

      for (kk = 0; kk < j; ++kk)
        res[kk] = R_F * BU[kk * work + (j - 1)];

      /* Update k to be the number of converged singular values. */
      convtests(j, dim, tol, svtol, Smax, svratio, res, &k, &converged, S);
      if (converged == 1)
      {
        iter++;
        break;
      }

      for (jj = 0; jj < j; ++jj)
        svratio[jj] = BS[jj];

      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, k, j, alpha, V, n,
                  BV, work, beta, V1, n);

      memmove(V, V1, n * k * sizeof(double));
      memmove(V + n * k, F, n * sizeof(double));

      memset(B, 0, work * work * sizeof(double));
      for (jj = 0; jj < k; ++jj)
      {
        B[jj * work + jj] = BS[jj];
        B[k * work + jj] = res[jj];
      }

      /*   Update the left approximate singular vectors */
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, k, j, alpha, W, m,
                  BU, work, beta, U1, m);

      memmove(W, U1, m * k * sizeof(double));
      iter++;
    }

    /* Results */
    memmove(s, BS, dim * sizeof(double)); /* Singular values */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, dim, work, alpha, W,
                m, BU, work, beta, U, m);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, dim, work, alpha, V,
                n, BV, work, beta, V1, n);
    memmove(V, V1, n * dim * sizeof(double));

    field<mat> out(3);
    out(0) = mat(U, m, dim);
    out(1) = vec(s, dim);
    out(2) = mat(V, n, dim);

    delete[] s;
    delete[] U;
    delete[] V;
    delete[] V1;
    delete[] U1;
    delete[] W;
    delete[] F;
    delete[] B;
    delete[] BU;
    delete[] BV;
    delete[] BS;
    delete[] BW;
    delete[] res;
    delete[] T;
    delete[] svratio;

    if (converged != 1)
    {
      stderr_printf(
          "IRLB did NOT converge! Try increasing the number of iterations\n");
    }

    return (orient_SVD(out));
  }

  //****************************************************************************************************************************************************************************
  // From: Xu Feng, Yuyang Xie, and Yaohang Li, "Fast Randomzied SVD for Sparse
  // Data," in Proc. the 10th Asian Conference on Machine Learning (ACML),
  // Beijing, China, Nov. 2018.
  //****************************************************************************************************************************************************************************
  field<mat> FengSVD(sp_mat &A, int dim, int iters, int seed = 0,
                     int verbose = 1)
  {
    int s = 5;

    int m = A.n_rows;
    int n = A.n_cols;

    dim = min(dim, min(m, n) + s - 1);

    if (verbose)
    {
      stdout_printf("Feng (sparse) -- A: %d x %d\n", A.n_rows, A.n_cols);
      FLUSH;
    }

    vec S;
    mat Q, L, U, V;
    field<mat> SVD_out;

    if (m < n)
    {
      randNorm(n, dim + s, seed);
      Q = A * Q;
      if (iters == 0)
      {
        SVD_out = eigSVD(Q);
        Q = SVD_out(0);
      }
      else
      {
        lu(L, U, Q);
        Q = L;
      }

      for (int i = 1; i <= iters; i++)
      {
        if (verbose)
        {
          stderr_printf("\r\tIteration %d/%d", i, iters);
          FLUSH;
        }

        if (i == iters)
        {
          SVD_out = eigSVD(A * (trans(A) * Q));
          Q = SVD_out(0);
        }
        else
        {
          lu(L, U, A * (trans(A) * Q));
          Q = L;
        }
      }

      SVD_out = eigSVD(trans(A) * Q);
      V = SVD_out(0);
      S = vec(SVD_out(1));
      U = SVD_out(2);

      U = Q * fliplr(U.cols(s, dim + s - 1));
      V = fliplr(V.cols(s, dim + s - 1));
      S = flipud(S(span(s, dim + s - 1)));
    }
    else
    {
      Q = randNorm(m, dim + s, seed);
      Q = trans(A) * Q;
      if (iters == 0)
      {
        SVD_out = eigSVD(Q);
        Q = SVD_out(0);
      }
      else
      {
        lu(L, U, Q);
        Q = L;
      }

      for (int i = 1; i <= iters; i++)
      {
        if (verbose)
        {
          stderr_printf("\r\tIteration %d/%d", i, iters);
          FLUSH;
        }

        if (i == iters)
        {
          SVD_out = eigSVD(trans(A) * (A * Q));
          Q = SVD_out(0);
        }
        else
        {
          lu(L, U, trans(A) * (A * Q));
          Q = L;
        }
      }
      if (verbose)
      {
        stdout_printf("\r\tIteration %d/%d", iters, iters);
        FLUSH;
      }

      SVD_out = eigSVD(A * Q);
      U = SVD_out(0);
      S = vec(SVD_out(1));
      V = SVD_out(2);

      U = fliplr(U.cols(s, dim + s - 1));
      V = Q * fliplr(V.cols(s, dim + s - 1));
      S = flipud(S(span(s, dim + s - 1)));
    }

    field<mat> out(3);
    out(0) = U;
    out(1) = S;
    out(2) = V;

    return (orient_SVD(out));
  }

  // From: Xu Feng, Yuyang Xie, and Yaohang Li, "Fast Randomzisped SVD for
  // Sparse Data," in Proc. the 10th Asian Conference on Machine Learning
  // (ACML), Beijing, China, Nov. 2018.
  field<mat> FengSVD(mat &A, int dim, int iters, int seed = 0, int verbose = 1)
  {
    int s = 5;

    int m = A.n_rows;
    int n = A.n_cols;

    dim = min(dim, min(m, n) + s - 1);

    if (verbose)
    {
      stdout_printf("Feng (dense) -- A: %d x %d\n", A.n_rows, A.n_cols);
      FLUSH;
    }

    vec S;
    mat Q, L, U, V;
    field<mat> SVD_out;

    if (m < n)
    {
      Q = randNorm(n, dim + s, seed);

      Q = A * Q;
      if (iters == 0)
      {
        SVD_out = eigSVD(Q);
        Q = SVD_out(0);
      }
      else
      {
        lu(L, U, Q);
        Q = L;
      }

      for (int i = 1; i <= iters; i++)
      {
        if (verbose)
        {
          stderr_printf("\r\tIteration %d/%d", i, iters);
          FLUSH;
        }
        if (i == iters)
        {
          SVD_out = eigSVD(A * (trans(A) * Q));
          Q = SVD_out(0);
        }
        else
        {
          lu(L, U, A * (trans(A) * Q));
          Q = L;
        }
        if (verbose)
          stdout_printf("done\n");
      }
      if (verbose)
      {
        stdout_printf("\r\tIteration %d/%d", iters, iters);
        FLUSH;
      }

      SVD_out = eigSVD(trans(A) * Q);
      V = SVD_out(0);
      S = vec(SVD_out(1));
      U = SVD_out(2);

      U = Q * fliplr(U.cols(s, dim + s - 1));
      V = fliplr(V.cols(s, dim + s - 1));
      S = flipud(S(span(s, dim + s - 1)));
    }
    else
    {
      Q = randNorm(m, dim + s, seed);
      Q = trans(A) * Q;
      if (iters == 0)
      {
        SVD_out = eigSVD(Q);
        Q = SVD_out(0);
      }
      else
      {
        lu(L, U, Q);
        Q = L;
      }

      for (int i = 1; i <= iters; i++)
      {
        if (verbose)
        {
          stderr_printf("\r\tIteration %d/%d", i, iters);
          FLUSH;
        }
        if (i == iters)
        {
          SVD_out = eigSVD(trans(A) * (A * Q));
          Q = SVD_out(0);
        }
        else
        {
          lu(L, U, trans(A) * (A * Q));
          Q = L;
        }
      }
      if (verbose)
      {
        stdout_printf("\r\tIteration %d/%d", iters, iters);
        FLUSH;
      }

      SVD_out = eigSVD(A * Q);
      U = SVD_out(0);
      S = vec(SVD_out(1));
      V = SVD_out(2);

      U = fliplr(U.cols(s, dim + s - 1));
      V = Q * fliplr(V.cols(s, dim + s - 1));
      S = flipud(S(span(s, dim + s - 1)));
    }

    field<mat> out(3);
    out(0) = U;
    out(1) = S;
    out(2) = V;

    return (orient_SVD(out));
  }

  //**************************************************************************************************************************************************************************************************
  // From: N Halko, P. G Martinsson, and J. A Tropp. Finding structure with
  // randomness: Probabilistic algorithms for constructing approximate matrix
  // decompositions. Siam Review, 53(2):217-288, 2011.
  //**************************************************************************************************************************************************************************************************
  field<mat> HalkoSVD(sp_mat &A, int dim, int iters, int seed = 0,
                      int verbose = 1)
  {
    field<mat> results(3);

    int m = A.n_rows;
    int n = A.n_cols;

    dim = min(dim, min(m, n) - 2);

    int l = dim + 2;

    vec s;
    mat R, Q;
    mat U, V, X;

    if (verbose)
    {
      stdout_printf("Halko (sparse) -- A: %d x %d\n", A.n_rows, A.n_cols);
      FLUSH;
    }

    if (m < n)
    {
      R = randNorm(l, m, seed);

      sp_mat At = A.t();
      Q = At * R.t();
    }
    else
    {
      R = randNorm(n, l, seed);

      Q = A * R;
    }

    // Form a matrix Q whose columns constitute a well-conditioned basis for the
    // columns of the earlier Q.
    gram_schmidt(Q);

    if (m < n)
    {
      // Conduct normalized power iterations.
      for (int it = 1; it <= iters; it++)
      {
        if (verbose)
        {
          stderr_printf("\r\tIteration %d/%d", it, iters);
          FLUSH;
        }

        Q = A * Q;
        gram_schmidt(Q);

        Q = A.t() * Q;
        gram_schmidt(Q);
      }
      if (verbose)
      {
        stdout_printf("\r\tIteration %d/%d", iters, iters);
        FLUSH;
      }

      X = mat(A * Q);
      svd_econ(U, s, V, X);
      V = Q * V;
    }
    else
    {
      // Conduct normalized power iterations.
      for (int it = 1; it <= iters; it++)
      {
        if (verbose)
        {
          stderr_printf("\r\tIteration %d/%d", it, iters);
          FLUSH;
        }

        Q = A.t() * Q;
        gram_schmidt(Q);

        Q = A * Q; // Apply A to a random matrix, obtaining Q.
        gram_schmidt(Q);
      }
      if (verbose)
      {
        stdout_printf("\r\tIteration %d/%d", iters, iters);
        FLUSH;
      }

      // SVD Q' applied to the centered A to obtain approximations to the
      // singular values and right singular vectors of the A;

      X = mat(Q.t() * A);
      svd_econ(U, s, V, X);
      U = Q * U;
    }

    U.shed_cols(dim, dim + 1);
    s = s(span(0, dim - 1));
    V.shed_cols(dim, dim + 1);

    results(0) = U;
    results(1) = s;
    results(2) = V;

    return (orient_SVD(results));
  }

  field<mat> HalkoSVD(mat &A, int dim, int iters, int seed = 0, int verbose = 1)
  {
    field<mat> results(3);

    int m = A.n_rows;
    int n = A.n_cols;

    dim = min(dim, min(m, n) - 2);

    int l = dim + 2;

    vec s;
    mat R, Q;
    mat U, V, X;

    if (verbose)
    {
      stdout_printf("Halko (dense) -- A: %d x %d\n", A.n_rows, A.n_cols);
      FLUSH;
    }

    if (m < n)
    {
      R = randNorm(l, m, seed);

      mat At = A.t();
      Q = At * R.t();
    }
    else
    {
      R = randNorm(n, l, seed);

      Q = A * R;
    }

    // Form a matrix Q whose columns constitute a well-conditioned basis for the
    // columns of the earlier Q.
    gram_schmidt(Q);

    if (m < n)
    {
      // Conduct normalized power iterations.=
      for (int it = 1; it <= iters; it++)
      {
        // stdout_printf("\tIteration %d\n", it);
        if (verbose)
        {
          stderr_printf("\r\tIteration %d/%d", it, iters);
          FLUSH;
        }

        Q = A * Q;
        gram_schmidt(Q);

        Q = A.t() * Q;
        gram_schmidt(Q);
      }
      if (verbose)
      {
        stdout_printf("\r\tIteration %d/%d\n", iters, iters);
        FLUSH;
      }

      X = mat(A * Q);
      svd_econ(U, s, V, X);
      V = Q * V;
    }
    else
    {
      // Conduct normalized power iterations.
      for (int it = 1; it <= iters; it++)
      {
        if (verbose)
        {
          stderr_printf("\r\tIteration %d/%d", it, iters);
          FLUSH;
        }

        Q = A.t() * Q;
        gram_schmidt(Q);

        Q = A * Q; // Apply A to a random matrix, obtaining Q.
        gram_schmidt(Q);
      }
      if (verbose)
      {
        stdout_printf("\r\tIteration %d/%d\n", iters, iters);
        FLUSH;
      }
      // SVD Q' applied to the centered A to obtain approximations to the
      // singular values and right singular vectors of the A;

      X = mat(Q.t() * A);
      svd_econ(U, s, V, X);
      U = Q * U;
    }

    U.shed_cols(dim, dim + 1);
    s = s(span(0, dim - 1));
    V.shed_cols(dim, dim + 1);

    results(0) = U;
    results(1) = s;
    results(2) = V;

    return (orient_SVD(results));
  }
} // namespace ACTIONet
