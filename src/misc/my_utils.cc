#include "my_utils.h"
#include "ACTIONet.h"

#include <cmath>
#include <cstdlib>
#include <limits>
#include <cholmod.h>

sp_mat &as_arma_sparse(cholmod_sparse *chol_A, sp_mat &A,
                       cholmod_common *chol_c);

void dsdmult(char transpose, int m, int n, void *a, double *b, double *c,
             cholmod_common *chol_cp);

cholmod_sparse *as_cholmod_sparse(sp_mat &A, cholmod_sparse *chol_A,
                                  cholmod_common *chol_c);

double r8_normal_01_cdf_inverse(double p);

#define A1 (-3.969683028665376e+01)
#define A2 2.209460984245205e+02
#define A3 (-2.759285104469687e+02)
#define A4 1.383577518672690e+02
#define A5 (-3.066479806614716e+01)
#define A6 2.506628277459239e+00

#define B1 (-5.447609879822406e+01)
#define B2 1.615858368580409e+02
#define B3 (-1.556989798598866e+02)
#define B4 6.680131188771972e+01
#define B5 (-1.328068155288572e+01)

#define C1 (-7.784894002430293e-03)
#define C2 (-3.223964580411365e-01)
#define C3 (-2.400758277161838e+00)
#define C4 (-2.549732539343734e+00)
#define C5 4.374664141464968e+00
#define C6 2.938163982698783e+00

#define D1 7.784695709041462e-03
#define D2 3.224671290700398e-01
#define D3 2.445134137142996e+00
#define D4 3.754408661907416e+00

#define P_LOW 0.02425
/* P_high = 1 - p_low*/
#define P_HIGH 0.97575

namespace ACTIONet
{
  /**
   * @brief Generate a single random number using the capped Tausworthe RNG
   *
   * @details
   * This generates random numbers according to the process described in [1]. As
   * an additional step, the resulting random number is capped to 0xFFFFFFFF
   * using a bitwise and. This is done to yield the range [0, 2^32-1]. On
   * return, the state variables are updated.
   *
   * [1]: @article{l1996maximally,
   *   title={Maximally equidistributed combined Tausworthe generators},
   *   author={L’ecuyer, Pierre},
   *   journal={Mathematics of Computation of the American Mathematical
   *   Society},
   *   volume={65},
   *   number={213},
   *   pages={203--213},
   *   year={1996}
   *   }
   *
   * @param[in,out] state 	pointer to current state array
   *
   * @return a generated random number
   */
  uint32_t lfsr113(uint64_t **state)
  {
    uint64_t z1, z2, z3, z4;
    uint64_t b;

    z1 = (*state)[0];
    z2 = (*state)[1];
    z3 = (*state)[2];
    z4 = (*state)[3];

    b = (((z1 << 6) ^ z1) >> 13);
    z1 = (((z1 & 4294967294) << 18) ^ b);

    b = (((z2 << 2) ^ z2) >> 27);
    z2 = (((z2 & 4294967288) << 2) ^ b);

    b = (((z3 << 13) ^ z3) >> 21);
    z3 = (((z3 & 4294967280) << 7) ^ b);

    b = (((z4 << 3) ^ z4) >> 12);
    z4 = (((z4 & 4294967168) << 13) ^ b);

    b = (z1 ^ z2 ^ z3 ^ z4);

    (*state)[0] = z1;
    (*state)[1] = z2;
    (*state)[2] = z3;
    (*state)[3] = z4;

    b = b & 0xFFFFFFFF;

    return ((uint32_t)b);
  }

  /**
   * @brief Seed the Tausworthe RNG using a seed value
   *
   * @details
   * This function seeds the state array using a supplied seed value. As noted
   * in [1] (see lfsr113()), the values of z1, z2, z3, and z4 should be larger
   * than 1, 7, 15, and 127 respectively.
   *
   * @param[in] seed 	user supplied seed value for the RNG
   * @param[out] state  	state of the RNG
   */
  void lfsr113_seed(uint32_t seed, uint64_t **state)
  {
    uint64_t z1 = 2, z2 = 8, z3 = 16, z4 = 128;

    z1 = (z1 * (seed + 1));
    z2 = (z2 * (seed + 1));
    z3 = (z3 * (seed + 1));
    z4 = (z4 * (seed + 1));

    z1 = (z1 > 1) ? z1 : z1 + 1;
    z2 = (z2 > 7) ? z2 : z2 + 7;
    z3 = (z3 > 15) ? z3 : z3 + 15;
    z4 = (z4 > 127) ? z4 : z4 + 127;

    (*state)[0] = z1;
    (*state)[1] = z2;
    (*state)[2] = z3;
    (*state)[3] = z4;
  }

  long double normsinv(long double p)
  {
    long double x;
    long double q, r, u, e;
    if ((0 < p) && (p < P_LOW))
    {
      q = sqrt(-2 * log(p));
      x = (((((C1 * q + C2) * q + C3) * q + C4) * q + C5) * q + C6) /
          ((((D1 * q + D2) * q + D3) * q + D4) * q + 1);
    }
    else
    {
      if ((P_LOW <= p) && (p <= P_HIGH))
      {
        q = p - 0.5;
        r = q * q;
        x = (((((A1 * r + A2) * r + A3) * r + A4) * r + A5) * r + A6) * q /
            (((((B1 * r + B2) * r + B3) * r + B4) * r + B5) * r + 1);
      }
      else
      {
        if ((P_HIGH < p) && (p < 1))
        {
          q = sqrt(-2 * log(1 - p));
          x = -(((((C1 * q + C2) * q + C3) * q + C4) * q + C5) * q + C6) /
              ((((D1 * q + D2) * q + D3) * q + D4) * q + 1);
        }
      }
    }

    /* If you are compiling this under UNIX OR LINUX, you may uncomment this block
  for better accuracy. if(( 0 < p)&&(p < 1)){ e = 0.5 * erfc(-x/sqrt(2)) - p; u
  = e * sqrt(2*M_PI) * exp(x*x/2); x = x - u/(1 + x*u/2);
  }
  */

    return x;
  }

  void randN_normsinv(double *values, int n)
  {
    for (int i = 0; i < n; i++)
    {
      long double u = rand() / (double)RAND_MAX;
      long double z = normsinv(u);
      values[i] = (double)z;
    }

    return;
  }

  // Marsaglia algorithm
  void randN_Marsaglia(double *values, int n, pcg32 rng)
  {
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);

    for (int i = 0; i < n; i += 2)
    {
      double x, y, rsq, f;
      do
      {
        x = 2.0 * uniform_dist(rng) - 1.0;
        y = 2.0 * uniform_dist(rng) - 1.0;
        rsq = x * x + y * y;
      } while (rsq >= 1. || rsq == 0.);
      f = sqrt(-2.0 * log(rsq) / rsq);

      values[i] = x * f;
      if (i < (n - 1))
      {
        values[i + 1] = y * f;
      }
    }

    return;
  }

  // Box_Muller Algorithm
  void randN_BM(double *values, int n, uint64_t **state)
  {
    static const double epsilon = std::numeric_limits<double>::min();
    static const double two_pi = 2.0 * 3.14159265358979323846;

    double kk = (1.0 / 4294967288.0);

    double z0, z1;
    for (int i = 0; i < n; i += 2)
    {
      double u1, u2;
      do
      {
        u1 = lfsr113(state) * kk;
        u2 = lfsr113(state) * kk;

      } while (u1 <= epsilon);

      z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
      z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);

      values[i] = z0;
      if (i < (n - 1))
      {
        values[i + 1] = z1;
      }
    }

    return;
  }

  mat sampleUnif(int l, int m, double a, double b, int seed)
  {
    std::default_random_engine gen(seed);
    std::uniform_real_distribution<double> unifDist(a, b);

    mat R(l, m);
    for (int j = 0; j < m; j++)
    {
      for (int i = 0; i < l; i++)
      {
        R(i, j) = unifDist(gen);
      }
    }
    return R;
  }

  mat randNorm(int l, int m, int seed)
  {
    std::default_random_engine gen(seed);
    std::normal_distribution<double> normDist(0.0, 1.0);

    mat R(l, m);
    for (int j = 0; j < m; j++)
    {
      for (int i = 0; i < l; i++)
      {
        R(i, j) = normDist(gen);
      }
    }
    return R;
  }

  void randNorm_inplace(int n, double *out, int seed = 0)
  {
    std::default_random_engine gen(seed);
    std::normal_distribution<double> normDist(0.0, 1.0);

    for (int i = 0; i < n; i++)
    {
      out[i] = normDist(gen);
    }
  }

  void gram_schmidt(mat &A)
  {
    for (uword i = 0; i < A.n_cols; ++i)
    {
      for (uword j = 0; j < i; ++j)
      {
        double r = dot(A.col(i), A.col(j));
        A.col(i) -= r * A.col(j);
      }

      double col_norm = norm(A.col(i), 2);

      if (col_norm < 1E-4)
      {
        for (uword k = i; k < A.n_cols; ++k)
          A.col(k).zeros();

        return;
      }
      A.col(i) /= col_norm;
    }
  }

  field<mat> eigSVD(mat A)
  {
    int n = A.n_cols;
    mat B = trans(A) * A;

    vec d;
    mat V;
    eig_sym(d, V, B);
    d = sqrt(d);

    // Compute U
    sp_mat S(n, n);
    S.diag() = 1 / d;
    mat U = (S * trans(V)) * trans(A);
    U = trans(U);

    field<mat> out(3);

    out(0) = U;
    out(1) = d;
    out(2) = V;

    return (out);
  }

  mat robust_zscore(mat &A, int dim, int thread_no)
  {
    int N = A.n_cols;
    if (dim != 0)
    {
      N = A.n_rows;
    }

    parallelFor(
        0, N,
        [&](size_t j)
        {
          vec v = A.col(j);
          if (dim == 0)
          {
            v = A.col(j);
          }
          else
          {
            v = A.row(j);
          }
          double med = arma::median(v);
          double mad = arma::median(arma::abs(v - med));

          vec z = (v - med) / mad;
          if (dim == 0)
          {
            A.col(j) = z;
          }
          else
          {
            A.row(j) = z;
          }
        },
        thread_no);
    A.replace(datum::nan, 0); // replace each NaN with 0

    return A;
  }

  mat zscore(mat &A, int dim, int thread_no)
  {
    int N = A.n_cols;
    if (dim != 0)
    {
      N = A.n_rows;
    }

    parallelFor(
        0, N,
        [&](size_t j)
        {
          vec v = A.col(j);
          if (dim == 0)
          {
            v = A.col(j);
          }
          else
          {
            v = A.row(j);
          }
          double mu = arma::mean(v);
          double sigma = arma::stddev(v);

          vec z = (v - mu) / sigma;
          if (dim == 0)
          {
            A.col(j) = z;
          }
          else
          {
            A.row(j) = z;
          }
        },
        thread_no);
    A.replace(datum::nan, 0); // replace each NaN with 0

    return A;
  }

  mat tzscoret(mat &A)
  {
    mat At = A.t();
    A = zscore(At);
    return (A.t());
  }

  mat normalize_mat(mat &X, int normalization, int dim)
  {
    mat X_norm = X;
    if (normalization == 1)
    {
      X_norm = normalise(X_norm, 1, dim);
    }
    if (normalization == 2)
    {
      X_norm = normalise(X_norm, 2, dim);
    }
    if (normalization == -1)
    {
      X_norm = zscore(X_norm, dim);
    }

    return (X_norm);
  }

  sp_mat normalize_mat(sp_mat &X, int normalization, int dim)
  {
    sp_mat X_norm = X;
    if (normalization == 1)
    {
      X_norm = normalise(X_norm, 1, dim);
    }
    if (normalization == 2)
    {
      X_norm = normalise(X_norm, 2, dim);
    }
    return (X_norm);
  }

  /* Rank a numeric vector giving ties their average rank */
  vec rank_vec(vec x, int method)
  {
    int n = x.n_elem;
    vec ranks(n);
    uvec indx = sort_index(x, "ascend");

    int ib = 0, i;
    double b = x[indx[0]];
    for (i = 1; i < n; ++i)
    {
      if (x[indx[i]] != b)
      { /* consecutive numbers differ */
        if (ib < i - 1)
        { /* average of sum of ranks */

          double rnk = method == 0 ? ((i - 1 + ib + 2) / 2.0) : (i - 1);
          for (int j = ib; j <= i - 1; ++j)
            ranks[indx[j]] = rnk;
        }
        else
        {
          ranks[indx[ib]] = (double)(ib + 1);
        }
        b = x[indx[i]];
        ib = i;
      }
    }
    /* now check leftovers */
    if (ib == i - 1) /* last two were unique */
      ranks[indx[ib]] = (double)i;
    else
    { /* ended with ties */
      double rnk = method == 0 ? ((i - 1 + ib + 2) / 2.0) : (i - 1);
      for (int j = ib; j <= i - 1; ++j)
        ranks[indx[j]] = rnk;
    }

    return ranks;
  }

  mat RIN_transform(mat &A, int thread_no)
  {
    int M = A.n_rows;
    int N = A.n_cols;

    mat Zr = zeros(M, N);
    parallelFor(
        0, N,
        [&](size_t i)
        {
          vec v = A.col(i);
          vec p = rank_vec(v) / (v.n_elem + 1);

          /*
          uvec row_perm_forward = stable_sort_index(v);
          uvec row_perm = stable_sort_index(row_perm_forward);

          vec p = (row_perm + ones(size(row_perm))) / (row_perm.n_elem + 1);
          */
          vec v_RINT = zeros(size(p));
          for (int j = 0; j < p.n_elem; j++)
          {
            double norm_inv = r8_normal_01_cdf_inverse(p(j));
            v_RINT(j) = norm_inv;
          }

          Zr.col(i) = v_RINT;
        },
        thread_no);
    Zr.replace(datum::nan, 0); // replace each NaN with 0

    return (Zr);
  }

  void convtests(int Bsz, int n, double tol, double svtol, double Smax,
                 double *svratio, double *residuals, int *k, int *converged,
                 double S)
  {
    int j, Len_res = 0;
    for (j = 0; j < Bsz; j++)
    {
      if ((fabs(residuals[j]) < tol * Smax) && (svratio[j] < svtol))
        Len_res++;
    }

    if (Len_res >= n || S == 0)
    {
      *converged = 1;
      return;
    }
    if (*k < n + Len_res)
      *k = n + Len_res;

    if (*k > Bsz - 3)
      *k = Bsz - 3;

    if (*k < 1)
      *k = 1;

    *converged = 0;

    return;
  }

  void orthog(double *X, double *Y, double *T, int xm, int xn, int yn)
  {
    double a = 1, b = 1;
    int inc = 1;
    memset(T, 0, xn * yn * sizeof(double));
    // T = t(X) * Y
    cblas_dgemv(CblasColMajor, CblasTrans, xm, xn, a, X, xm, Y, inc, b, T, inc);
    // Y = Y - X * T
    a = -1.0;
    b = 1.0;
    cblas_dgemv(CblasColMajor, CblasNoTrans, xm, xn, a, X, xm, T, inc, b, Y, inc);
  }

} // namespace ACTIONet
