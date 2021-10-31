#include <ACTIONet.h>
#include <RCpp_util.h>
#include <RcppArmadillo.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

#define ARMA_USE_CXX11_RNG
#define DYNSCHED

// [[Rcpp::export]]
vec roll_var(vec &X) {
  const uword n_max = X.n_elem;
  double xbar = 0, M = 0;
  vec out(n_max);
  double *x = X.begin(), *o = out.begin();

  for (uword n = 1; n <= n_max; ++n, ++x, ++o) {
    double tmp = (*x - xbar);
    xbar += (*x - xbar) / n;
    M += tmp * (*x - xbar);
    if (n > 1L) *o = M / (n - 1.);
  }

  if (n_max > 0) out[0] = std::numeric_limits<double>::quiet_NaN();

  return out;
}

// Adapted from
// https://github.com/GreenleafLab/MPAL-Single-Cell-2019/blob/master/scRNA_02_Cluster_Disease_w_Reference_v1.R
// [[Rcpp::export]]
Rcpp::NumericVector computeSparseRowVariances(IntegerVector j,
                                              NumericVector val,
                                              NumericVector rm, int n) {
  const int nv = j.size();
  const int nm = rm.size();
  Rcpp::NumericVector rv(nm);
  Rcpp::NumericVector rit(nm);
  int current;
  // Calculate RowVars Initial
  for (int i = 0; i < nv; ++i) {
    current = j(i) - 1;
    rv(current) = rv(current) + (val(i) - rm(current)) * (val(i) - rm(current));
    rit(current) = rit(current) + 1;
  }
  // Calculate Remainder Variance
  for (int i = 0; i < nm; ++i) {
    rv(i) = rv(i) + (n - rit(i)) * rm(i) * rm(i);
  }
  rv = rv / (n - 1);
  return (rv);
}
