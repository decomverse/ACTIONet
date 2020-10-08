#ifndef RCPP_UTIL_H
#define RCPP_UTIL_H

namespace ACTIONet {
  Rcpp::NumericVector arma2vec(const T& x);
  arma::vec roll_var(arma::vec &X);
  Rcpp::NumericVector fast_row_sums(SEXP &A);
  Rcpp::NumericVector fast_column_sums(SEXP &A);
  Rcpp::NumericVector computeSparseRowVariances(IntegerVector j, NumericVector val, NumericVector rm, int n);
}

#endif
