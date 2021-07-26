#ifndef RCPP_UTIL_H
#define RCPP_UTIL_H

#include <RcppArmadillo.h>

#define stdout_printf Rprintf
#define stderr_printf REprintf
#define stderr_stop stop
#define FLUSH R_FlushConsole()

namespace ACTIONet
{
    template <typename T>
    Rcpp::NumericVector arma2vec(const T &x);
    arma::vec roll_var(arma::vec &X);
    Rcpp::NumericVector fast_row_sums(SEXP &A);
    Rcpp::NumericVector fast_column_sums(SEXP &A);
    Rcpp::NumericVector computeSparseRowVariances(Rcpp::IntegerVector j,
                                                  Rcpp::NumericVector val,
                                                  Rcpp::NumericVector rm, int n);
    arma::sp_mat bind_sparse_mats(sp_mat &A, sp_mat &B, int dim);
} // namespace ACTIONet

#endif
