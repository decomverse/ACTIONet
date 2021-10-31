#ifndef RCPP_UTIL_H
#define RCPP_UTIL_H

#include <RcppArmadillo.h>

#define stdout_printf Rprintf
#define stderr_printf REprintf
#define stderr_stop stop
#define FLUSH R_FlushConsole()

namespace ACTIONet
{
    arma::vec roll_var(arma::vec &X);
    Rcpp::NumericVector computeSparseRowVariances(Rcpp::IntegerVector j, Rcpp::NumericVector val, Rcpp::NumericVector rm, int n);
} // namespace ACTIONet

#endif
