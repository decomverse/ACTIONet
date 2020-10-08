#include <RcppArmadillo.h>
#include <RCpp_util.h>
#include <ACTIONet.h>


#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

#define ARMA_USE_CXX11_RNG
#define DYNSCHED


template <typename T>
Rcpp::NumericVector arma2vec(const T& x) {
    return Rcpp::NumericVector(x.begin(), x.end());
}

// [[Rcpp::export]]
vec roll_var(vec &X){
    const uword n_max = X.n_elem;
    double xbar = 0, M = 0;
    vec out(n_max);
    double *x = X.begin(), *o = out.begin();

    for(uword n = 1; n <= n_max; ++n, ++x, ++o){
      double tmp = (*x - xbar);
      xbar += (*x - xbar) / n;
      M += tmp * (*x - xbar);
      if(n > 1L)
        *o = M / (n - 1.);
    }

    if(n_max > 0)
      out[0] = std::numeric_limits<double>::quiet_NaN();

    return out;
  }


// [[Rcpp::export]]
NumericVector fast_row_sums(SEXP &A) {	
	vec sum_vec;
    if (Rf_isS4(A)) {
		sp_mat X = as<arma::sp_mat>(A);		
		sum_vec = zeros(X.n_rows);
		
		sp_mat::const_iterator it     = X.begin();
		sp_mat::const_iterator it_end = X.end();
		for(; it != it_end; ++it) {
			sum_vec[it.row()] += (*it);
		}
    } else {
		mat X = as<arma::mat>(A);				
		sum_vec = sum(X, 1);		
    }     

	return(arma2vec(sum_vec));
}

// [[Rcpp::export]]
NumericVector fast_column_sums(SEXP &A) {	
	vec sum_vec;
    if (Rf_isS4(A)) {
		sp_mat X = as<arma::sp_mat>(A);		
		sum_vec = zeros(X.n_rows);
		
		sp_mat::const_iterator it     = X.begin();
		sp_mat::const_iterator it_end = X.end();
		for(; it != it_end; ++it) {
			sum_vec[it.col()] += (*it);
		}
    } else {
		mat X = as<arma::mat>(A);				
		sum_vec = trans(sum(X, 0));
    }     

	return(arma2vec(sum_vec));
}
