#include <RcppArmadillo.h>
#include <ACTIONet.h>

#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;

#define ARMA_USE_CXX11_RNG

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List FengSVD(SEXP A, int dim, int iters = 5, int seed = 0) {	
	field<mat> SVD_out;

    if (Rf_isS4(A)) {
		sp_mat tmp = as<arma::sp_mat>(A);
		SVD_out = ACTIONet::FengSVD(tmp, dim, iters, seed);            
    } else {
		mat tmp = as<arma::mat>(A);
		SVD_out = ACTIONet::FengSVD(tmp, dim, iters, seed);            
    } 
        
	List res;
	res["U"] = SVD_out(0);	
	res["sigma"] = SVD_out(1);	
	res["V"] = SVD_out(2);	
		
	return res;
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List HalkoSVD(SEXP A, int dim, int iters = 5, int seed = 0) {	
	field<mat> SVD_out;
	
    if (Rf_isS4(A)) {
		sp_mat tmp = as<arma::sp_mat>(A);
		SVD_out = ACTIONet::HalkoSVD(tmp, dim, iters, seed);            
    } else {
		mat tmp = as<arma::mat>(A);
		SVD_out = ACTIONet::HalkoSVD(tmp, dim, iters, seed);            
    } 
	
	List res;
	res["U"] = SVD_out(0);	
	res["sigma"] = SVD_out(1);	
	res["V"] = SVD_out(2);	
		
	return res;
}

