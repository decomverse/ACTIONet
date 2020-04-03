#include <RcppArmadillo.h>
#include <ACTIONet.h>

#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;

#define ARMA_USE_CXX11_RNG

//' Computes SVD decomposition
//'
//' This is direct implementation of the randomized SVD algorithm for sparse matrices:
//' Xu Feng, Yuyang Xie, and Yaohang Li, "Fast Randomzied SVD for Sparse Data," in Proc. the 10th Asian Conference on Machine Learning (ACML), Beijing, China, Nov. 2018.
//' 
//' @param A Input matrix (either a "matrix" or "sparseMatrix")
//' @param dim Dimension of SVD decomposition
//' @param iters Number of 
//' @return A named list with U, sigma, and V components
//' @export
//' @examples
//' A = randn(100, 20)
//' SVD.out = FengSVD(A, dim = 2)
//' U = SVD.out$U
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



//' Computes SVD decomposition
//'
//' This is direct implementation of the randomized SVD algorithm:
//' XFrom: N Halko, P. G Martinsson, and J. A Tropp. Finding structure with randomness: Probabilistic algorithms for constructing approximate matrix decompositions. Siam Review, 53(2):217-288, 2011.
//' 
//' @param A Input matrix (either a "matrix" or "sparseMatrix")
//' @param dim Dimension of SVD decomposition
//' @param iters Number of 
//' @return A named list with U, sigma, and V components
//' @export
//' @examples
//' A = randn(100, 20)
//' SVD.out = HalkoSVD(A, dim = 2)
//' U = SVD.out$U
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

