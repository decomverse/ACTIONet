#ifndef ARMA_BASE_H
#define ARMA_BASE_H


//#undef FC_LEN_T


#include <cblas.h>

//#undef ARMA_USE_MKL_TYPES
//#define ARMA_64BIT_WORD
//#define ARMA_BLAS_LONG_LONG
//#define ARMA_USE_FORTRAN_HIDDEN_ARGS


#define ARMA_DONT_USE_WRAPPER
#undef ARMA_BLAS_CAPITALS
#define ARMA_BLAS_UNDERSCORE


#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace std;


#endif
