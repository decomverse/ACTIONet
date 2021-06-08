#ifndef ARMA_BASE_H
#define ARMA_BASE_H

#include <my_cblas.h>

//#undef ARMA_USE_MKL_TYPES
//#define ARMA_64BIT_WORD
//#define ARMA_BLAS_LONG_LONG
//#define ARMA_USE_FORTRAN_HIDDEN_ARGS

#define ARMA_DONT_USE_WRAPPER
#undef ARMA_BLAS_CAPITALS
#define ARMA_BLAS_UNDERSCORE

#include <arma/armadillo>

using namespace arma;
using namespace std;

#endif
