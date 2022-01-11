#ifndef ACTIONet_H
#define ACTIONet_H

#define stdout_printf printf
#define stderr_printf printf
#define FLUSH fflush(stdout)

#define ARMA_DONT_USE_WRAPPER
#undef ARMA_BLAS_CAPITALS
#define ARMA_BLAS_UNDERSCORE

#include <arma/armadillo>
#include <my_cblas.h>

using namespace arma;
using namespace std;

#include <base.h>

#endif
