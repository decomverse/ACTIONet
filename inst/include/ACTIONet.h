#ifndef ACTIONet_H
#define ACTIONet_H

#define stdout_printf Rprintf
#define stderr_printf REprintf
#define FLUSH R_FlushConsole()

#define ARMA_DONT_USE_WRAPPER
#undef ARMA_BLAS_CAPITALS
#define ARMA_BLAS_UNDERSCORE

#include <Rinterface.h>
#include <RcppArmadillo.h>
#include <my_cblas.h>

using namespace arma;
using namespace std;

#include <base.h>

#endif
