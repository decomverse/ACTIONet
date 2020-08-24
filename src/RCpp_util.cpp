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


// [[Rcpp::export]]
arma::vec roll_var(arma::vec &X){
    const arma::uword n_max = X.n_elem;
    double xbar = 0, M = 0;
    arma::vec out(n_max);
    double *x = X.begin(), *o = out.begin();

    for(arma::uword n = 1; n <= n_max; ++n, ++x, ++o){
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
