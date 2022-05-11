#include <ACTIONet.h>
#include <RcppArmadillo.h>


/*
#include "Rforceatlas_types.h"
#include "params.hpp"
#include "graph.hpp"
#include "work.hpp"
*/

#include "uwot/coords.h"
#include "uwot/epoch.h"
#include "uwot/gradient.h"
#include "uwot/optimize.h"
#include "uwot/sampler.h"
#include "uwot/rng.h"
#include "uwot/rparallel.h"

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

#define ARMA_USE_CXX11_RNG
#define DYNSCHED

template <typename T>
Rcpp::NumericVector arma2vec(const T &x)
{
  return Rcpp::NumericVector(x.begin(), x.end());
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List run_ACTION_muV(const List &S, int k_min, int k_max, vec alpha, double lambda = 1, int AA_iters = 50, int Opt_iters = 0, int thread_no = 0)
{

  int n_list = S.size();
  vector<mat> cell_signatures(n_list);
  for (int i = 0; i < n_list; i++)
  {
    cell_signatures[i] = (as<mat>(S[i]));
  }

  full_trace run_trace = ACTIONet::runACTION_muV(cell_signatures, k_min, k_max, alpha, lambda, AA_iters, Opt_iters, thread_no);

  List res;

  List H_consensus(k_max);
  for (int kk = k_min; kk <= k_max; kk++)
  {
    H_consensus[kk - 1] = run_trace.H_consensus[kk];
  }
  res["H_consensus"] = H_consensus;

  char ds_name[128];
  for (int i = 0; i < n_list; i++)
  {
    List individual_trace;

    List H_primary(k_max);
    for (int kk = k_min; kk <= k_max; kk++)
    {
      H_primary[kk - 1] = run_trace.indiv_trace[kk].H_primary[i];
    }
    individual_trace["H_primary"] = H_primary;

    List H_secondary(k_max);
    for (int kk = k_min; kk <= k_max; kk++)
    {
      H_secondary[kk - 1] = run_trace.indiv_trace[kk].H_secondary[i];
    }
    individual_trace["H_secondary"] = H_secondary;

    List C_primary(k_max);
    for (int kk = k_min; kk <= k_max; kk++)
    {
      C_primary[kk - 1] = run_trace.indiv_trace[kk].C_primary[i];
    }
    individual_trace["C_primary"] = C_primary;

    List C_consensus(k_max);
    for (int kk = k_min; kk <= k_max; kk++)
    {
      C_consensus[kk - 1] = run_trace.indiv_trace[kk].C_consensus[i];
    }
    individual_trace["C_consensus"] = C_consensus;

    sprintf(ds_name, "View%d_trace", i + 1);
    res[ds_name] = individual_trace;
  }

  return res;
}

// set seed
// [[Rcpp::export]]
void set_seed(double seed)
{
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(std::floor(std::fabs(seed)));
}

//' Computes SVD decomposition
//'
//' This is direct implementation of the randomized SVD algorithm:
//' From: IRLBA R Package
//'
//' @param A Input matrix ("sparseMatrix")
//' @param dim Dimension of SVD decomposition
//' @param iters Number of iterations (default=5)
//' @param seed Random seed (default=0)
//'
//' @return A named list with U, sigma, and V components
//'
//' @examples
//' A = randn(100, 20)
//' SVD.out = IRLBA_SVD(A, dim = 2)
//' U = SVD.out$U
// [[Rcpp::export]]
List IRLB_SVD(sp_mat &A, int dim, int iters = 1000, int seed = 0, int verbose = 1)
{
  field<mat> SVD_out = ACTIONet::IRLB_SVD(A, dim, iters, seed, verbose);

  List res;

  res["u"] = SVD_out(0);
  res["d"] = SVD_out(1);
  res["v"] = SVD_out(2);

  return res;
}

//' Computes SVD decomposition
//'
//' This is direct implementation of the randomized SVD algorithm:
//' From: IRLBA R Package
//'
//' @param A Input matrix ("sparseMatrix")
//' @param dim Dimension of SVD decomposition
//' @param iters Number of iterations (default=5)
//' @param seed Random seed (default=0)
//'
//' @return A named list with U, sigma, and V components
//'
//' @examples
//' A = randn(100, 20)
//' SVD.out = IRLBA_SVD_full(A, dim = 2)
//' U = SVD.out$U
// [[Rcpp::export]]
List IRLB_SVD_full(mat &A, int dim, int iters = 1000, int seed = 0, int verbose = 1)
{
  field<mat> SVD_out = ACTIONet::IRLB_SVD(A, dim, iters, seed, verbose);

  List res;

  res["u"] = SVD_out(0);
  res["d"] = SVD_out(1);
  res["v"] = SVD_out(2);

  return res;
}

//' Computes SVD decomposition
//'
//' This is direct implementation of the randomized SVD algorithm for sparse
// matrices: ' Xu Feng, Yuyang Xie, and Yaohang Li, "Fast Randomzied SVD for
// Sparse Data," in Proc. the 10th Asian Conference on Machine Learning (ACML),
// Beijing, China, Nov. 2018.
//'
//' @param A Input matrix ("sparseMatrix")
//' @param dim Dimension of SVD decomposition
//' @param iters Number of iterations (default=5)
//' @param seed Random seed (default=0)
//'
//' @return A named list with U, sigma, and V components
//'
//' @examples
//' A = randn(100, 20)
//' SVD.out = FengSVD(A, dim = 2)
//' U = SVD.out$U
// [[Rcpp::export]]
List FengSVD(sp_mat &A, int dim, int iters = 5, int seed = 0, int verbose = 1)
{
  field<mat> SVD_out = ACTIONet::FengSVD(A, dim, iters, seed, verbose);

  List res;

  res["u"] = SVD_out(0);
  res["d"] = SVD_out(1);
  res["v"] = SVD_out(2);

  return res;
}
//' Computes SVD decomposition
//'
//' This is direct implementation of the randomized SVD algorithm for sparse
// matrices: ' Xu Feng, Yuyang Xie, and Yaohang Li, "Fast Randomzied SVD for
// Sparse Data," in Proc. the 10th Asian Conference on Machine Learning (ACML),
// Beijing, China, Nov. 2018.
//'
//' @param A Input matrix ("matrix")
//' @param dim Dimension of SVD decomposition
//' @param iters Number of iterations (default=5)
//' @param seed Random seed (default=0)
//'
//' @return A named list with U, sigma, and V components
//'
//' @examples
//' A = randn(100, 20)
//' SVD.out = FengSVD(A, dim = 2)
//' U = SVD.out$U
// [[Rcpp::export]]
List FengSVD_full(mat &A, int dim, int iters = 5, int seed = 0, int verbose = 1)
{
  field<mat> SVD_out = ACTIONet::FengSVD(A, dim, iters, seed, verbose);

  List res;

  res["u"] = SVD_out(0);
  res["d"] = SVD_out(1);
  res["v"] = SVD_out(2);

  return res;
}

//' Computes SVD decomposition
//'
//' This is direct implementation of the randomized SVD algorithm:
//' From: N Halko, P. G Martinsson, and J. A Tropp. Finding structure with
// randomness: Probabilistic algorithms for constructing approximate matrix
// decompositions. Siam Review, 53(2):217-288, 2011.
//'
//' @param A Input matrix ("sparseMatrix")
//' @param dim Dimension of SVD decomposition
//' @param iters Number of iterations (default=5)
//' @param seed Random seed (default=0)
//'
//' @return A named list with U, sigma, and V components
//'
//' @examples
//' A = randn(100, 20)
//' SVD.out = HalkoSVD(A, dim = 2)
//' U = SVD.out$U
// [[Rcpp::export]]
List HalkoSVD(sp_mat &A, int dim, int iters = 5, int seed = 0, int verbose = 1)
{
  field<mat> SVD_out = ACTIONet::HalkoSVD(A, dim, iters, seed, verbose);

  List res;

  res["u"] = SVD_out(0);
  res["d"] = SVD_out(1);
  res["v"] = SVD_out(2);

  return res;
}

//' Computes SVD decomposition
//'
//' This is direct implementation of the randomized SVD algorithm:
//' From: N Halko, P. G Martinsson, and J. A Tropp. Finding structure with
// randomness: Probabilistic algorithms for constructing approximate matrix
// decompositions. Siam Review, 53(2):217-288, 2011.
//'
//' @param A Input matrix ("matrix")
//' @param dim Dimension of SVD decomposition
//' @param iters Number of iterations (default=5)
//' @param seed Random seed (default=0)
//'
//' @return A named list with U, sigma, and V components
//'
//' @examples
//' A = randn(100, 20)
//' SVD.out = HalkoSVD(A, dim = 2)
//' U = SVD.out$U
// [[Rcpp::export]]
List HalkoSVD_full(mat &A, int dim, int iters = 5, int seed = 0, int verbose = 1)
{
  field<mat> SVD_out = ACTIONet::HalkoSVD(A, dim, iters, seed, verbose);

  List res;

  res["u"] = SVD_out(0);
  res["d"] = SVD_out(1);
  res["v"] = SVD_out(2);

  return res;
}

//' Computes reduced kernel matrix for a given (single-cell) profile
//'
//' @param S Input matrix ("sparseMatrix")
//' @param reduced_dim Dimension of the reduced kernel matrix (default=50)
//' @param iters Number of SVD iterations (default=5)
//' @param seed Random seed (default=0)
//' @param reduction_algorithm Kernel reduction algorithm. Currently only ACTION
// method (1) is implemented (default=1) ' @param SVD_algorithm SVD algorithm to
// use. Currently supported methods are Halko (1) and Feng (2) (default=1)
//'
//' @return A named list with S_r, V, lambda, and exp_var. \itemize{
//' \item S_r: reduced kernel matrix of size reduced_dim x #samples.
//' \item V: Associated left singular-vectors (useful for reconstructing
// discriminative scores for features, such as genes). ' \item lambda, exp_var:
// Summary statistics of the sigular-values. ' }
//'
//' @examples
//' S = logcounts(sce)
//' reduction.out = reduce(S, reduced_dim = 50)
//' S_r = reduction.out$S_r
// [[Rcpp::export]]
List reduce_kernel(sp_mat &S, int reduced_dim = 50, int iter = 5, int seed = 0,
                   int SVD_algorithm = 0, bool prenormalize = false, int verbose = 1)
{
  field<mat> reduction = ACTIONet::reduce_kernel(S, reduced_dim, iter, seed,
                                                 SVD_algorithm, prenormalize, verbose);

  List res;
  res["V"] = reduction(0);

  vec sigma = reduction(1).col(0);
  res["sigma"] = sigma;

  mat V = reduction(2);
  // printf("%d x %d\n", V.n_rows, V.n_cols);
  for (int i = 0; i < V.n_cols; i++)
  {
    vec v = V.col(i) * sigma(i);
    V.col(i) = v;
  }
  V = trans(V);
  res["S_r"] = V.eval();

  res["A"] = reduction(3);
  res["B"] = reduction(4);

  return res;
}

//' Computes reduced kernel matrix for a given (single-cell) profile
//'
//' @param S Input matrix ("matrix")
//' @param reduced_dim Dimension of the reduced kernel matrix (default=50)
//' @param iters Number of SVD iterations (default=5)
//' @param seed Random seed (default=0)
//' @param reduction_algorithm Kernel reduction algorithm. Currently only ACTION
// method (1) is implemented (default=1) ' @param SVD_algorithm SVD algorithm to
// use. Currently supported methods are Halko (1) and Feng (2) (default=1)
//'
//' @return A named list with S_r, V, lambda, and exp_var. \itemize{
//' \item S_r: reduced kernel matrix of size reduced_dim x #samples.
//' \item V: Associated left singular-vectors (useful for reconstructing
// discriminative scores for features, such as genes). ' \item lambda, exp_var:
// Summary statistics of the sigular-values. ' }
//'
//' @examples
//' S = logcounts(sce)
//' reduction.out = reduce(S, reduced_dim = 50)
//' S_r = reduction.out$S_r
// [[Rcpp::export]]
List reduce_kernel_full(mat &S, int reduced_dim = 50, int iter = 5,
                        int seed = 0, int SVD_algorithm = 0,
                        bool prenormalize = false, int verbose = 1)
{
  field<mat> reduction = ACTIONet::reduce_kernel(S, reduced_dim, iter, seed,
                                                 SVD_algorithm, prenormalize, verbose);

  List res;
  res["V"] = reduction(0);

  vec sigma = reduction(1).col(0);
  res["sigma"] = sigma;

  mat V = reduction(2);
  // printf("%d x %d\n", V.n_rows, V.n_cols);
  for (int i = 0; i < V.n_cols; i++)
  {
    vec v = V.col(i) * sigma(i);
    V.col(i) = v;
  }
  V = trans(V);
  res["S_r"] = V.eval();

  res["A"] = reduction(3);
  res["B"] = reduction(4);

  return res;
}

//' Solves min_{X} (|| AX - B ||) s.t. simplex constraint
//'
//' @param A Input matrix
//' @param B Input matrix
//'
//' @return X Solution
//'
//' @examples
//' C = ACTION.out$C[[10]]
//' A = S_r %*% C
//' B = S_r
//' H = run_simplex_regression(A, B)
// [[Rcpp::export]]
mat run_simplex_regression(mat &A, mat &B, bool computeXtX = false)
{
  mat X = ACTIONet::run_simplex_regression(A, B, computeXtX);

  return X;
}

// [[Rcpp::export]]
mat run_simplex_regression_FW(mat &A, mat &B, int max_iter = -1, double min_diff = 0.01)
{
  mat X = ACTIONet::run_simplex_regression_FW(A, B, max_iter, min_diff);

  return X;
}


//' Runs Successive Projection Algorithm (SPA) to solve separable NMF
//'
//' @param A Input matrix
//' @param k Number of columns to select
//'
//' @return A named list with entries 'selected_columns' and 'norms'
//' @examples
//' H = run_SPA(S_r, 10)
// [[Rcpp::export]]
List run_SPA(mat &A, int k)
{
  ACTIONet::SPA_results res = ACTIONet::run_SPA(A, k);
  uvec selected_columns = res.selected_columns;

  vec cols(k);
  for (int i = 0; i < k; i++)
  {
    cols[i] = selected_columns[i] + 1;
  }

  List out;
  out["selected_columns"] = cols;
  out["norms"] = res.column_norms;

  return out;
}

//' Runs Successive Projection Algorithm (SPA) to solve separable NMF
//'
//' @param A Input matrix
//' @param k Number of columns to select
//'
//' @return A named list with entries 'selected_columns' and 'norms'
//' @examples
//' H = run_SPA(S_r, 10)
// [[Rcpp::export]]
List run_SPA_rows_sparse(sp_mat &A, int k)
{
  ACTIONet::SPA_results res = ACTIONet::run_SPA_rows_sparse(A, k);
  uvec selected_columns = res.selected_columns;

  vec cols(k);
  for (int i = 0; i < k; i++)
  {
    cols[i] = selected_columns[i] + 1;
  }

  List out;
  out["selected_rows"] = cols;
  out["norms"] = res.column_norms;

  return out;
}

//' Runs multi-level ACTION decomposition method
//'
//' @param S_r Reduced kernel matrix
//' @param k_min Minimum number of archetypes to consider (default=2)
//' @param k_max Maximum number of archetypes to consider, or "depth" of
// decomposition (default=30) ' @param thread_no Number of parallel threads
//(default = 0) ' @param max_it,min_delta Convergence parameters for archetypal
// analysis
//'
//' @return A named list with entries 'C' and 'H', each a list for different
// values of k ' @examples ' ACTION.out = run_ACTION(S_r, k_max = 10) ' H8 =
// ACTION.out$H[[8]] ' cell.assignments = apply(H8, 2, which.max)
// [[Rcpp::export]]
List run_ACTION(mat &S_r, int k_min = 2, int k_max = 30, int thread_no = 0,
                int max_it = 100, double min_delta = 1e-6)
{
  ACTIONet::ACTION_results trace =
      ACTIONet::run_ACTION(S_r, k_min, k_max, thread_no, max_it, min_delta);

  List res;

  List C(k_max);
  for (int i = k_min; i <= k_max; i++)
  {
    mat cur_C = trace.C[i];
    //cur_C.transform( [](double val) { return (val < 1e-5? 0:val); } );
    //cur_C = round(cur_C*1e5)/1e-5;
    //cur_C = normalise(cur_C, 1);
    C[i - 1] = cur_C;
  }
  res["C"] = C;

  List H(k_max);
  for (int i = k_min; i <= k_max; i++)
  {
    mat cur_H = trace.H[i];
    //cur_H.transform( [](double val) { return (val < 1e-5? 0:val); } );
    //cur_H = normalise(cur_H, 1);
    H[i - 1] = cur_H;
  }
  res["H"] = H;

  return res;
}

//' Runs multi-level ACTION decomposition method
//'
//' @param S_r Reduced kernel matrix
//' @param k_min Minimum number of archetypes to consider (default=2)
//' @param k_max Maximum number of archetypes to consider, or "depth" of
// decomposition (default=30) ' @param max_it,min_delta Convergence parameters
// for archetypal analysis ' @param max_trial Maximum number of trials before
// termination
//'
//' @return A named list with entries 'C' and 'H', each a list for different
// values of k ' @examples ' ACTION.out = run_ACTION_plus(S_r, k_max = 10) ' H8
// = ACTION.out$H[[8]] ' cell.assignments = apply(H8, 2, which.max)
// [[Rcpp::export]]
List run_ACTION_plus(mat &S_r, int k_min = 2, int k_max = 30, int max_it = 100,
                     double min_delta = 1e-6, int max_trial = 3)
{
  ACTIONet::ACTION_results trace = ACTIONet::run_ACTION_plus(
      S_r, k_min, k_max, max_it, min_delta, max_trial);

  List res;

  List C(trace.H.n_elem - 1);
  for (int i = k_min; i < trace.H.n_elem; i++)
  {
    C[i - 1] = trace.C[i];
  }
  res["C"] = C;

  List H(trace.H.n_elem - 1);
  for (int i = k_min; i < trace.H.n_elem; i++)
  {
    H[i - 1] = trace.H[i];
  }
  res["H"] = H;

  return res;
}

//' Runs basic archetypal analysis
//'
//' @param A Inpu matrix
//' @param W0 Starting archetypes
//' @param max_it,min_delta Convergence parameters for archetypal analysis
//'
//' @return A named list with entries 'C' and 'H', each a list for different
// values of k ' @examples ' S_r = t(reducedDims(ace)$ACTION) ' SPA.out =
// run_SPA(S_r, 10) ' W0 = S_r[, SPA.out$selected_columns] ' AA.out =
// run_AA(S_r, W0) ' H = AA.out$H ' cell.assignments = apply(H, 2, which.max)
// [[Rcpp::export]]
List run_AA(mat &A, mat &W0, int max_it = 100, double min_delta = 1e-6)
{
  field<mat> res = ACTIONet::run_AA(A, W0, max_it, min_delta);

  List out;
  out["C"] = res(0);
  out["H"] = res(1);
  out["W"] = A * res(0);

  return out;
}

//' Runs multi-level Online ACTION decomposition method (under development)
//'
//' @param S_r Reduced kernel matrix
//' @param k_min Minimum number of archetypes to consider (default=2)
//' @param k_max Maximum number of archetypes to consider, or "depth" of
// decomposition (default=30) ' @param samples List of sampled cells to use for
// updating archetype decomposition ' @param thread_no Number of parallel
// threads (default = 0)
//'
//' @return A named list with entries 'C' and 'H', each a list for different
// values of k ' @examples ' ACTION.out = run_online_ACTION(S_r, k_max = 10)
// [[Rcpp::export]]
List run_online_ACTION(mat &S_r, field<uvec> samples, int k_min = 2,
                       int k_max = 30, int thread_no = 0)
{
  ACTIONet::Online_ACTION_results trace =
      ACTIONet::run_online_ACTION(S_r, samples, k_min, k_max, thread_no);

  List res;

  List A(k_max);
  for (int i = k_min; i <= k_max; i++)
  {
    A[i - 1] = trace.A[i];
  }
  res["A"] = A;

  List B(k_max);
  for (int i = k_min; i <= k_max; i++)
  {
    B[i - 1] = trace.B[i];
  }
  res["B"] = B;

  List C(k_max);
  for (int i = k_min; i <= k_max; i++)
  {
    C[i - 1] = trace.C[i];
  }
  res["C"] = C;

  List D(k_max);
  for (int i = k_min; i <= k_max; i++)
  {
    D[i - 1] = trace.D[i];
  }
  res["D"] = D;

  return res;
}

//' Runs multi-level weighted ACTION decomposition method (under development)
//'
//' @param S_r Reduced kernel matrix
//' @param w Weight vector for each observation
//' @param k_min Minimum number of archetypes to consider (default=2)
//' @param k_max Maximum number of archetypes to consider, or "depth" of
// decomposition (default=30) ' @param thread_no Number of parallel threads
//(default=0)
//'
//' @return A named list with entries 'C' and 'H', each a list for different
// values of k ' @examples ' ACTION.out = run_weighted_ACTION(S_r, w, k_max =
// 20)
// [[Rcpp::export]]
List run_weighted_ACTION(mat &S_r, vec w, int k_min = 2, int k_max = 30,
                         int thread_no = 0, int max_it = 50,
                         double min_delta = 1e-16)
{
  ACTIONet::ACTION_results trace = ACTIONet::run_weighted_ACTION(
      S_r, w, k_min, k_max, thread_no, max_it, min_delta);

  List res;

  List C(k_max);
  for (int i = k_min; i <= k_max; i++)
  {
    C[i - 1] = trace.C[i];
  }
  res["C"] = C;

  List H(k_max);
  for (int i = k_min; i <= k_max; i++)
  {
    H[i - 1] = trace.H[i];
  }
  res["H"] = H;

  return res;
}

//' Filters multi-level archetypes and concatenate filtered archetypes.
//' (Pre-ACTIONet archetype processing)
//'
//' @param C_trace,H_trace Output of ACTION
//' @param min_specificity_z_threshold Defines the stringency of pruning
// nonspecific archetypes. ' The larger the value, the more archetypes will be
// filtered out (default=-1)
//'
//' @return A named list: \itemize{
//' \item selected_archs: List of final archetypes that passed the
// filtering/pruning step. ' \item C_stacked,H_stacked: Horizontal/Vertical
// concatenation of filtered C and H matrices, respectively. ' } ' @examples ' S
//= logcounts(sce) ' reduction.out = reduce(S, reduced_dim = 50) ' S_r =
// reduction.out$S_r ' ACTION.out = run_ACTION(S_r, k_max = 10) '
// reconstruction.out = reconstruct_archetypes(S, ACTION.out$C, ACTION.out$H)
// [[Rcpp::export]]
List prune_archetypes(const List &C_trace, const List &H_trace,
                      double min_specificity_z_threshold = -3,
                      int min_cells = 3)
{
  int n_list = H_trace.size();
  field<mat> C_trace_vec(n_list + 1);
  field<mat> H_trace_vec(n_list + 1);
  for (int i = 0; i < n_list; i++)
  {
    if (Rf_isNull(H_trace[i]))
    {
      continue;
    }
    C_trace_vec[i + 1] = (as<mat>(C_trace[i]));
    H_trace_vec[i + 1] = (as<mat>(H_trace[i]));
  }

  ACTIONet::multilevel_archetypal_decomposition results =
      ACTIONet::prune_archetypes(C_trace_vec, H_trace_vec,
                                 min_specificity_z_threshold, min_cells);

  List out_list;

  for (int i = 0; i < results.selected_archs.n_elem; i++)
    results.selected_archs[i]++;
  out_list["selected_archs"] = results.selected_archs;

  out_list["C_stacked"] = results.C_stacked;
  out_list["H_stacked"] = results.H_stacked;

  return out_list;
}

//' Identifies and aggregates redundant archetypes into equivalent classes
//' (Post-ACTIONet archetype processing)
//'
//' @param G Adjacency matrix of the ACTIONet graph
//' @param S_r Reduced kernel profile
//' @param archetypes Archetype profile (S*C)
//' @param C_stacked,H_stacked Output of reconstruct_archetypes()
//' @param minPoints, minClusterSize, outlier_threshold HDBSCAN parameters
//' @param reduced_dim Kernel reduction
//'
//' @return A named list: \itemize{
//' \item archetype_groups: Equivalent classes of archetypes (non-redundant)
//' \item C_unified,H_unified: C and H matrices of unified archetypes
//' \item sample_assignments: Assignment of samples/cells to unified archetypes
//' }
//' @examples
//' prune.out = prune_archetypes(ACTION.out$C, ACTION.out$H)
//'	G = buildNetwork(prune.out$H_stacked)
//' unification.out = unify_archetypes(G, S_r, prune.out$C_stacked,
// prune.out$H_stacked) ' cell.clusters = unification.out$sample_assignments
// [[Rcpp::export]]
List unify_archetypes(mat &S_r, mat &C_stacked, mat &H_stacked,
                      double violation_threshold = 0.0,
                      double backbone_density = 0.5, int outlier_threshold = 1,
                      int thread_no = 0)
{
  ACTIONet::unification_results results = ACTIONet::unify_archetypes(S_r, C_stacked, H_stacked, backbone_density, outlier_threshold, thread_no);

  List out_list;

  for (int i = 0; i < results.selected_archetypes.n_elem; i++)
    results.selected_archetypes[i]++;
  out_list["selected_archetypes"] = results.selected_archetypes;

  out_list["C_unified"] = results.C_unified;
  out_list["H_unified"] = results.H_unified;

  for (int i = 0; i < results.assigned_archetypes.n_elem; i++)
    results.assigned_archetypes[i]++;

  out_list["assigned_archetypes"] = results.assigned_archetypes;
  out_list["arch_membership_weights"] = results.arch_membership_weights;

  out_list["ontology"] = results.dag_adj;
  out_list["ontology_node_attributes"] = results.dag_node_annotations;

  return out_list;
}

//' Builds an interaction network from the multi-level archetypal decompositions
//'
//' @param H_stacked Output of the prune_archetypes() function.
//' @param density Overall density of constructed graph. The higher the density,
// the more edges are retained (default = 1.0). ' @param thread_no Number of
// parallel threads (default = 0). ' @param mutual_edges_only Symmetrization
// strategy for nearest-neighbor edges. ' If it is true, only mutual
// nearest-neighbors are returned (default=TRUE).
//'
//' @return G Adjacency matrix of the ACTIONet graph.
//'
//' @examples
//' prune.out = prune_archetypes(ACTION.out$C, ACTION.out$H)
//'	G = buildNetwork(prune.out$H_stacked)
// [[Rcpp::export]]
sp_mat buildNetwork(mat H, string algorithm = "k*nn", string distance_metric = "jsd", double density = 1.0, int thread_no = 0,
                    bool mutual_edges_only = true, int k = 10, int ef_construction = 200, int ef = 200)
{

  double M = 16;
  sp_mat G = ACTIONet::buildNetwork(H, algorithm, distance_metric, density, thread_no, M,
                                    ef_construction, ef, mutual_edges_only, k);

  return G;
}

// [[Rcpp::export]]
sp_mat build_knn(mat H, string distance_metric = "jsd", double k = 10, int thread_no = 0,
                 bool mutual_edges_only = true)
{

  double M = 16, ef_construction = 200, ef = 200;

  sp_mat G = ACTIONet::buildNetwork(H, "knn", distance_metric, 1, thread_no, M,
                                    ef_construction, ef, mutual_edges_only, k);

  return G;
}

//' Performs stochastic force-directed layout on the input graph (ACTIONet)
//'
//' @param G Adjacency matrix of the ACTIONet graph
//' @param S_r Reduced kernel matrix (is used for reproducible initialization).
//' @param compactness_level A value between 0-100, indicating the compactness
// of ACTIONet layout (default=50) ' @param n_epochs Number of epochs for SGD
// algorithm (default=100). ' @param thread_no Number of threads (default = 0).
//'
//' @return A named list \itemize{
//' \item coordinates 2D coordinates of vertices.
//' \item coordinates_3D 3D coordinates of vertices.
//' \item colors De novo color of nodes inferred from their 3D embedding.
//' }
//'
//' @examples
//'	G = buildNetwork(prune.out$H_stacked)
//'	vis.out = layoutNetwrok(G, S_r)
// [[Rcpp::export]]
List layoutNetwrok(sp_mat &G, mat initial_position, string alg_name = "umap", float spread = 1.0, float min_dist = 0.01, string opt_name = "adam", 
                                unsigned int n_epochs = 1000,
                                int seed = 0, int thread_no = 0) {

  field<mat> res = ACTIONet::layoutNetwork_xmap(G, initial_position, alg_name, spread, min_dist, opt_name, n_epochs, seed, thread_no);

  List out_list;
  out_list["coordinates"] = res(0);
  out_list["coordinates_3D"] = res(1);
  out_list["colors"] = res(2);

  return out_list;
}

//' Encrypts a set of given input ids
//'
//' @param ids List of input string ids
//' @param pass Pass phrase to use for encryption
//'
//' @return A string array of encoded ids
//'
// [[Rcpp::export]]
vector<string> encode_ids(vector<string> ids, string pass)
{
  vector<string> encoded_ids(ids.size());

  cryptor::set_key(pass);
  for (int i = 0; i < ids.size(); i++)
  {
    auto enc = cryptor::encrypt(ids[i]);
    encoded_ids[i] = enc;
  }

  return encoded_ids;
}

//' Decrypts a set of given encrypted ids
//'
//' @param encoded_ids List of encrypted string ids
//' @param pass Pass phrase to use for decryption
//'
//' @return A string array of decrypted ids
//'
// [[Rcpp::export]]
vector<string> decode_ids(vector<string> encoded_ids, string pass)
{
  vector<string> decoded_ids(encoded_ids.size());

  cryptor::set_key(pass);
  for (int i = 0; i < encoded_ids.size(); i++)
  {
    auto dec = cryptor::decrypt(encoded_ids[i]);
    decoded_ids[i] = dec;
  }

  return decoded_ids;
}

//' Aggregate matrix wiithin groups
//'
//' @param S matrix of type "dgCMatrix"
//' @param sample_assignments Vector of column groupings. Group labels must be continuous integers or coercible to such.
//'
//' @return S matrix with columns of values aggregated within each group of sample_assignments
//'
// [[Rcpp::export]]
mat compute_grouped_rowsums(sp_mat &S,
                                    arma::Col<unsigned long long> sample_assignments)
{
  mat pb = ACTIONet::compute_grouped_rowsums(S, sample_assignments);

  return pb;
}

//' Aggregate matrix wiithin groups
//'
//' @param S matrix
//' @param sample_assignments Vector of column groupings. Group labels must be continuous integers or coercible to such.
//'
//' @return S matrix with columns of values aggregated within each group of sample_assignments
//'
// [[Rcpp::export]]
mat compute_grouped_rowsums_full(mat &S,
                                         arma::Col<unsigned long long> sample_assignments)
{
  mat pb = ACTIONet::compute_grouped_rowsums(S, sample_assignments);

  return pb;
}

//' Average matrix wiithin groups
//'
//' @param S matrix of type "dgCMatrix"
//' @param sample_assignments Vector of column groupings. Group labels must be continuous integers or coercible to such.
//'
//' @return S matrix with columns of values average within each group of sample_assignments
//'
// [[Rcpp::export]]
mat compute_grouped_rowmeans(sp_mat &S, arma::Col<unsigned long long> sample_assignments)
{
  mat pb = ACTIONet::compute_grouped_rowmeans(S, sample_assignments);

  return pb;
}

//' Average matrix wiithin groups
//'
//' @param S matrix
//' @param sample_assignments Vector of column groupings. Group labels must be continuous integers or coercible to such.
//'
//' @return S matrix with columns of values average within each group of sample_assignments
//'
// [[Rcpp::export]]
mat compute_grouped_rowmeans_full(mat &S, arma::Col<unsigned long long> sample_assignments)
{
  mat pb = ACTIONet::compute_grouped_rowmeans(S, sample_assignments);

  return pb;
}

// [[Rcpp::export]]
mat compute_grouped_rowvars(sp_mat &S, arma::Col<unsigned long long> sample_assignments)
{
  mat pb = ACTIONet::compute_grouped_rowvars(S, sample_assignments);

  return pb;
}

// [[Rcpp::export]]
mat compute_grouped_rowvars_full(mat &S, arma::Col<unsigned long long> sample_assignments)
{
  mat pb = ACTIONet::compute_grouped_rowvars(S, sample_assignments);

  return pb;
}

// [[Rcpp::export]]
mat compute_pseudo_bulk_per_archetype(sp_mat &S, mat &H)
{
  mat pb = ACTIONet::compute_pseudo_bulk_per_archetype(S, H);

  return pb;
}

// [[Rcpp::export]]
mat compute_pseudo_bulk_per_archetype_full(mat &S,
                                           mat &H)
{
  mat pb = ACTIONet::compute_pseudo_bulk_per_archetype(S, H);

  return pb;
}

// [[Rcpp::export]]
field<mat> compute_pseudo_bulk_per_archetype_and_ind(
    sp_mat &S, mat &H,
    arma::Col<unsigned long long> individuals)
{
  field<mat> pbs_list =
      ACTIONet::compute_pseudo_bulk_per_archetype_and_ind(S, H, individuals);

  return pbs_list;
}

// [[Rcpp::export]]
field<mat> compute_pseudo_bulk_per_archetype_and_ind_full(
    mat &S, mat &H,
    arma::Col<unsigned long long> individuals)
{
  field<mat> pbs_list =
      ACTIONet::compute_pseudo_bulk_per_archetype_and_ind(S, H, individuals);

  return pbs_list;
}

//' Renormalized input matrix to minimize differences in means
//'
//' @param S Input matrix
//' @param sample_assignments Any primary grouping - typically based on cell
// type/state (it has to be in {1, ..., k1})
//'
//' @return A list with the first entry being the renormalized input matrix
//'
//' @examples
//' prune.out = prune_archetypes(ACTION.out$C, ACTION.out$H)
//'	G = buildNetwork(prune.out$H_stacked)
//' unification.out = unify_archetypes(G, S_r, prune.out$C_stacked,
// prune.out$H_stacked) ' cell.clusters = unification.out$sample_assignments '
// S.norm = renormalize_input_matrix(S, cell.clusters)
// [[Rcpp::export]]
sp_mat renormalize_input_matrix(
    sp_mat &S, arma::Col<unsigned long long> sample_assignments)
{
  sp_mat S_norm = ACTIONet::renormalize_input_matrix(S, sample_assignments);

  return (S_norm);
}

//' Renormalized input matrix to minimize differences in means
//'
//' @param S Input matrix ("matrix" type)
//' @param sample_assignments Any primary grouping - typically based on cell
// type/state (it has to be in {1, ..., k1})
//'
//' @return A list with the first entry being the renormalized input matrix
//'
//' @examples
//' prune.out = prune_archetypes(ACTION.out$C, ACTION.out$H)
//'	G = buildNetwork(prune.out$H_stacked)
//' unification.out = unify_archetypes(G, S_r, prune.out$C_stacked,
// prune.out$H_stacked) ' cell.clusters = unification.out$sample_assignments '
// S.norm = renormalize_input_matrix(S, cell.clusters)
// [[Rcpp::export]]
mat renormalize_input_matrix_full(
    mat &S, arma::Col<unsigned long long> sample_assignments)
{
  mat S_norm = ACTIONet::renormalize_input_matrix(S, sample_assignments);

  return (S_norm);
}

//' Compute feature specificity (from archetype footprints and binary input)
//'
//' @param S Input matrix (sparseMatrix - binary)
//' @param H A soft membership matrix - Typically H_unified from the
// unify_archetypes() function.
//'
//' @return A list with the over/under-logPvals
//'
//' @examples
//'	logPvals.list = compute_archetype_feature_specificity_bin(S.bin,
// unification.out$H_unified) ' specificity.scores =
// logPvals.list$upper_significance
// [[Rcpp::export]]
List compute_archetype_feature_specificity_bin(sp_mat &S, mat &H, int thread_no = 0)
{
  field<mat> res = ACTIONet::compute_feature_specificity_bin(S, H, thread_no);

  List out_list;
  out_list["archetypes"] = res(0);
  out_list["upper_significance"] = res(1);
  out_list["lower_significance"] = res(2);

  return (out_list);
}

//' Compute feature specificity (from archetype footprints)
//'
//' @param S Input matrix (sparseMatrix)
//' @param H A soft membership matrix - Typically H_unified from the
// unify_archetypes() function.
//'
//' @return A list with the over/under-logPvals
//'
//' @examples
//' prune.out = prune_archetypes(ACTION.out$C, ACTION.out$H)
//'	G = buildNetwork(prune.out$H_stacked)
//' unification.out = unify_archetypes(G, S_r, prune.out$C_stacked,
// prune.out$H_stacked) ' cell.clusters = unification.out$sample_assignments '
// S.norm = renormalize_input_matrix(S, cell.clusters) '	logPvals.list =
// compute_archetype_feature_specificity(S.norm, unification.out$H_unified) '
// specificity.scores = logPvals.list$upper_significance
// [[Rcpp::export]]
List compute_archetype_feature_specificity(sp_mat &S, mat &H, int thread_no = 0)
{
  field<mat> res = ACTIONet::compute_feature_specificity(S, H, thread_no);

  List out_list;
  out_list["archetypes"] = res(0);
  out_list["upper_significance"] = res(1);
  out_list["lower_significance"] = res(2);

  return (out_list);
}

//' Compute feature specificity (from archetype footprints)
//'
//' @param S Input matrix ("matrix" type)
//' @param H A soft membership matrix - Typically H_unified from the
// unify_archetypes() function.
//'
//' @return A list with the over/under-logPvals
//'
//' @examples
//' prune.out = prune_archetypes(ACTION.out$C, ACTION.out$H)
//'	G = buildNetwork(prune.out$H_stacked)
//' unification.out = unify_archetypes(G, S_r, prune.out$C_stacked,
// prune.out$H_stacked) ' cell.clusters = unification.out$sample_assignments '
// S.norm = renormalize_input_matrix(S, cell.clusters) '	logPvals.list =
// compute_archetype_feature_specificity(S.norm, unification.out$H_unified) '
// specificity.scores = logPvals.list$upper_significance
// [[Rcpp::export]]
List compute_archetype_feature_specificity_full(mat &S, mat &H, int thread_no = 0)
{
  field<mat> res = ACTIONet::compute_feature_specificity(S, H, thread_no);

  List out_list;
  out_list["archetypes"] = res(0);
  out_list["upper_significance"] = res(1);
  out_list["lower_significance"] = res(2);

  return (out_list);
}

//' Compute feature specificity (from cluster assignments)
//'
//' @param S Input matrix ("sparseMatrix")
//' @param sample_assignments Vector of cluster assignments
//'
//' @return A list with the over/under-logPvals
//'
//' @examples
//' prune.out = prune_archetypes(ACTION.out$C, ACTION.out$H)
//'	G = buildNetwork(prune.out$H_stacked)
//' unification.out = unify_archetypes(G, S_r, prune.out$C_stacked,
// prune.out$H_stacked) ' cell.clusters = unification.out$sample_assignments '
// S.norm = renormalize_input_matrix(S, cell.clusters) '	logPvals.list =
// compute_cluster_feature_specificity(S.norm, cell.clusters) '
// specificity.scores = logPvals.list$upper_significance
// [[Rcpp::export]]
List compute_cluster_feature_specificity(sp_mat &S, uvec sample_assignments, int thread_no = 0)
{
  field<mat> res = ACTIONet::compute_feature_specificity(S, sample_assignments, thread_no);

  List out_list;
  out_list["average_profile"] = res(0);
  out_list["upper_significance"] = res(1);
  out_list["lower_significance"] = res(2);

  return (out_list);
}

//' Compute feature specificity (from cluster assignments)
//'
//' @param S Input matrix ("matrix")
//' @param sample_assignments Vector of cluster assignments
//'
//' @return A list with the over/under-logPvals
//'
//' @examples
//' prune.out = prune_archetypes(ACTION.out$C, ACTION.out$H)
//'	G = buildNetwork(prune.out$H_stacked)
//' unification.out = unify_archetypes(G, S_r, prune.out$C_stacked,
// prune.out$H_stacked) ' cell.clusters = unification.out$sample_assignments '
// S.norm = renormalize_input_matrix(S, cell.clusters) '	logPvals.list =
// compute_cluster_feature_specificity(S.norm, cell.clusters) '
// specificity.scores = logPvals.list$upper_significance
// [[Rcpp::export]]
List compute_cluster_feature_specificity_full(mat &S, uvec sample_assignments, int thread_no = 0)
{
  field<mat> res = ACTIONet::compute_feature_specificity(S, sample_assignments, thread_no);

  List out_list;
  out_list["average_profile"] = res(0);
  out_list["upper_significance"] = res(1);
  out_list["lower_significance"] = res(2);

  return (out_list);
}

//' Compute coreness of graph vertices
//'
//' @param G Input graph
//'
//' @return cn core-number of each graph node
//'
//' @examples
//' G = colNets(ace)$ACTIONet
//' cn = compute_core_number(G)
// [[Rcpp::export]]
uvec compute_core_number(sp_mat &G)
{
  uvec core_num = ACTIONet::compute_core_number(G);

  return (core_num);
}

//' Compute coreness of subgraph vertices induced by each archetype
//'
//' @param G Input graph
//' @param sample_assignments Archetype discretization (output of
// unify_archetypes())
//'
//' @return cn core-number of each graph node
//'
//' @examples
//' G = colNets(ace)$ACTIONet
//' assignments = ace$archetype.assignment
//' connectivity = compute_core_number(G, assignments)
// [[Rcpp::export]]
vec compute_archetype_core_centrality(sp_mat &G, uvec sample_assignments)
{
  vec conn = ACTIONet::compute_archetype_core_centrality(G, sample_assignments);

  return (conn);
}

//' Computes network diffusion over a given network, starting with an arbitrarty
// set of initial scores
//'
//' @param G Input graph
//' @param X0 Matrix of initial values per diffusion (ncol(G) == nrow(G) ==
// ncol(X0)) ' @param thread_no Number of parallel threads (default=0) ' @param
// alpha Random-walk depth ( between [0, 1] ) ' @param max_it PageRank
// iterations
//'
//' @return Matrix of diffusion scores
//'
//' @examples
//' G = colNets(ace)$ACTIONet
//' gene.expression = Matrix::t(logcounts(ace))[c("CD19", "CD14", "CD16"), ]
//' smoothed.expression = compute_network_diffusion(G, gene.expression)
// [[Rcpp::export]]
mat compute_network_diffusion_fast(sp_mat &G, sp_mat &X0, int thread_no = 0,
                                   double alpha = 0.85, int max_it = 3)
{
  mat Diff =
      ACTIONet::compute_network_diffusion_fast(G, X0, thread_no, alpha, max_it);

  return (Diff);
}

/*
mat compute_network_diffusion_SFMULT(sp_mat &G, sp_mat &X0, double alpha = 0.85, int max_it = 3)
{
  mat Diff =
      ACTIONet::compute_network_diffusion_SFMULT(G, X0, alpha, max_it);

  return (Diff);
}
*/

//' Computes network diffusion over a given network, starting with an arbitrarty
// set of initial scores (direct approach)
//'
//' @param G Input graph
//' @param X0 Matrix of initial values per diffusion (ncol(G) == nrow(G) ==
// ncol(X0)) ' @param thread_no Number of parallel threads (default=0) ' @param
// alpha Random-walk depth ( between [0, 1] ) ' @param max_it PageRank
// iterations
//'
//' @return Matrix of diffusion scores
//'
//' @examples
//' G = colNets(ace)$ACTIONet
//' gene.expression = Matrix::t(logcounts(ace))[c("CD19", "CD14", "CD16"), ]
//' smoothed.expression = compute_network_diffusion_direct(G, gene.expression)
// [[Rcpp::export]]
mat compute_network_diffusion_direct(sp_mat &G, sp_mat &X0, int thread_no = 0,
                                     double alpha = 0.85)
{
  mat Diff =
      ACTIONet::compute_network_diffusion_direct(G, X0, thread_no, alpha);

  return (Diff);
}

//' Computes feature enrichment wrt a given annotation
//'
//' @param scores Specificity scores of features
//' @param associations Binary matrix of annotations
//' @param L Length of the top-ranked scores to scan
//'
//' @return Matrix of log-pvalues
//'
//' @examples
//' data("gProfilerDB_human")
//' G = colNets(ace)$ACTIONet
//' associations = gProfilerDB_human$SYMBOL$REAC
//' common.genes = intersect(rownames(ace), rownames(associations))
//' specificity_scores = rowFactors(ace)[["H_unified_upper_significance"]]
//' logPvals = compute_feature_specificity(specificity_scores[common.genes, ],
// annotations[common.genes, ]) ' rownames(logPvals) =
// colnames(specificity_scores) ' colnames(logPvals) = colnames(annotations)
// [[Rcpp::export]]
List assess_enrichment(mat &scores, sp_mat &associations, int thread_no = 0)
{
  field<mat> res = ACTIONet::assess_enrichment(scores, associations, thread_no);

  List out_list;
  out_list["logPvals"] = res(0);
  out_list["thresholds"] = res(1);

  return (out_list);
}

//' Computes disjoint clusters for vertices of G.
//' (It uses an adjusted DBSCAN procedure)
//'
//' @param G Adjacency matrix of the input graph
//' @param minPts, eps DBSCAN parameters
//' @param alpha Diffusion parameter for initial node ordering
//'
//' @return Matrix of log-pvalues
//'
//' @examples
//' G = colNets(ace)$ACTIONet
//' clusters = NetDBSCAN(G)
// [[Rcpp::export]]
vec NetDBSCAN(SEXP G, int minPts = 10, double eps = 0.5, double alpha = 0.85)
{
  sp_mat Adj;
  if (Rf_isS4(G))
  {
    Adj = as<arma::sp_mat>(G);
  }
  else
  {
    Adj = sp_mat(as<arma::mat>(G));
  }

  vec clusters = ACTIONet::NetDBSCAN(Adj, minPts, eps, alpha);

  return (clusters);
}

//' Clusters data points using the hierarchical DBSCAN algorithm.
//'
//' @param X Input data matrix with each row being a data point
//'
//' @return A list with \itemize{
//' \item labels
//' \item membershipProbabilities
//' \item outlierScores
//'}
//'
//' @examples
//' S_r = t(reducedDims(ace)[["S_r"]])
//' W_r = S_r %*% trace$pruning.out$C_stacked
//' X = Matrix::t(W_r)
//' HDBSCAN.out = run_HDBSCAN(X)
//' clusters = HDBSCAN.out$labels
// [[Rcpp::export]]
List run_HDBSCAN(mat &X, int minPoints = 5, int minClusterSize = 5)
{
  field<vec> res = ACTIONet::run_HDBSCAN(X, minPoints, minClusterSize);

  List out_list;
  out_list["labels"] = res(0);
  out_list["membershipProbabilities"] = res(1);
  out_list["outlierScores"] = res(2);

  return (out_list);
}

//' Computes the maximum-weight bipartite graph matching
//'
//' @param G Adjacency matrix of the input graph
//'
//' @return G_matched An adjacency matrix with a maximum of one nonzero entry on
// rows/columns
//'
//' @examples
//' G_matched = MWM_hungarian(G)
// [[Rcpp::export]]
mat MWM_hungarian(mat &G)
{
  mat G_matched = ACTIONet::MWM_hungarian(G);

  return G_matched;
}

//' Computes graph clustering using Leiden algorith over signed graphs
//'
//' @param G Adjacency matrix of the input graph
//' @param resolution_parameter Granularity of clustering. Larger values result
// in more clusters (default = 1.0) ' @param initial_clusters_ Initialization
// vector for clusters (if available) ' @param seed Random seed
//'
//' @return clusters Assignment vector of samples to clusters
//'
//' @examples
//' clusters = signed_cluster(G_signed)
// [[Rcpp::export]]
vec signed_cluster(sp_mat A, double resolution_parameter = 1.0,
                   Nullable<IntegerVector> initial_clusters_ = R_NilValue,
                   int seed = 0)
{
  set_seed(seed);

  uvec initial_clusters_uvec(A.n_rows);
  if (initial_clusters_.isNotNull())
  {
    NumericVector initial_clusters(initial_clusters_);

    for (int i = 0; i < A.n_rows; i++)
      initial_clusters_uvec(i) = initial_clusters(i);
  }
  else
  {
    for (int i = 0; i < A.n_rows; i++)
      initial_clusters_uvec(i) = i;
  }

  vec clusters = ACTIONet::signed_cluster(A, resolution_parameter,
                                          initial_clusters_uvec, seed);

  return clusters;
}

// [[Rcpp::export]]
mat unsigned_cluster_batch(
    sp_mat A, vec resolutions,
    Nullable<IntegerVector> initial_clusters_ = R_NilValue, int seed = 0)
{
  set_seed(seed);

  uvec initial_clusters_uvec(A.n_rows);
  if (initial_clusters_.isNotNull())
  {
    NumericVector initial_clusters(initial_clusters_);

    for (int i = 0; i < A.n_rows; i++)
      initial_clusters_uvec(i) = initial_clusters(i);
  }
  else
  {
    for (int i = 0; i < A.n_rows; i++)
      initial_clusters_uvec(i) = i;
  }

  mat clusters = ACTIONet::unsigned_cluster_batch(A, resolutions,
                                                  initial_clusters_uvec, seed);

  return clusters;
}

//' Computes graph clustering using Leiden algorith over unsigned graphs
//'
//' @param G Adjacency matrix of the input graph
//' @param resolution_parameter Granularity of clustering. Larger values result
// in more clusters (default = 1.0) ' @param initial_clusters_ Initialization
// vector for clusters (if available) ' @param seed Random seed
//'
//' @return clusters Assignment vector of samples to clusters
//'
//' @examples
//' clusters = unsigned_cluster(G)
// [[Rcpp::export]]
vec unsigned_cluster(sp_mat A, double resolution_parameter = 1.0,
                     Nullable<IntegerVector> initial_clusters_ = R_NilValue,
                     int seed = 0)
{
  set_seed(seed);

  uvec initial_clusters_uvec(A.n_rows);
  if (initial_clusters_.isNotNull())
  {
    NumericVector initial_clusters(initial_clusters_);

    for (int i = 0; i < A.n_rows; i++)
      initial_clusters_uvec(i) = initial_clusters(i);
  }
  else
  {
    for (int i = 0; i < A.n_rows; i++)
      initial_clusters_uvec(i) = i;
  }

  vec clusters = ACTIONet::unsigned_cluster(A, resolution_parameter,
                                            initial_clusters_uvec, seed);

  return clusters;
}

// [[Rcpp::export]]
mat sgd2_layout_weighted(sp_mat &G, mat S_r, int t_max = 30, double eps = .01,
                         int seed = 0)
{
  int n = S_r.n_cols;
  G.diag().zeros();

  int m = G.n_nonzero;
  int *I = new int[m];
  int *J = new int[m];
  double *V = new double[m];

  sp_mat::const_iterator it = G.begin();
  sp_mat::const_iterator it_end = G.end();
  int idx = 0;
  for (; it != it_end; ++it)
  {
    I[idx] = it.row();
    J[idx] = it.col();
    V[idx] = (*it);
    idx++;
  }

  mat X(2, n);
  X = S_r.rows(0, 1);
  layout_weighted(n, X.memptr(), m, I, J, V, t_max, eps, seed);

  delete[] I;
  delete[] J;
  delete[] V;

  return (trans(X));
}

// [[Rcpp::export]]
mat sgd2_layout_weighted_convergent(sp_mat &G, mat S_r, int t_max = 30,
                                    double eps = 0.01, double delta = 0.03,
                                    int t_maxmax = 200, int seed = 0)
{
  int n = S_r.n_cols;
  G.diag().zeros();

  int m = G.n_nonzero;
  int *I = new int[m];
  int *J = new int[m];
  double *V = new double[m];

  sp_mat::const_iterator it = G.begin();
  sp_mat::const_iterator it_end = G.end();
  int idx = 0;
  for (; it != it_end; ++it)
  {
    I[idx] = it.row();
    J[idx] = it.col();
    V[idx] = (*it);
    idx++;
  }

  mat X(2, n);
  X = S_r.rows(0, 1);
  layout_weighted_convergent(n, X.memptr(), m, I, J, V, t_max, eps, delta,
                             t_maxmax, seed);

  delete[] I;
  delete[] J;
  delete[] V;

  return (trans(X));
}

// [[Rcpp::export]]
mat sgd2_layout_sparse_weighted(sp_mat &G, mat S_r, int p = 200, int t_max = 30,
                                double eps = 0.01, int seed = 0)
{
  int n = S_r.n_cols;
  G.diag().zeros();

  int m = G.n_nonzero;
  int *I = new int[m];
  int *J = new int[m];
  double *V = new double[m];

  sp_mat::const_iterator it = G.begin();
  sp_mat::const_iterator it_end = G.end();
  int idx = 0;
  for (; it != it_end; ++it)
  {
    I[idx] = it.row();
    J[idx] = it.col();
    V[idx] = (*it);
    idx++;
  }

  mat X(2, n);
  X = S_r.rows(0, 1);
  layout_sparse_weighted(n, X.memptr(), m, I, J, V, p, t_max, eps, seed);

  delete[] I;
  delete[] J;
  delete[] V;

  return (trans(X));
}

//' Computes a coreset for archetypal analysis
//' Ref: Coresets for Archetypal Analysis
//(http://papers.neurips.cc/paper/8945-coresets-for-archetypal-analysis)
//'
//' @param S Input matrix (e.g., gene x cell)
//' @param m Number of samples (or 0, to be automatically identified)
//' @param seed Random seed
//'
//' @return clusters Assignment vector of samples to clusters
//'
//' @examples
//' coreset = compute_AA_coreset(S, 1000)
// [[Rcpp::export]]
List compute_AA_coreset(sp_mat &S, int m = 0)
{
  ACTIONet::Coreset coreset = ACTIONet::compute_AA_coreset(S, m);

  List out_list;
  out_list["S_coreset"] = coreset.S_coreset;
  out_list["w_coreset"] = coreset.w_coreset;

  uvec index = coreset.index + 1;
  out_list["index"] = index;

  return (out_list);
}

//' Computes reduced kernel matrix for a given (single-cell) profile and prior
// SVD
//'
//' @param S Input matrix ("sparseMatrix")
//' @param U Left singular vectors
//' @param s signular values
//' @param V Right singular vectors
//'
//' @return A named list with S_r, V, lambda, and exp_var. \itemize{
//' \item S_r: reduced kernel matrix of size reduced_dim x #samples.
//' \item V: Associated left singular-vectors (useful for reconstructing
// discriminative scores for features, such as genes). ' \item lambda, exp_var:
// Summary statistics of the sigular-values. ' }
//'
//' @examples
//' S = logcounts(sce)
//' irlba.out = irlba::irlba(S, nv = 50)
//' red.out = SVD2ACTIONred_full(S, irlba.out$u, as.matrix(irlba.out$d),
// irlba.out$v) ' Sr = red.out$S_r
// [[Rcpp::export]]
List SVD2ACTIONred(sp_mat &S, mat u, vec d, mat v)
{
  if (1 < d.n_cols)
    d = d.diag();

  field<mat> SVD_results(3);
  SVD_results(0) = u;
  SVD_results(1) = d;
  SVD_results(2) = v;

  field<mat> reduction = ACTIONet::SVD2ACTIONred(S, SVD_results);

  List res;
  res["V"] = reduction(0);

  vec sigma = reduction(1).col(0);
  res["sigma"] = sigma;

  mat V = reduction(2);
  for (int i = 0; i < V.n_cols; i++)
  {
    V.col(i) *= sigma(i);
  }
  res["S_r"] = trans(V);

  res["A"] = reduction(3);
  res["B"] = reduction(4);

  return res;
}

//' Computes reduced kernel matrix for a given (single-cell) profile and prior
// SVD
//'
//' @param S Input matrix ("sparseMatrix")
//' @param U Left singular vectors
//' @param s signular values
//' @param V Right singular vectors
//'
//' @return A named list with S_r, V, lambda, and exp_var. \itemize{
//' \item S_r: reduced kernel matrix of size reduced_dim x #samples.
//' \item V: Associated left singular-vectors (useful for reconstructing
// discriminative scores for features, such as genes). ' \item lambda, exp_var:
// Summary statistics of the sigular-values. ' }
//'
//' @examples
//' S = logcounts(sce)
//' irlba.out = irlba::irlba(S, nv = 50)
//' red.out = SVD2ACTIONred_full(S, irlba.out$u, as.matrix(irlba.out$d),
// irlba.out$v) ' Sr = red.out$S_r
// [[Rcpp::export]]
List SVD2ACTIONred_full(mat &S, mat u, vec d, mat v)
{
  if (1 < d.n_cols)
    d = d.diag();

  field<mat> SVD_results(3);
  SVD_results(0) = u;
  SVD_results(1) = d;
  SVD_results(2) = v;

  field<mat> reduction = ACTIONet::SVD2ACTIONred(S, SVD_results);

  List res;
  res["V"] = reduction(0);

  vec sigma = reduction(1).col(0);
  res["sigma"] = sigma;

  mat V = reduction(2);
  for (int i = 0; i < V.n_cols; i++)
  {
    V.col(i) *= sigma(i);
  }
  res["S_r"] = trans(V);

  res["A"] = reduction(3);
  res["B"] = reduction(4);

  return res;
}

//' Computes reduced kernel matrix for a given (single-cell) profile and prior
// SVD
//'
//' @param S Input matrix ("sparseMatrix")
//' @param U Left singular vectors
//' @param s signular values
//' @param V Right singular vectors
//'
//' @return A named list with S_r, V, lambda, and exp_var. \itemize{
//' \item S_r: reduced kernel matrix of size reduced_dim x #samples.
//' \item V: Associated left singular-vectors (useful for reconstructing
// discriminative scores for features, such as genes). ' \item lambda, exp_var:
// Summary statistics of the sigular-values. ' }
//'
//' @examples
//' S = logcounts(sce)
//' irlba.out = irlba::prcomp_irlba(S, n = 50, retx = TRUE, center = T)
//' red.out = PCA2ACTIONred_full(S, irlba.out$x, irlba.out$rotation,
// as.matrix(irlba.out$sdev)) ' Sr = red.out$S_r
// [[Rcpp::export]]
List PCA2ACTIONred(sp_mat &S, mat x, vec sdev, mat rotation)
{
  field<mat> SVD_results(3);

  vec d = sdev * sqrt(x.n_rows - 1);
  mat U = x;
  for (int i = 0; i < U.n_cols; i++)
  {
    U.col(i) /= d(i);
  }

  SVD_results(0) = U;
  SVD_results(1) = d;
  SVD_results(2) = rotation;

  field<mat> reduction = ACTIONet::PCA2ACTIONred(S, SVD_results);

  List res;
  res["V"] = reduction(0);

  vec sigma = reduction(1).col(0);
  res["sigma"] = sigma;

  mat V = reduction(2);
  for (int i = 0; i < V.n_cols; i++)
  {
    V.col(i) *= sigma(i);
  }
  res["S_r"] = trans(V);

  res["A"] = reduction(3);
  res["B"] = reduction(4);

  return res;
}

//' Computes reduced kernel matrix for a given (single-cell) profile and prior
// SVD
//'
//' @param S Input matrix ("sparseMatrix")
//' @param U Left singular vectors
//' @param s signular values
//' @param V Right singular vectors
//'
//' @return A named list with S_r, V, lambda, and exp_var. \itemize{
//' \item S_r: reduced kernel matrix of size reduced_dim x #samples.
//' \item V: Associated left singular-vectors (useful for reconstructing
// discriminative scores for features, such as genes). ' \item lambda, exp_var:
// Summary statistics of the sigular-values. ' }
//'
//' @examples
//' S = logcounts(sce)
//' irlba.out = irlba::prcomp_irlba(S, n = 50, retx = TRUE, center = T)
//' red.out = PCA2ACTIONred_full(S, irlba.out$x, irlba.out$rotation,
// as.matrix(irlba.out$sdev)) ' Sr = red.out$S_r
// [[Rcpp::export]]
List PCA2ACTIONred_full(mat &S, mat x, vec sdev, mat rotation)
{
  field<mat> SVD_results(3);

  vec d = sdev * sqrt(x.n_rows - 1);
  mat U = x;
  for (int i = 0; i < U.n_cols; i++)
  {
    U.col(i) /= d(i);
  }

  SVD_results(0) = U;
  SVD_results(1) = d;
  SVD_results(2) = rotation;

  field<mat> reduction = ACTIONet::PCA2ACTIONred(S, SVD_results);

  List res;
  res["V"] = reduction(0);

  vec sigma = reduction(1).col(0);
  res["sigma"] = sigma;

  mat V = reduction(2);
  for (int i = 0; i < V.n_cols; i++)
  {
    V.col(i) *= sigma(i);
  }
  res["S_r"] = trans(V);

  res["A"] = reduction(3);
  res["B"] = reduction(4);

  return res;
}

// [[Rcpp::export]]
List PCA2SVD(sp_mat &S, mat x, vec sdev, mat rotation)
{
  field<mat> PCA_results(3);
  vec d = sdev * sqrt(x.n_rows - 1);

  mat U = x;
  for (int i = 0; i < U.n_cols; i++)
  {
    U.col(i) /= d(i);
  }
  PCA_results(0) = U;
  PCA_results(1) = d;
  PCA_results(2) = rotation;

  field<mat> SVD_results = ACTIONet::PCA2SVD(S, PCA_results);

  List res;
  res["u"] = SVD_results(0);
  res["d"] = SVD_results(1);
  res["v"] = SVD_results(2);

  return res;
}

// [[Rcpp::export]]
List PCA2SVD_full(mat &S, mat x, vec sdev, mat rotation)
{
  field<mat> PCA_results(3);
  vec d = sdev * sqrt(x.n_rows - 1);

  mat U = x;
  for (int i = 0; i < U.n_cols; i++)
  {
    U.col(i) /= d(i);
  }
  PCA_results(0) = U;
  PCA_results(1) = d;
  PCA_results(2) = rotation;

  field<mat> SVD_results = ACTIONet::PCA2SVD(S, PCA_results);

  List res;
  res["u"] = SVD_results(0);
  res["d"] = SVD_results(1);
  res["v"] = SVD_results(2);

  return res;
}

// [[Rcpp::export]]
List SVD2PCA(sp_mat &S, mat u, vec d, mat v)
{
  if (1 < d.n_cols)
    d = d.diag();

  field<mat> SVD_results(3);
  SVD_results(0) = u;
  SVD_results(1) = d;
  SVD_results(2) = v;

  field<mat> PCA_results = ACTIONet::SVD2PCA(S, SVD_results);

  List res;
  vec s = PCA_results(1).col(0);

  mat X = PCA_results(0);
  for (int i = 0; i < X.n_cols; i++)
  {
    X.col(i) *= s(i);
  }
  res["x"] = X;
  res["rotation"] = PCA_results(2);
  res["sdev"] = s / sqrt(X.n_rows - 1);

  return res;
}

// [[Rcpp::export]]
List SVD2PCA_full(mat &S, mat u, vec d, mat v)
{
  if (1 < d.n_cols)
    d = d.diag();

  field<mat> SVD_results(3);
  SVD_results(0) = u;
  SVD_results(1) = d;
  SVD_results(2) = v;

  field<mat> PCA_results = ACTIONet::SVD2PCA(S, SVD_results);

  List res;
  vec s = PCA_results(1).col(0);

  mat X = PCA_results(0);
  for (int i = 0; i < X.n_cols; i++)
  {
    X.col(i) *= s(i);
  }
  res["x"] = X;
  res["rotation"] = PCA_results(2);
  res["sdev"] = s / sqrt(X.n_rows - 1);

  return res;
}

// [[Rcpp::export]]
List perturbedSVD(mat u, vec d, mat v, mat A, mat B)
{
  if (1 < d.n_cols)
    d = d.diag();

  field<mat> SVD_results(3);
  SVD_results(0) = u;
  SVD_results(1) = d;
  SVD_results(2) = v;

  field<mat> perturbed_SVD = ACTIONet::perturbedSVD(SVD_results, A, B);

  List res;
  res["u"] = perturbed_SVD(0);
  res["d"] = perturbed_SVD(1).col(0);
  res["v"] = perturbed_SVD(2);

  return res;
}

// [[Rcpp::export]]
mat computeFullSim(mat &H, int thread_no = 0)
{
  mat G = ACTIONet::computeFullSim(H, thread_no);

  return (G);
}

// [[Rcpp::export]]
List run_subACTION(mat &S_r, mat &W_parent, mat &H_parent, int kk, int k_min,
                   int k_max, int thread_no, int max_it = 50,
                   double min_delta = 1e-16)
{
  ACTIONet::ACTION_results trace =
      ACTIONet::run_subACTION(S_r, W_parent, H_parent, kk - 1, k_min, k_max,
                              thread_no, max_it, min_delta);

  List res;

  List C(k_max);
  for (int i = k_min; i <= k_max; i++)
  {
    C[i - 1] = trace.C[i];
  }
  res["C"] = C;

  List H(k_max);
  for (int i = k_min; i <= k_max; i++)
  {
    H[i - 1] = trace.H[i];
  }
  res["H"] = H;

  return res;
}

// [[Rcpp::export]]
List deflate_reduction(mat &old_S_r, mat &old_V, mat &old_A, mat &old_B,
                       vec &old_sigma, mat &A, mat &B)
{
  field<mat> SVD_results(5);

  SVD_results(0) = old_V;
  SVD_results(1) = old_sigma;
  SVD_results(2) = old_S_r;
  for (int i = 0; i < old_sigma.n_elem; i++)
  {
    SVD_results(2).col(i) /= old_sigma(i);
  }
  SVD_results(3) = old_A;
  SVD_results(4) = old_B;

  field<mat> deflated_reduction =
      ACTIONet::deflate_reduction(SVD_results, A, B);

  List res;
  res["V"] = deflated_reduction(0);

  vec sigma = deflated_reduction(1).col(0);
  res["sigma"] = sigma;

  mat V = deflated_reduction(2);
  for (int i = 0; i < V.n_cols; i++)
  {
    V.col(i) *= sigma(i);
  }
  res["S_r"] = trans(V);

  res["A"] = deflated_reduction(3);
  res["B"] = deflated_reduction(4);

  return res;
}

// [[Rcpp::export]]
List orthogonalize_batch_effect(sp_mat &S, mat &old_S_r, mat &old_V, mat &old_A,
                                mat &old_B, vec &old_sigma, mat &design)
{
  field<mat> SVD_results(5);

  SVD_results(0) = old_V;
  SVD_results(1) = old_sigma;
  SVD_results(2) = old_S_r;
  for (int i = 0; i < old_sigma.n_elem; i++)
  {
    SVD_results(2).col(i) /= old_sigma(i);
  }
  SVD_results(3) = old_A;
  SVD_results(4) = old_B;

  field<mat> orthogonalized_reduction =
      ACTIONet::orthogonalize_batch_effect(S, SVD_results, design);

  List res;
  res["V"] = orthogonalized_reduction(0);

  vec sigma = orthogonalized_reduction(1).col(0);
  res["sigma"] = sigma;

  mat V = orthogonalized_reduction(2);
  for (int i = 0; i < V.n_cols; i++)
  {
    V.col(i) *= sigma(i);
  }
  res["S_r"] = trans(V);

  res["A"] = orthogonalized_reduction(3);
  res["B"] = orthogonalized_reduction(4);

  return res;
}

//[[Rcpp::export]]
List orthogonalize_batch_effect_full(mat &S, mat &old_S_r, mat &old_V,
                                     mat &old_A, mat &old_B, vec &old_sigma,
                                     mat &design)
{
  field<mat> SVD_results(5);

  SVD_results(0) = old_V;
  SVD_results(1) = old_sigma;
  SVD_results(2) = old_S_r;
  for (int i = 0; i < old_sigma.n_elem; i++)
  {
    SVD_results(2).col(i) /= old_sigma(i);
  }
  SVD_results(3) = old_A;
  SVD_results(4) = old_B;

  field<mat> orthogonalized_reduction =
      ACTIONet::orthogonalize_batch_effect(S, SVD_results, design);

  List res;
  res["V"] = orthogonalized_reduction(0);

  vec sigma = orthogonalized_reduction(1).col(0);
  res["sigma"] = sigma;

  mat V = orthogonalized_reduction(2);
  for (int i = 0; i < V.n_cols; i++)
  {
    V.col(i) *= sigma(i);
  }
  res["S_r"] = trans(V);

  res["A"] = orthogonalized_reduction(3);
  res["B"] = orthogonalized_reduction(4);

  return res;
}

// [[Rcpp::export]]
List orthogonalize_basal(sp_mat &S, mat &old_S_r, mat &old_V, mat &old_A,
                         mat &old_B, vec &old_sigma, mat &basal)
{
  field<mat> SVD_results(5);

  SVD_results(0) = old_V;
  SVD_results(1) = old_sigma;
  SVD_results(2) = old_S_r;
  for (int i = 0; i < old_sigma.n_elem; i++)
  {
    SVD_results(2).col(i) /= old_sigma(i);
  }
  SVD_results(3) = old_A;
  SVD_results(4) = old_B;

  field<mat> orthogonalized_reduction =
      ACTIONet::orthogonalize_basal(S, SVD_results, basal);

  List res;
  res["V"] = orthogonalized_reduction(0);

  vec sigma = orthogonalized_reduction(1).col(0);
  res["sigma"] = sigma;

  mat V = orthogonalized_reduction(2);
  for (int i = 0; i < V.n_cols; i++)
  {
    V.col(i) *= sigma(i);
  }
  res["S_r"] = trans(V);

  res["A"] = orthogonalized_reduction(3);
  res["B"] = orthogonalized_reduction(4);

  return res;
}

//[[Rcpp::export]]
List orthogonalize_basal_full(mat &S, mat &old_S_r, mat &old_V,
                              mat &old_A, mat &old_B, vec &old_sigma,
                              mat &basal)
{
  field<mat> SVD_results(5);

  SVD_results(0) = old_V;
  SVD_results(1) = old_sigma;
  SVD_results(2) = old_S_r;
  for (int i = 0; i < old_sigma.n_elem; i++)
  {
    SVD_results(2).col(i) /= old_sigma(i);
  }
  SVD_results(3) = old_A;
  SVD_results(4) = old_B;

  field<mat> orthogonalized_reduction =
      ACTIONet::orthogonalize_basal(S, SVD_results, basal);

  List res;
  res["V"] = orthogonalized_reduction(0);

  vec sigma = orthogonalized_reduction(1).col(0);
  res["sigma"] = sigma;

  mat V = orthogonalized_reduction(2);
  for (int i = 0; i < V.n_cols; i++)
  {
    V.col(i) *= sigma(i);
  }
  res["S_r"] = trans(V);

  res["A"] = orthogonalized_reduction(3);
  res["B"] = orthogonalized_reduction(4);

  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
umat MWM_rank1(vec u, vec v, double u_threshold = 0, double v_threshold = 0)
{
  umat pairs = ACTIONet::MWM_rank1(u, v, u_threshold, v_threshold);

  pairs = pairs + 1;

  return (pairs);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat NetEnh(mat A)
{
  mat A_enh = ACTIONet::NetEnh(A);

  return (A_enh);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::NumericVector run_LPA(sp_mat &G, vec labels, double lambda = 1, int iters = 3, double sig_threshold = 3, Nullable<IntegerVector> fixed_labels_ = R_NilValue, int thread_no = 0)
{
  uvec fixed_labels_vec;
  if (fixed_labels_.isNotNull())
  {
    NumericVector fixed_labels(fixed_labels_);
    fixed_labels_vec.set_size(fixed_labels.size());
    for (int i = 0; i < fixed_labels.size(); i++)
    {
      fixed_labels_vec(i) = fixed_labels(i) - 1;
    }
  }

  mat new_labels = ACTIONet::LPA(G, labels, lambda, iters, sig_threshold, fixed_labels_vec, thread_no);
  // vec labels_out = arma::conv_to<arma::vec>::from(new_labels);

  return (arma2vec(new_labels));
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List run_AA_with_batch_correction(mat &Z, mat &W0, vec batch, int max_it = 100, int max_correction_rounds = 10, double lambda = 1, double min_delta = 1e-6)
{

  field<mat> res = ACTIONet::run_AA_with_batch_correction(Z, W0, batch, max_it, max_correction_rounds, lambda, min_delta);

  List out;
  out["C"] = res(0);
  out["H"] = res(1);
  out["Z_cor"] = res(2);
  out["W"] = res(2) * res(0);

  return (out);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List run_ACTION_with_batch_correction(mat &S_r, vec batch, int k_min, int k_max, int thread_no,
                                      int max_it = 100, int max_correction_rounds = 10, double lambda = 1, double min_delta = 1e-6)
{

  ACTIONet::ACTION_results trace = ACTIONet::run_ACTION_with_batch_correction(S_r, batch, k_min, k_max, thread_no, max_it, max_correction_rounds, lambda, min_delta);

  List res;

  List C(k_max);
  for (int i = k_min; i <= k_max; i++)
  {
    mat cur_C = trace.C[i];
    C[i - 1] = cur_C;
  }
  res["C"] = C;

  List H(k_max);
  for (int i = k_min; i <= k_max; i++)
  {
    mat cur_H = trace.H[i];
    H[i - 1] = cur_H;
  }
  res["H"] = H;

  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat compute_marker_aggregate_stats(sp_mat &G, sp_mat &S, sp_mat &marker_mat, double alpha = 0.85, int max_it = 5, int thread_no = 0, bool ignore_baseline_expression = false)
{
  mat stats = ACTIONet::compute_marker_aggregate_stats(G, S, marker_mat, alpha, max_it, thread_no, ignore_baseline_expression);

  return (stats);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
sp_mat LSI(sp_mat &X, double size_factor = 100000)
{
  sp_mat TFIDF = ACTIONet::LSI(X, size_factor);

  return (TFIDF);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat compute_marker_aggregate_stats_TFIDF_sum_smoothed(sp_mat &G, sp_mat &S, sp_mat &marker_mat, double alpha = 0.85, int max_it = 5, int perm_no = 100, int thread_no = 0, int normalization = 1)
{
  mat stats = ACTIONet::compute_marker_aggregate_stats_TFIDF_sum_smoothed(G, S, marker_mat, alpha, max_it, perm_no, thread_no, normalization);

  return (stats);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List autocorrelation_Geary(sp_mat G, mat scores, int normalization_method = 1, int perm_no = 30, int thread_no = 0)
{
  field<vec> out = ACTIONet::autocorrelation_Geary(G, scores, normalization_method, perm_no, thread_no);

  List res;
  res["Geary_C"] = out[0];
  res["zscore"] = out[1];
  res["mu"] = out[2];
  res["sigma"] = out[3];

  return (res);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List autocorrelation_Geary_full(mat G, mat scores, int normalization_method = 1, int perm_no = 30, int thread_no = 0)
{
  field<vec> out = ACTIONet::autocorrelation_Geary(G, scores, normalization_method, perm_no, thread_no);

  List res;
  res["Geary_C"] = out[0];
  res["zscore"] = out[1];
  res["mu"] = out[2];
  res["sigma"] = out[3];

  return (res);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List autocorrelation_Moran(sp_mat G, mat scores, int normalization_method = 1, int perm_no = 30, int thread_no = 0)
{
  field<vec> out = ACTIONet::autocorrelation_Moran(G, scores, normalization_method, perm_no, thread_no);

  List res;
  res["Moran_I"] = out[0];
  res["zscore"] = out[1];
  res["mu"] = out[2];
  res["sigma"] = out[3];

  return (res);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List autocorrelation_Moran_full(mat G, mat scores, int normalization_method = 1, int perm_no = 30, int thread_no = 0)
{
  field<vec> out = ACTIONet::autocorrelation_Moran(G, scores, normalization_method, perm_no, thread_no);

  List res;
  res["Moran_I"] = out[0];
  res["zscore"] = out[1];
  res["mu"] = out[2];
  res["sigma"] = out[3];

  return (res);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
vec spmat_vec_product(sp_mat &A, vec &x)
{
  vec res = ACTIONet::spmat_vec_product(A, x);
  return (res);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat spmat_mat_product(sp_mat &A, mat &B)
{
  mat res = ACTIONet::spmat_mat_product(A, B);
  return (res);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
sp_mat spmat_spmat_product(sp_mat &A, sp_mat &B)
{
  sp_mat res = ACTIONet::spmat_spmat_product(A, B);

  return (res);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat spmat_mat_product_parallel(sp_mat &A, mat &B, int thread_no)
{
  mat res = ACTIONet::spmat_mat_product_parallel(A, B, thread_no);

  return (res);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat mat_mat_product_parallel(mat &A, mat &B, int thread_no)
{
  mat res = ACTIONet::mat_mat_product_parallel(A, B, thread_no);

  return (res);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List transform_layout(sp_mat &G, sp_mat &inter_graph, mat reference_coordinates, int compactness_level = 50,
                      unsigned int n_epochs = 500,
                      int layout_alg = 0, int thread_no = 0,
                      int seed = 0)
{
  field<mat> res = ACTIONet::transform_layout(G, inter_graph, reference_coordinates, compactness_level,
                                              n_epochs, layout_alg, thread_no, seed);

  List out_list;
  out_list["coordinates"] = res(0);
  out_list["coordinates_3D"] = res(1);
  out_list["colors"] = res(2);

  return (out_list);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
sp_mat normalize_adj(sp_mat &G, int norm_type = 0)
{
  sp_mat P = ACTIONet::normalize_adj(G, norm_type);

  return (P);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat compute_network_diffusion_Chebyshev(sp_mat &P, mat &X0, int thread_no = 0, double alpha = 0.85, int max_it = 5, double res_threshold = 1e-8)
{
  if (P.n_rows != X0.n_rows)
  {
    fprintf(stderr, "Dimnsion mismatch: P (%dx%d) and X0 (%dx%d)\n", P.n_rows, P.n_cols, X0.n_rows, X0.n_cols);
    return (mat());
  }
  mat X = ACTIONet::compute_network_diffusion_Chebyshev(P, X0, thread_no, alpha, max_it, res_threshold);

  return (X);
}

//' Computes network diffusion over a given network, starting with an arbitrarty
// set of initial scores
//'
//' @param G Input graph
//' @param X0 Matrix of initial values per diffusion (ncol(G) == nrow(G) ==
// ncol(X0)) ' @param thread_no Number of parallel threads (default=0) ' @param
// alpha Random-walk depth ( between [0, 1] ) ' @param max_it PageRank
// iterations
//'
//' @return Matrix of diffusion scores
//'
//' @examples
//' G = colNets(ace)$ACTIONet
//' gene.expression = Matrix::t(logcounts(ace))[c("CD19", "CD14", "CD16"), ]
//' smoothed.expression = compute_network_diffusion_approx(G, gene.expression)
// [[Rcpp::export]]
mat compute_network_diffusion_approx(sp_mat &G, mat &X0, int thread_no = 0, double alpha = 0.85, int max_it = 5, double res_threshold = 1e-8, int norm_type = 0)
{
  if (G.n_rows != X0.n_rows)
  {
    stderr_printf("Dimension mismatch: G (%dx%d) and X0 (%dx%d)\n", G.n_rows, G.n_cols, X0.n_rows, X0.n_cols);
    return (mat());
  }

  sp_mat P = normalize_adj(G, norm_type);
  mat X = ACTIONet::compute_network_diffusion_Chebyshev(P, X0, thread_no, alpha, max_it, res_threshold);

  return (X);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat compute_marker_aggregate_stats_nonparametric(mat &S, sp_mat &marker_mat, int thread_no = 0)
{
  mat X = ACTIONet::compute_marker_aggregate_stats_nonparametric(S, marker_mat, thread_no);
  return (X);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat compute_markers_eigengene(mat &S, sp_mat &marker_mat, int normalization = 0, int thread_no = 0)
{
  mat X = ACTIONet::compute_markers_eigengene(S, marker_mat, normalization, thread_no);
  return (X);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
vec sweepcut(sp_mat &A, vec s, int min_size = 5, int max_size = -1)
{
  vec cond = ACTIONet::sweepcut(A, s, min_size, max_size);

  return (cond);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat aggregate_genesets_mahalanobis_2archs(sp_mat &G, sp_mat &S, sp_mat &marker_mat, int network_normalization_method = 0, int expression_normalization_method = 0, int gene_scaling_method = 0, double pre_alpha = 0.85, double post_alpha = 0.85, int thread_no = 0) {
  mat stats = ACTIONet::aggregate_genesets_mahalanobis_2archs(G, S, marker_mat, network_normalization_method, expression_normalization_method, gene_scaling_method, pre_alpha, post_alpha, thread_no);

  return (stats);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat aggregate_genesets_mahalanobis_2gmm(sp_mat &G, sp_mat &S, sp_mat &marker_mat, int network_normalization_method = 0, int expression_normalization_method = 0, int gene_scaling_method = 0, double pre_alpha = 0.85, double post_alpha = 0.85, int thread_no = 0) {
  mat stats = ACTIONet::aggregate_genesets_mahalanobis_2gmm(G, S, marker_mat, network_normalization_method, expression_normalization_method, gene_scaling_method, pre_alpha, post_alpha, thread_no);

  return (stats);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat aggregate_genesets_weighted_enrichment(sp_mat &G, sp_mat &S, sp_mat &marker_mat, int network_normalization_method = 0, int expression_normalization_method = 0, int gene_scaling_method = 3, double pre_alpha = 0.85, double post_alpha = 0.85, int thread_no = 0) {
  mat stats = ACTIONet::aggregate_genesets_weighted_enrichment(G, S, marker_mat, network_normalization_method, expression_normalization_method, gene_scaling_method, pre_alpha, post_alpha, thread_no);

  return (stats);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat aggregate_genesets_weighted_enrichment_permutation(sp_mat &G, sp_mat &S, sp_mat &marker_mat, int network_normalization_method = 0, int expression_normalization_method = 0, int gene_scaling_method = 3, double pre_alpha = 0.85, double post_alpha = 0.85, int thread_no = 0, int perm_no = 30) {
  mat stats = ACTIONet::aggregate_genesets_weighted_enrichment_permutation(G, S, marker_mat, network_normalization_method, expression_normalization_method, gene_scaling_method, pre_alpha, post_alpha, thread_no, perm_no);

  return (stats);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List recursiveNMU(mat M, int dim = 100, int max_SVD_iter = 1000, int max_iter_inner = 100) {
  field<mat> stats = ACTIONet::recursiveNMU(M, dim, max_SVD_iter, max_iter_inner);

  List res;
  res["W"] = stats[0];
  res["H"] = trans(stats[1]);
  res["factor_weights"] = stats[2];

  return (res);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List recursiveNMU_mine(mat M, int dim = 100, int max_SVD_iter = 1000, int max_iter_inner = 100) {
  field<mat> stats = ACTIONet::recursiveNMU_mine(M, dim, max_SVD_iter, max_iter_inner);

  List res;
  res["W"] = stats[0];
  res["H"] = trans(stats[1]);
  res["factor_weights"] = stats[2];

  return (res);
}


  sp_mat smoothKNN(sp_mat D, int thread_no = -1)
  {
    double epsilon = 1e-6;

    int nV = D.n_cols;
    sp_mat G = D;

    //#pragma omp parallel for num_threads(thread_no)
    for (int i = 0; i < nV; i++)
    {
      //  ParallelFor(0, nV, thread_no, [&](size_t i, size_t threadId) {
      sp_mat v = D.col(i);
      vec vals = nonzeros(v);
      if (vals.n_elem > 0)
      {
        double rho = min(vals);
        vec negated_shifted_vals = -(vals - rho);
        double target = log2(vals.n_elem);

        // Binary search to find optimal sigma
        double sigma = 1.0;
        double lo = 0.0;
        double hi = DBL_MAX;

        int j;
        for (j = 0; j < 64; j++)
        {
          double obj = sum(exp(negated_shifted_vals / sigma));

          if (abs(obj - target) < epsilon)
          {
            break;
          }

          if (target < obj)
          {
            hi = sigma;
            sigma = 0.5 * (lo + hi);
          }
          else
          {
            lo = sigma;
            if (hi == DBL_MAX)
            {
              sigma *= 2;
            }
            else
            {
              sigma = 0.5 * (lo + hi);
            }
          }
        }

        double obj = sum(exp(negated_shifted_vals / sigma));

        for (sp_mat::col_iterator it = G.begin_col(i); it != G.end_col(i); ++it)
        {
          *it = max(1e-16, exp(-max(0.0, (*it) - rho) / sigma));
        }
      }
    }

    return (G);
  }

template <typename T>
auto lget(List list, const std::string &name, T default_value) -> T {
  auto key = name.c_str();
  if (!list.containsElementNamed(key)) {
    return default_value;
  } else {
    return list[key];
  }
}

// Template class specialization to handle different rng/batch combinations
template <bool DoBatch = true> struct BatchRngFactory {
  using PcgFactoryType = batch_pcg_factory;
  using TauFactoryType = batch_tau_factory;
};
template <> struct BatchRngFactory<false> {
  using PcgFactoryType = pcg_factory;
  using TauFactoryType = tau_factory;
};

struct UmapFactory {
  bool move_other;
  bool pcg_rand;
  std::vector<float> &head_embedding;
  std::vector<float> &tail_embedding;
  const std::vector<unsigned int> &positive_head;
  const std::vector<unsigned int> &positive_tail;
  const std::vector<unsigned int> &positive_ptr;
  unsigned int n_epochs;
  unsigned int n_head_vertices;
  unsigned int n_tail_vertices;
  const std::vector<float> &epochs_per_sample;
  float initial_alpha;
  List opt_args;
  float negative_sample_rate;
  bool batch;
  std::size_t n_threads;
  std::size_t grain_size;
  bool verbose;
  std::mt19937_64 engine;


  UmapFactory(bool move_other, bool pcg_rand,
              std::vector<float> &head_embedding,
              std::vector<float> &tail_embedding,
              const std::vector<unsigned int> &positive_head,
              const std::vector<unsigned int> &positive_tail,
              const std::vector<unsigned int> &positive_ptr,
              unsigned int n_epochs, unsigned int n_head_vertices,
              unsigned int n_tail_vertices,
              const std::vector<float> &epochs_per_sample, float initial_alpha,
              List opt_args, float negative_sample_rate, bool batch,
              std::size_t n_threads, std::size_t grain_size, bool verbose, std::mt19937_64 &engine)
      : move_other(move_other), pcg_rand(pcg_rand),
        head_embedding(head_embedding), tail_embedding(tail_embedding),
        positive_head(positive_head), positive_tail(positive_tail),
        positive_ptr(positive_ptr), n_epochs(n_epochs),
        n_head_vertices(n_head_vertices), n_tail_vertices(n_tail_vertices),
        epochs_per_sample(epochs_per_sample), initial_alpha(initial_alpha),
        opt_args(opt_args), negative_sample_rate(negative_sample_rate),
        batch(batch), n_threads(n_threads), grain_size(grain_size),
        verbose(verbose), engine(engine) {}

  template <typename Gradient> void create(const Gradient &gradient) {
    if (move_other) {
      create_impl<true>(gradient, pcg_rand, batch);
    } else {
      create_impl<false>(gradient, pcg_rand, batch);
    }
  }

  template <bool DoMove, typename Gradient>
  void create_impl(const Gradient &gradient, bool pcg_rand, bool batch) {
    if (batch) {
      create_impl<BatchRngFactory<true>, DoMove>(gradient, pcg_rand, batch);
    } else {
      create_impl<BatchRngFactory<false>, DoMove>(gradient, pcg_rand, batch);
    }
  }

  template <typename BatchRngFactory, bool DoMove, typename Gradient>
  void create_impl(const Gradient &gradient, bool pcg_rand, bool batch) {
    if (pcg_rand) {
      create_impl<typename BatchRngFactory::PcgFactoryType, DoMove>(gradient,
                                                                    batch);
    } else {
      create_impl<typename BatchRngFactory::TauFactoryType, DoMove>(gradient,
                                                                    batch);
    }
  }

  auto create_adam(List opt_args) -> uwot::Adam {
    float alpha = lget(opt_args, "alpha", 1.0);
    float beta1 = lget(opt_args, "beta1", 0.9);
    float beta2 = lget(opt_args, "beta2", 0.999);
    float eps = lget(opt_args, "eps", 1e-7);
    if (verbose) {
      Rcerr << "Optimizing with Adam"
            << " alpha = " << alpha << " beta1 = " << beta1
            << " beta2 = " << beta2 << " eps = " << eps << std::endl;
    }

    return uwot::Adam(alpha, beta1, beta2, eps, head_embedding.size());
  }

  auto create_sgd(List opt_args) -> uwot::Sgd {
    float alpha = lget(opt_args, "alpha", 1.0);
    if (verbose) {
      Rcerr << "Optimizing with SGD"
            << " alpha = " << alpha << std::endl;
    }

    return uwot::Sgd(alpha);
  }

  template <typename RandFactory, bool DoMove, typename Gradient>
  void create_impl(const Gradient &gradient, bool batch) {
    if (batch) {
      std::string opt_name = opt_args["method"];
      if (opt_name == "adam") {
        auto opt = create_adam(opt_args);
        create_impl_batch_opt<decltype(opt), RandFactory, DoMove, Gradient>(
            gradient, opt, batch);
      } else if (opt_name == "sgd") {
        auto opt = create_sgd(opt_args);
        create_impl_batch_opt<decltype(opt), RandFactory, DoMove, Gradient>(
            gradient, opt, batch);
      } else {
        stop("Unknown optimization method");
      }
    } else {
      const std::size_t ndim = head_embedding.size() / n_head_vertices;
      uwot::Sampler sampler(epochs_per_sample, negative_sample_rate);
      uwot::InPlaceUpdate<DoMove> update(head_embedding, tail_embedding,
                                         initial_alpha);
      uwot::EdgeWorker<Gradient, decltype(update), RandFactory> worker(
          gradient, update, positive_head, positive_tail, sampler, ndim,
          n_tail_vertices, n_threads, engine);
      create_impl(worker, gradient);
    }
  }

  template <typename Opt, typename RandFactory, bool DoMove, typename Gradient>
  void create_impl_batch_opt(const Gradient &gradient, Opt &opt, bool batch) {
    uwot::Sampler sampler(epochs_per_sample, negative_sample_rate);
    const std::size_t ndim = head_embedding.size() / n_head_vertices;
    uwot::BatchUpdate<DoMove, decltype(opt)> update(
        head_embedding, tail_embedding, opt);
    uwot::NodeWorker<Gradient, decltype(update), RandFactory> worker(
        gradient, update, positive_head, positive_tail, positive_ptr, sampler,
        ndim, n_tail_vertices, engine);
    create_impl(worker, gradient);
  }

  template <typename Worker, typename Gradient>
  void create_impl(Worker &worker, const Gradient &gradient) {

    if (n_threads > 0) {
      RParallel parallel(n_threads, grain_size);
      create_impl(worker, gradient, parallel);
    } else {
      RSerial serial;
      create_impl(worker, gradient, serial);
    }
  }

  template <typename Worker, typename Gradient,
            typename Parallel>
  void create_impl(Worker &worker, const Gradient &gradient,
                   Parallel &parallel) {
    uwot::optimize_layout(worker, n_epochs, parallel);
  }
};

auto r_to_coords(NumericMatrix head_embedding,
                 Nullable<NumericMatrix> tail_embedding) -> uwot::Coords {
  auto head_vec = as<std::vector<float>>(head_embedding);
  if (tail_embedding.isNull()) {
    return uwot::Coords(head_vec);
  } else {
    auto tail_vec = as<std::vector<float>>(tail_embedding);
    return uwot::Coords(head_vec, tail_vec);
  }
}

auto r_to_coords(NumericMatrix head_embedding) -> uwot::Coords {
  auto head_vec = as<std::vector<float>>(head_embedding);
  return uwot::Coords(head_vec);
}

void validate_args(List method_args,
                   const std::vector<std::string> &arg_names) {
  for (auto &arg_name : arg_names) {
    if (!method_args.containsElementNamed(arg_name.c_str())) {
      stop("Missing embedding method argument: " + arg_name);
    }
  }
}

void create_umap(UmapFactory &umap_factory, List method_args) {
  std::vector<std::string> arg_names = {"a", "b", "gamma", "approx_pow"};
  validate_args(method_args, arg_names);

  float a = method_args["a"];
  float b = method_args["b"];
  float gamma = method_args["gamma"];
  bool approx_pow = method_args["approx_pow"];
  if (approx_pow) {
    const uwot::apumap_gradient gradient(a, b, gamma);
    umap_factory.create(gradient);
  } else {
    const uwot::umap_gradient gradient(a, b, gamma);
    umap_factory.create(gradient);
  }
}

void create_tumap(UmapFactory &umap_factory, List) {
  const uwot::tumap_gradient gradient;
  umap_factory.create(gradient);
}

void create_umapai(UmapFactory &umap_factory, List method_args) {
  std::vector<std::string> arg_names = {"ai", "b", "ndim"};
  validate_args(method_args, arg_names);

  std::vector<float> ai = method_args["ai"];
  float b = method_args["b"];
  std::size_t ndim = method_args["ndim"];
  const uwot::umapai_gradient gradient(ai, b, ndim);
  umap_factory.create(gradient);
}

void create_umapai2(UmapFactory &umap_factory, List method_args) {
  std::vector<std::string> arg_names = {"ai", "aj", "b", "ndim"};
  validate_args(method_args, arg_names);

  std::vector<float> ai = method_args["ai"];
  std::vector<float> aj = method_args["ai"];
  float b = method_args["b"];
  std::size_t ndim = method_args["ndim"];
  const uwot::umapai2_gradient gradient(ai, aj, b, ndim);
  umap_factory.create(gradient);
}

void create_pacmap(UmapFactory &umap_factory, List method_args) {
  std::vector<std::string> arg_names = {"a", "b"};
  validate_args(method_args, arg_names);

  float a = method_args["a"];
  float b = method_args["b"];
  const uwot::pacmap_gradient gradient(a, b);
  umap_factory.create(gradient);
}

void create_largevis(UmapFactory &umap_factory, List method_args) {
  std::vector<std::string> arg_names = {"gamma"};
  validate_args(method_args, arg_names);

  float gamma = method_args["gamma"];
  const uwot::largevis_gradient gradient(gamma);
  umap_factory.create(gradient);
}

// [[Rcpp::export]]
List optimize_layout_interface_v2(sp_mat &G, mat &initial_position,
    unsigned int n_epochs, const std::string &method,
    List method_args, float initial_alpha, List opt_args, float negative_sample_rate,
    bool pcg_rand = true, bool batch = false, std::size_t n_threads = 0,
    std::size_t grain_size = 1, bool move_other = true, bool verbose = false, int seed = 0) {

    std::mt19937_64 engine(seed);

    mat init_coors = initial_position.rows(0, 2);

    sp_mat H = G;
    /*
    H.for_each([](sp_mat::elem_type &val)
               { val = 1 / val; });
    H = smoothKNN(H, n_threads);
  */

    sp_mat Ht = trans(H);
    Ht.sync();

    unsigned int nV = H.n_rows;
    unsigned int nE = H.n_nonzero;
    vector<unsigned int> positive_head(nE);
    vector<unsigned int> positive_tail(nE);
    vector<float> epochs_per_sample(nE);

    std::vector<unsigned int> positive_ptr(Ht.n_cols + 1);

    int i = 0;
    double w_max = max(max(H));
    if(batch == false ) {
      for (sp_mat::iterator it = H.begin(); it != H.end(); ++it)
      {
        epochs_per_sample[i] = w_max / (*it);
        positive_head[i] = it.row();
        positive_tail[i] = it.col();
        i++;
      }
    }
    else {
      for (sp_mat::iterator it = Ht.begin(); it != Ht.end(); ++it)
      {
        epochs_per_sample[i] = w_max / (*it);
        positive_tail[i] = it.row();
        positive_head[i] = it.col();
        i++;
      }      
      for (int k = 0; k < Ht.n_cols + 1; k++)
      {
          positive_ptr[k] = Ht.col_ptrs[k];
      }      
    }

    // Initial coordinates of vertices (0-simplices)
    vector<float> head_embedding(init_coors.n_cols * 2);
    fmat sub_coor = conv_to<fmat>::from(init_coors.rows(0, 1));
    float *ptr = sub_coor.memptr();
    memcpy(head_embedding.data(), ptr, sizeof(float) * head_embedding.size());
    vector<float> tail_embedding(head_embedding);

    uwot::Coords coords = uwot::Coords(head_embedding);

  //auto coords = r_to_coords(head_embedding, tail_embedding);
  const std::size_t ndim = head_embedding.size() / nV;


  UmapFactory umap_factory(move_other, pcg_rand, coords.get_head_embedding(),
                           coords.get_tail_embedding(), positive_head,
                           positive_tail, positive_ptr, n_epochs,
                           nV, nV, epochs_per_sample,
                           initial_alpha, opt_args, negative_sample_rate, batch,
                           n_threads, grain_size, verbose, engine);

/*
  UmapFactory umap_factory(move_other, pcg_rand, head_embedding,
                           tail_embedding, positive_head,
                           positive_tail, positive_ptr, n_epochs,
                           nV, nV, epochs_per_sample,
                           initial_alpha, opt_args, negative_sample_rate, batch,
                           n_threads, grain_size, verbose, engine);
*/

  if (verbose) {
    Rcerr << "Using method '" << method << "'" << std::endl;
  }
  if (method == "umap") {
    create_umap(umap_factory, method_args);
  } else if (method == "tumap") {
    create_tumap(umap_factory, method_args);
  } else if (method == "largevis") {
    create_largevis(umap_factory, method_args);
  } else if (method == "pacmap") {
    create_pacmap(umap_factory, method_args);
  } else if (method == "leopold") {
    create_umapai(umap_factory, method_args);
  } else if (method == "leopold2") {
    create_umapai2(umap_factory, method_args);
  } else {
    stop("Unknown method: '" + method + "'");
  }

/*
  

  return NumericMatrix(ndim, nV,
                       head_embedding.begin());
*/

/*
  fmat coordinates_float(coords.get_head_embedding().begin(), ndim, nV);
  mat coordinates = trans(conv_to<mat>::from(coordinates_float));
  */

  List res;
  res["coordinates"] = NumericMatrix(ndim, nV, coords.get_head_embedding().begin());
  res["positive_head"] = positive_head;
  res["positive_tail"] = positive_tail;
  res["epochs_per_sample"] = epochs_per_sample;  
  res["positive_ptr"] = positive_ptr;
  res["head_embedding"] = head_embedding;
  res["tail_embedding"] = tail_embedding;

  return (res);



}


template<typename Float>
std::pair<Float, Float> find_ab(Float spread, Float min_dist, Float grid = 300, Float limit = 0.5, int iter = 50, Float tol = 1e-6) {
    Float x_half = std::log(limit) * -spread + min_dist;
    Float d_half = limit / -spread;

    // Compute the x and y coordinates of the expected distance curve.
    std::vector<Float> grid_x(grid), grid_y(grid), log_x(grid);
    const Float delta = spread * 3 / grid;
    for (int g = 0; g < grid; ++g) {
        grid_x[g] = (g + 1) * delta; // +1 to avoid meaningless least squares result at x = 0, where both curves have y = 1 (and also the derivative w.r.t. b is not defined).
        log_x[g] = std::log(grid_x[g]);
        grid_y[g] = (grid_x[g] <= min_dist ? 1 : std::exp(- (grid_x[g] - min_dist) / spread));
    }

    // Starting estimates.
    Float b = - d_half * x_half / (1 / limit - 1) / (2 * limit * limit);
    Float a = (1 / limit - 1) / std::pow(x_half, 2 * b);

    std::vector<Float> observed_y(grid), xpow(grid);
    auto compute_ss = [&](Float A, Float B) -> Float {
        for (int g = 0; g < grid; ++g) {
            xpow[g] = std::pow(grid_x[g], 2 * B);
            observed_y[g] = 1 / (1 + A * xpow[g]);
        }

        Float ss = 0;
        for (int g = 0; g < grid; ++g) {
            ss += (grid_y[g] - observed_y[g]) * (grid_y[g] - observed_y[g]);
        }

        return ss;
    };
    Float ss = compute_ss(a, b);

    for (int it = 0; it < iter; ++it) {
        // Computing the first and second derivatives of the sum of squared differences.
        Float da = 0, db = 0, daa = 0, dab = 0, dbb = 0;
        for (int g = 0; g < grid; ++g) {
            const Float& x = grid_x[g];
            const Float& gy = grid_y[g];
            const Float& oy = observed_y[g];

            const Float& x2b = xpow[g];
            const Float logx2 = log_x[g] * 2;
            const Float delta = oy - gy;

            // -(2 * (x^(2 * b)/(1 + a * x^(2 * b))^2 * (1/(1 + a * x^(2 * b)) - y)))
            da += -2 * x2b * oy * oy * delta;

            // -(2 * (a * (x^(2 * b) * (log(x) * 2))/(1 + a * x^(2 * b))^2 * (1/(1 + a * x^(2 * b)) - y)))
            db += -2 * a * x2b * logx2 * oy * oy * delta;

            // 2 * (
            //     x^(2 * b)/(1 + a * x^(2 * b))^2 * (x^(2 * b)/(1 + a * x^(2 * b))^2) 
            //     + x^(2 * b) * (2 * (x^(2 * b) * (1 + a * x^(2 * b))))/((1 + a * x^(2 * b))^2)^2 * (1/(1 + a * x^(2 * b)) - y)
            // ) 
            daa += 2 * (
                x2b * oy * oy * x2b * oy * oy
                + x2b * 2 * x2b * oy * oy * oy * delta
            );

            //-(2 * 
            //    (
            //        (
            //            (x^(2 * b) * (log(x) * 2))/(1 + a * x^(2 * b))^2 
            //            - a * (x^(2 * b) * (log(x) * 2)) * (2 * (x^(2 * b) * (1 + a * x^(2 * b))))/((1 + a * x^(2 * b))^2)^2
            //        ) 
            //        * (1/(1 + a * x^(2 * b)) - y) 
            //        - a * (x^(2 * b) * (log(x) * 2))/(1 + a * x^(2 * b))^2 * (x^(2 * b)/(1 + a * x^(2 * b))^2)
            //    )
            //)
            dab += -2 * (
                (
                    x2b * logx2 * oy * oy
                    - a * x2b * logx2 * 2 * x2b * oy * oy * oy
                ) * delta
                - a * x2b * logx2 * oy * oy * x2b * oy * oy
            );

            // -(2 * 
            //     (
            //         (
            //             a * (x^(2 * b) * (log(x) * 2) * (log(x) * 2))/(1 + a * x^(2 * b))^2 
            //             - a * (x^(2 * b) * (log(x) * 2)) * (2 * (a * (x^(2 * b) * (log(x) * 2)) * (1 + a * x^(2 * b))))/((1 + a * x^(2 * b))^2)^2
            //         ) 
            //         * (1/(1 + a * x^(2 * b)) - y) 
            //         - a * (x^(2 * b) * (log(x) * 2))/(1 + a * x^(2 * b))^2 * (a * (x^(2 * b) * (log(x) * 2))/(1 + a * x^(2 * b))^2)
            //     )
            // ) 
            dbb += -2 * (
                (
                    (a * x2b * logx2 * logx2 * oy * oy)
                    - (a * x2b * logx2 * 2 * a * x2b * logx2 * oy * oy * oy)
                ) * delta 
                - a * x2b * logx2 * oy * oy * a * x2b * logx2 * oy * oy
            );
        }

        // Applying the Newton iterations with damping.
        Float determinant = daa * dbb - dab * dab;
        const Float delta_a = (da * dbb - dab * db) / determinant;
        const Float delta_b = (- da * dab + daa * db) / determinant; 

        Float ss_next = 0;
        Float factor = 1;
        for (int inner = 0; inner < 10; ++inner, factor /= 2) {
            ss_next = compute_ss(a - factor * delta_a, b - factor * delta_b);
            if (ss_next < ss) {
                break;
            }
        }

        if (ss && 1 - ss_next/ss > tol) {
            a -= factor * delta_a;
            b -= factor * delta_b;
            ss = ss_next;
        } else {
            break;
        }
    }

    return std::make_pair(a, b);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List decomp_G(sp_mat &G, mat &initial_position, double a_param, double b_param, int n_epochs = 200, int thread_no = 0) {

    mat init_coors = initial_position.rows(0, 1);

    sp_mat H = G;
    H.for_each([](sp_mat::elem_type &val)
               { val = 1 / val; });
    H = ACTIONet::smoothKNN(H, 0);

    unsigned int nV = H.n_rows;

    // linearized list of edges (1-simplices)
    unsigned int nE = H.n_nonzero;
    vector<unsigned int> positive_head(nE);
    vector<unsigned int> positive_tail(nE);
    vector<float> epochs_per_sample(nE);

    int i = 0;
    double w_max = max(max(H));
    for (sp_mat::iterator it = H.begin(); it != H.end(); ++it)
    {
      epochs_per_sample[i] = w_max / (*it);
      positive_head[i] = it.row();
      positive_tail[i] = it.col();
      i++;
    }

    // Initial coordinates of vertices (0-simplices)
    vector<float> head_embedding(init_coors.n_cols * 2);
    fmat sub_coor = conv_to<fmat>::from(init_coors.rows(0, 1));
    float *ptr = sub_coor.memptr();
    memcpy(head_embedding.data(), ptr, sizeof(float) * head_embedding.size());
    vector<float> tail_embedding(head_embedding);

    sp_mat H2 = trans(H);
    H2.sync();
    std::vector<unsigned int> positive_ptr(H2.n_cols + 1);
    for (int k = 0; k < H2.n_cols + 1; k++)
    {
        positive_ptr[k] = H2.col_ptrs[k];
    }

/*
  vector<float> result;
  mat coordinates(nV, 2);
  result = optimize_layout_interface(
    head_embedding, tail_embedding,
    positive_head,
    positive_tail,
    positive_ptr, n_epochs,
    nV, nV,
    epochs_per_sample, "umap",
    1, a_param, b_param, 1, false, 5,
    true, true, thread_no);

  fmat coordinates_float(result.data(), 2, nV);
  coordinates = trans(conv_to<mat>::from(coordinates_float));


*/

  List res;
  //res["coordinates"] = coordinates;
  res["positive_head"] = positive_head;
  res["positive_tail"] = positive_tail;
  res["epochs_per_sample"] = epochs_per_sample;  
  res["positive_ptr"] = positive_ptr;
  res["head_embedding"] = head_embedding;
  res["tail_embedding"] = tail_embedding;

  return (res);
}    