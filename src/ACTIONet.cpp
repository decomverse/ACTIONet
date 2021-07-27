#include <ACTIONet.h>
#include <RcppArmadillo.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

#define ARMA_USE_CXX11_RNG
#define DYNSCHED

// set seed
// [[Rcpp::export]]
void set_seed(double seed)
{
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(std::floor(std::fabs(seed)));
}

template <class T1, class T2>
bool kv_pair_less(const std::pair<T1, T2> &x, const std::pair<T1, T2> &y)
{
  return x.first < y.first;
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
                      double min_specificity_z_threshold = -1,
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
//'	G = build_ACTIONet(prune.out$H_stacked)
//' unification.out = unify_archetypes(G, S_r, prune.out$C_stacked,
// prune.out$H_stacked) ' cell.clusters = unification.out$sample_assignments
// [[Rcpp::export]]
List unify_archetypes(mat &S_r, mat &C_stacked, mat &H_stacked,
                      double violation_threshold = 0.0,
                      int thread_no = 0)
{
  ACTIONet::unification_results results = ACTIONet::unify_archetypes(S_r, C_stacked, H_stacked, violation_threshold, thread_no);

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
//'	G = build_ACTIONet(prune.out$H_stacked)
// [[Rcpp::export]]
sp_mat build_ACTIONet(mat H_stacked, double density = 1.0, int thread_no = 0,
                      bool mutual_edges_only = true, string distance_metric = "jsd", string nn_approach = "k*nn", int k = 10)
{

  double M = 16, ef_construction = 200, ef = 50;
  sp_mat G = ACTIONet::build_ACTIONet(H_stacked, density, thread_no, M,
                                      ef_construction, ef, mutual_edges_only, distance_metric, nn_approach, k);

  return G;
}

// [[Rcpp::export]]
sp_mat build_knn(mat H_stacked, double k = 10, int thread_no = 0,
                 bool mutual_edges_only = true, string distance_metric = "jsd")
{

  double M = 16, ef_construction = 200, ef = 10;

  sp_mat G = ACTIONet::build_ACTIONet(H_stacked, 1, thread_no, M,
                                      ef_construction, ef, mutual_edges_only, distance_metric, "knn", k);

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
//'	G = build_ACTIONet(prune.out$H_stacked)
//'	vis.out = layout_ACTIONet(G, S_r)
// [[Rcpp::export]]
List layout_ACTIONet(sp_mat &G, mat S_r, int compactness_level = 50,
                     unsigned int n_epochs = 500, int layout_alg = 0, int thread_no = 0, int seed = 0)
{
  field<mat> res =
      ACTIONet::layout_ACTIONet(G, S_r, compactness_level, n_epochs, layout_alg, thread_no, seed);

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
//' @examples
//'	encoded.ids = encode_ids(colnames(sce))
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
//' @examples
//'	ids = decode_ids(encoded.ids)
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

//' Computes pseudobulk profiles
//'
//' @param S Input matrix ("sparseMatrix")
//' @param sample_assignments Any sample clustering/annotation (it has to be in
//{1, ..., max_class_num})
//'
//' @return S matrix aggregated within each class of sample_assignments
//'
//' @examples
//' prune.out = prune_archetypes(ACTION.out$C, ACTION.out$H)
//'	G = build_ACTIONet(prune.out$H_stacked)
//' unification.out = unify_archetypes(G, S_r, prune.out$C_stacked,
// prune.out$H_stacked) ' cell.clusters = unification.out$sample_assignments '
// pbs = compute_pseudo_bulk(S, cell.clusters)
// [[Rcpp::export]]
mat compute_pseudo_bulk_per_cluster(sp_mat &S,
                                    arma::Col<unsigned long long> sample_assignments)
{
  mat pb = ACTIONet::compute_pseudo_bulk_per_cluster(S, sample_assignments);

  return pb;
}

//' Computes pseudobulk profiles
//'
//' @param S Input matrix ("matrix")
//' @param sample_assignments Any sample clustering/annotation (it has to be in
//{1, ..., max_class_num})
//'
//' @return S matrix aggregated within each class of sample_assignments
//'
//' @examples
//' prune.out = prune_archetypes(ACTION.out$C, ACTION.out$H)
//'	G = build_ACTIONet(prune.out$H_stacked)
//' unification.out = unify_archetypes(G, S_r, prune.out$C_stacked,
// prune.out$H_stacked) ' cell.clusters = unification.out$sample_assignments '
// pbs = compute_pseudo_bulk(S, cell.clusters)
// [[Rcpp::export]]
mat compute_pseudo_bulk_per_cluster_full(mat &S,
                                         arma::Col<unsigned long long> sample_assignments)
{
  mat pb = ACTIONet::compute_pseudo_bulk_per_cluster(S, sample_assignments);

  return pb;
}

//' Computes pseudobulk profiles (groups[k1] x individuals[k2])
//'
//' @param S Input matrix ("sparseMatrix")
//' @param sample_assignments Any primary grouping - typically based on cell
// type/state (it has to be in {1, ..., k1}) ' @param individuals Any Secondary
// grouping - typically corresponds to individuals (it has to be in {1, ...,
// k2})
//'
//' @return A list of pseudobulk profile, where each entry is matrix
// corresponding to one cell type/state
//'
//' @examples
//' prune.out = prune_archetypes(ACTION.out$C, ACTION.out$H)
//'	G = build_ACTIONet(prune.out$H_stacked)
//' unification.out = unify_archetypes(G, S_r, prune.out$C_stacked,
// prune.out$H_stacked) ' cell.clusters = unification.out$sample_assignments '
// pbs.list = compute_pseudo_bulk(S, cell.clusters, sce$individuals)
// [[Rcpp::export]]
field<mat> compute_pseudo_bulk_per_cluster_and_ind(
    sp_mat &S, arma::Col<unsigned long long> sample_assignments,
    arma::Col<unsigned long long> individuals)
{
  field<mat> pbs_list =
      ACTIONet::compute_pseudo_bulk_per_cluster_and_ind(S, sample_assignments, individuals);

  return pbs_list;
}

//' Computes pseudobulk profiles (groups[k1] x individuals[k2])
//'
//' @param S Input matrix ("matrix")
//' @param sample_assignments Any primary grouping - typically based on cell
// type/state (it has to be in {1, ..., k1}) ' @param individuals Any Secondary
// grouping - typically corresponds to individuals (it has to be in {1, ...,
// k2})
//'
//' @return A list of pseudobulk profile, where each entry is matrix
// corresponding to one cell type/state
//'
//' @examples
//' prune.out = prune_archetypes(ACTION.out$C, ACTION.out$H)
//'	G = build_ACTIONet(prune.out$H_stacked)
//' unification.out = unify_archetypes(G, S_r, prune.out$C_stacked,
// prune.out$H_stacked) ' cell.clusters = unification.out$sample_assignments '
// pbs.list = compute_pseudo_bulk(S, cell.clusters, sce$individuals)
// [[Rcpp::export]]
field<mat> compute_pseudo_bulk_per_cluster_and_ind_full(
    mat &S, arma::Col<unsigned long long> sample_assignments,
    arma::Col<unsigned long long> individuals)
{
  field<mat> pbs_list =
      ACTIONet::compute_pseudo_bulk_per_cluster_and_ind(S, sample_assignments, individuals);

  return pbs_list;
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
//'	G = build_ACTIONet(prune.out$H_stacked)
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
//'	G = build_ACTIONet(prune.out$H_stacked)
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
List compute_archetype_feature_specificity_bin(sp_mat &S, mat &H)
{
  field<mat> res = ACTIONet::compute_feature_specificity_bin(S, H);

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
//'	G = build_ACTIONet(prune.out$H_stacked)
//' unification.out = unify_archetypes(G, S_r, prune.out$C_stacked,
// prune.out$H_stacked) ' cell.clusters = unification.out$sample_assignments '
// S.norm = renormalize_input_matrix(S, cell.clusters) '	logPvals.list =
// compute_archetype_feature_specificity(S.norm, unification.out$H_unified) '
// specificity.scores = logPvals.list$upper_significance
// [[Rcpp::export]]
List compute_archetype_feature_specificity(sp_mat &S, mat &H)
{
  field<mat> res = ACTIONet::compute_feature_specificity(S, H);

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
//'	G = build_ACTIONet(prune.out$H_stacked)
//' unification.out = unify_archetypes(G, S_r, prune.out$C_stacked,
// prune.out$H_stacked) ' cell.clusters = unification.out$sample_assignments '
// S.norm = renormalize_input_matrix(S, cell.clusters) '	logPvals.list =
// compute_archetype_feature_specificity(S.norm, unification.out$H_unified) '
// specificity.scores = logPvals.list$upper_significance
// [[Rcpp::export]]
List compute_archetype_feature_specificity_full(mat &S, mat &H)
{
  field<mat> res = ACTIONet::compute_feature_specificity(S, H);

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
//'	G = build_ACTIONet(prune.out$H_stacked)
//' unification.out = unify_archetypes(G, S_r, prune.out$C_stacked,
// prune.out$H_stacked) ' cell.clusters = unification.out$sample_assignments '
// S.norm = renormalize_input_matrix(S, cell.clusters) '	logPvals.list =
// compute_cluster_feature_specificity(S.norm, cell.clusters) '
// specificity.scores = logPvals.list$upper_significance
// [[Rcpp::export]]
List compute_cluster_feature_specificity(sp_mat &S, uvec sample_assignments)
{
  field<mat> res = ACTIONet::compute_feature_specificity(S, sample_assignments);

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
//'	G = build_ACTIONet(prune.out$H_stacked)
//' unification.out = unify_archetypes(G, S_r, prune.out$C_stacked,
// prune.out$H_stacked) ' cell.clusters = unification.out$sample_assignments '
// S.norm = renormalize_input_matrix(S, cell.clusters) '	logPvals.list =
// compute_cluster_feature_specificity(S.norm, cell.clusters) '
// specificity.scores = logPvals.list$upper_significance
// [[Rcpp::export]]
List compute_cluster_feature_specificity_full(mat &S, uvec sample_assignments)
{
  field<mat> res = ACTIONet::compute_feature_specificity(S, sample_assignments);

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
mat compute_network_diffusion(sp_mat &G, sp_mat &X0, int thread_no = 0,
                              double alpha = 0.85, int max_it = 3)
{
  mat Diff =
      ACTIONet::compute_network_diffusion(G, X0, thread_no, alpha, max_it);

  return (Diff);
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

//' Computes sparse network diffusion over a given network, starting with an
// arbitrarty set of initial scores
//'
//' @param G Input graph
//' @param X0 Matrix of initial values per diffusion (ncol(G) == nrow(G) ==
// ncol(X0)) ' @param alpha Random-walk depth ( between [0, 1] ) ' @param rho
// Sparsity controling parameter ' @param epsilon,max_it Conditions on the
// length of diffusion
//'
//' @return Matrix of sparse diffusion scores
//'
//' @examples
//' G = colNets(ace)$ACTIONet
//' gene.expression = Matrix::t(logcounts(ace))[c("CD19", "CD14", "CD16"), ]
//' smoothed.expression = compute_sparse_network_diffusion(G, gene.expression)
// [[Rcpp::export]]
sp_mat compute_sparse_network_diffusion(sp_mat &G, sp_mat &X0,
                                        double alpha = 0.85, double rho = 1e-4,
                                        double epsilon = 0.001,
                                        int max_iter = 20)
{
  sp_mat U = X0.transform([](double val)
                          { return (val < 0 ? 0 : val); });
  U = normalise(U, 1, 0);

  sp_mat scores = ACTIONet::compute_sparse_network_diffusion(G, U, alpha, rho,
                                                             epsilon, max_iter);
  /*
          sp_mat scores(size(X0));
          for(int i = 0; i < X0.n_cols; i++) {
                  vec v = vec(X0.col(i));
                  vec pr = ACTIONet::solve_fista(alpha, rho, epsilon, G, v, 20);
                  scores.col(i) = pr;
          }
  */

  return (scores);
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
mat Prune_PageRank(mat &U, double density = 1.0)
{
  mat G_matched = ACTIONet::Prune_PageRank(U, density);

  return G_matched;
}

// [[Rcpp::export]]
List transform_layout(sp_mat &W, mat coor2D, mat coor3D, mat colRGB,
                      int compactness_level = 50, unsigned int n_epochs = 500,
                      int thread_no = 0, int seed = 0)
{
  field<mat> res = ACTIONet::transform_layout(
      W, coor2D, coor3D, colRGB, compactness_level, n_epochs, thread_no, seed);

  List out_list;
  out_list["coordinates"] = res(0);
  out_list["coordinates_3D"] = res(1);
  out_list["colors"] = res(2);

  return out_list;
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
void csr_sort_indices_inplace(IntegerVector &Ap, IntegerVector &Aj,
                              NumericVector &Ax)
{
  int n_row = Ap.size() - 1;
  std::vector<std::pair<int, double>> temp;

  for (int i = 0; i < n_row; i++)
  {
    int row_start = (int)Ap[i];
    int row_end = (int)Ap[i + 1];
    int len = row_end - row_start;

    temp.resize(len);
    bool is_sorted = true;
    for (int jj = row_start, n = 0; jj < row_end; jj++, n++)
    {
      temp[n].first = (int)Aj(jj);
      temp[n].second = Ax(jj);
      if ((jj < (row_end - 1)) && (Aj(jj + 1) < Aj(jj)))
      {
        is_sorted = false;
      }
    }
    if (is_sorted)
      continue;

    std::sort(temp.begin(), temp.begin() + len, kv_pair_less<int, double>);
    for (int jj = row_start, n = 0; jj < row_end; jj++, n++)
    {
      Aj(jj) = temp[n].first;
      Ax(jj) = temp[n].second;
    }
  }
}

// [[Rcpp::export]]
void csc_sort_indices_inplace(IntegerVector &Ap, IntegerVector &Ai,
                              NumericVector &Ax)
{
  int n_col = Ap.size() - 1;

  std::vector<std::pair<int, double>> temp;
  for (int i = 0; i < n_col; i++)
  {
    int col_start = (int)Ap[i];
    int col_end = (int)Ap[i + 1];
    int len = col_end - col_start;

    temp.resize(len);
    bool is_sorted = true;
    for (int jj = col_start, n = 0; jj < col_end; jj++, n++)
    {
      temp[n].first = (int)Ai(jj);
      temp[n].second = Ax(jj);
      if ((jj < (col_end - 1)) && (Ai(jj + 1) < Ai(jj)))
      {
        is_sorted = false;
      }
    }
    if (is_sorted)
      continue;

    std::sort(temp.begin(), temp.begin() + len, kv_pair_less<int, double>);
    for (int jj = col_start, n = 0; jj < col_end; jj++, n++)
    {
      Ai(jj) = temp[n].first;
      Ax(jj) = temp[n].second;
    }
  }
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
      ACTIONet::orthogonalize_batch_effect(S, SVD_results, basal);

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
      ACTIONet::orthogonalize_batch_effect(S, SVD_results, basal);

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
vec run_LPA(sp_mat &G, vec labels, double lambda = 1, int iters = 3, double sig_threshold = 3, Nullable<IntegerVector> fixed_labels_ = R_NilValue)
{
  uvec fixed_labels_vec;
  if (fixed_labels_.isNotNull())
  {
    NumericVector fixed_labels(fixed_labels_);
    uvec fixed_labels_vec(fixed_labels.size());
    for (int i = 0; i < fixed_labels.size(); i++)
    {
      fixed_labels_vec(i) = fixed_labels(i) - 1;
    }
  }
  return (ACTIONet::LPA(G, labels, lambda, iters, sig_threshold, fixed_labels_vec));
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
mat compute_marker_aggregate_stats_basic_sum(sp_mat &S, sp_mat &marker_mat)
{
  mat stats = ACTIONet::compute_marker_aggregate_stats_basic_sum(S, marker_mat);

  return (stats);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat compute_marker_aggregate_stats_basic_sum_perm(sp_mat &S, sp_mat &marker_mat, int perm_no = 100, int thread_no = 0)
{
  mat stats = ACTIONet::compute_marker_aggregate_stats_basic_sum_perm(S, marker_mat, perm_no, thread_no);

  return (stats);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat compute_marker_aggregate_stats_basic_sum_perm_smoothed(sp_mat &G, sp_mat &S, sp_mat &marker_mat, double alpha = 0.85, int max_it = 5, int perm_no = 100, int thread_no = 0)
{
  mat stats = ACTIONet::compute_marker_aggregate_stats_basic_sum_perm_smoothed(G, S, marker_mat, alpha, max_it, perm_no, thread_no);

  return (stats);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat compute_marker_aggregate_stats_basic_sum_perm_smoothed_v2(sp_mat &G, sp_mat &S, sp_mat &marker_mat, double alpha = 0.85, int max_it = 5, int perm_no = 100, int thread_no = 0)
{
  mat stats = ACTIONet::compute_marker_aggregate_stats_basic_sum_perm_smoothed(G, S, marker_mat, alpha, max_it, perm_no, thread_no);

  return (stats);
}
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat compute_marker_aggregate_stats_basic_sum_smoothed(sp_mat &G, sp_mat &S, sp_mat &marker_mat, double alpha = 0.85, int max_it = 5, int perm_no = 100, int thread_no = 0)
{
  mat stats = ACTIONet::compute_marker_aggregate_stats_basic_sum_smoothed(G, S, marker_mat, alpha, max_it, perm_no, thread_no);

  return (stats);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat compute_marker_aggregate_stats_basic_sum_smoothed_normalized(sp_mat &G, sp_mat &S, sp_mat &marker_mat, double alpha = 0.85, int max_it = 5, int perm_no = 100, int thread_no = 0)
{
  mat stats = ACTIONet::compute_marker_aggregate_stats_basic_sum_smoothed_normalized(G, S, marker_mat, alpha, max_it, perm_no, thread_no);

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
