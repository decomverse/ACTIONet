#include <numpy/npy_common.h>
#include <pybind11/stl.h>
#include <stdint.h>
#include <string.h>
#include <functional>
#include <string>
#include <utility>
#include "include/ACTIONet.h"
#include "include/arma_wrapper.h"

using aw::dcube;
using aw::dmat;
using aw::dvec;
using aw::intmat;
using aw::npdouble;
using aw::npint;

namespace py = pybind11;
using namespace py::literals;

// Computes SVD decomposition
//
// This is direct implementation of the randomized SVD algorithm:
// From: IRLBA R Package
//
// @param A Input matrix ("sparseMatrix")
// @param dim Dimension of SVD decomposition
// @param iters Number of iterations (default=1000)
// @param seed Random seed (default=0)
//
// @return A named dictionary with U, sigma, and V components
py::dict IRLB_SVD(arma::SpMat<npdouble> &A, int dim, int iters = 1000,
                  int seed = 0, int verbose = 1) {
  arma::field<arma::Mat<npdouble>> SVD_out =
      ACTIONet::IRLB_SVD(A, dim, iters, seed, verbose);
  py::dict res;
  res["u"] = SVD_out(0);
  res["d"] = SVD_out(1);
  res["v"] = SVD_out(2);
  return res;
}

py::dict IRLB_SVD_full(arma::Mat<npdouble> &A, int dim, int iters = 1000,
                       int seed = 0, int verbose = 1) {
  arma::field<arma::Mat<npdouble>> SVD_out =
      ACTIONet::IRLB_SVD(A, dim, iters, seed, verbose);
  py::dict res;
  res["u"] = SVD_out(0);
  res["d"] = SVD_out(1);
  res["v"] = SVD_out(2);
  return res;
}

// Computes SVD decomposition
//
// This is direct implementation of the randomized SVD algorithm for sparse
// matrices: Xu Feng, Yuyang Xie, and Yaohang Li, "Fast Randomzied SVD for
// Sparse Data," in Proc. the 10th Asian Conference on Machine Learning (ACML),
// Beijing, China, Nov. 2018.
//
// @param A Input matrix ("sparseMatrix")
// @param dim Dimension of SVD decomposition
// @param iters Number of iterations (default=5)
// @param seed Random seed (default=0)
//
// @return A named dictionary with U, sigma, and V components
py::dict FengSVD(arma::SpMat<npdouble> &A, int dim, int iters = 5, int seed = 0,
                 int verbose = 1) {
  field<arma::Mat<npdouble>> SVD_out =
      ACTIONet::FengSVD(A, dim, iters, seed, verbose);
  py::dict res;
  res["u"] = SVD_out(0);
  res["d"] = SVD_out(1);
  res["v"] = SVD_out(2);
  return res;
}
py::dict FengSVD_full(arma::Mat<npdouble> &A, int dim, int iters = 5,
                      int seed = 0, int verbose = 1) {
  field<arma::Mat<npdouble>> SVD_out =
      ACTIONet::FengSVD(A, dim, iters, seed, verbose);
  py::dict res;
  res["u"] = SVD_out(0);
  res["d"] = SVD_out(1);
  res["v"] = SVD_out(2);
  return res;
}

// Computes SVD decomposition
//
// This is direct implementation of the randomized SVD algorithm:
// From: N Halko, P. G Martinsson, and J. A Tropp. Finding structure with
// randomness: Probabilistic algorithms for constructing approximate matrix
// decompositions. Siam Review, 53(2):217-288, 2011.
//
// @param A Input matrix ("sparseMatrix")
// @param dim Dimension of SVD decomposition
// @param iters Number of iterations (default=5)
// @param seed Random seed (default=0)
//
// @return A named dictionary with U, sigma, and V components
py::dict HalkoSVD(arma::SpMat<npdouble> &A, int dim, int iters = 5,
                  int seed = 0, int verbose = 1) {
  field<arma::Mat<npdouble>> SVD_out =
      ACTIONet::HalkoSVD(A, dim, iters, seed, verbose);
  py::dict res;
  res["u"] = SVD_out(0);
  res["d"] = SVD_out(1);
  res["v"] = SVD_out(2);
  return res;
}
py::dict HalkoSVD_full(arma::Mat<npdouble> &A, int dim, int iters = 5,
                       int seed = 0, int verbose = 1) {
  field<arma::Mat<npdouble>> SVD_out =
      ACTIONet::HalkoSVD(A, dim, iters, seed, verbose);
  py::dict res;
  res["u"] = SVD_out(0);
  res["d"] = SVD_out(1);
  res["v"] = SVD_out(2);
  return res;
}

// Solves min_{X} (|| AX - B ||) s.t. simplex constraint
//
// @param A Input matrix (dense)
// @param B Input matrix (dense)
//
// @return X Solution
arma::Mat<npdouble> run_simplex_regression(arma::Mat<npdouble> &A,
                                           arma::Mat<npdouble> &B,
                                           bool computeXtX = false) {
  arma::Mat<npdouble> X = ACTIONet::run_simplex_regression(A, B, computeXtX);

  return X;
}

// Computes reduced kernel matrix for a given (single-cell) profile
//
// @param S Input matrix (dense)
// @param reduced_dim Dimension of the reduced kernel matrix (default=50)
// @param iters Number of SVD iterations (default=5)
// @param seed Random seed (default=0)
// @param reduction_algorithm Kernel reduction algorithm. Currently only ACTION
// method (1) is implemented (default=1)
// @param SVD_algorithm SVD algorithm to use. Currently supported methods are
// Halko (1) and Feng (2) (default=1)
//
// @return A named list with S_r, V, lambda, and exp_var. \itemize{
// \item S_r: reduced kernel matrix of size reduced_dim x #samples.
// \item V: Associated left singular-vectors (useful for reconstructing
// discriminative scores for features, such as genes). \item lambda, exp_var:
// Summary statistics of the sigular-values.
// }
py::dict reduce_kernel_full(arma::Mat<npdouble> &S, int reduced_dim = 50,
                            int iters = 5, int seed = 0, int SVD_algorithm = 0,
                            bool prenormalize = false, int verbose = 1) {
  field<arma::Mat<npdouble>> reduction = ACTIONet::reduce_kernel(
      S, reduced_dim, iters, seed, SVD_algorithm, prenormalize, verbose);

  py::dict res;

  res["V"] = reduction(0);

  vec sigma = reduction(1).col(0);
  res["sigma"] = sigma;

  double epsilon = 0.01 / sqrt(reduction(2).n_rows);
  arma::Mat<npdouble> V = round(reduction(2) / epsilon) * epsilon;

  for (int i = 0; i < V.n_cols; i++) {
    vec v = V.col(i) * sigma(i);
    V.col(i) = v;
  }
  V = trans(V);
  res["S_r"] = V.eval();

  res["A"] = reduction(3);
  res["B"] = reduction(4);

  return res;
}
py::dict reduce_kernel(arma::SpMat<npdouble> &S, int reduced_dim = 50,
                       int iters = 5, int seed = 0, int SVD_algorithm = 0,
                       bool prenormalize = false, int verbose = 1) {
  field<arma::Mat<npdouble>> reduction = ACTIONet::reduce_kernel(
      S, reduced_dim, iters, seed, SVD_algorithm, prenormalize, verbose);

  py::dict res;

  res["V"] = reduction(0);

  vec sigma = reduction(1).col(0);
  res["sigma"] = sigma;

  double epsilon = 0.01 / sqrt(reduction(2).n_rows);
  arma::Mat<npdouble> V = round(reduction(2) / epsilon) * epsilon;

  for (int i = 0; i < V.n_cols; i++) {
    vec v = V.col(i) * sigma(i);
    V.col(i) = v;
  }
  V = trans(V);
  res["S_r"] = V.eval();

  res["A"] = reduction(3);
  res["B"] = reduction(4);

  return res;
}

// Runs Successive Projection Algorithm (SPA) to solve separable NMF
//
// @param A Input matrix (dense)
// @param k Number of columns to select
//
// @return A named list with entries 'selected_columns' and 'norms'
py::dict run_SPA(arma::Mat<npdouble> &A, int k) {
  ACTIONet::SPA_results res = ACTIONet::run_SPA(A, k);
  uvec selected_columns = res.selected_columns;

  vec cols(k);
  for (int i = 0; i < k; i++) {
    cols[i] = selected_columns[i] + 1;
  }

  py::dict out;
  out["selected_columns"] = cols;
  out["norms"] = res.column_norms;

  return out;
}

// Runs Successive Projection Algorithm (SPA) to solve separable NMF
//
// @param A Input matrix (dense)
// @param k Number of columns to select
//
// @return A named list with entries 'selected_columns' and 'norms'
py::dict run_SPA_rows_sparse(arma::SpMat<npdouble> &A, int k) {
  ACTIONet::SPA_results res = ACTIONet::run_SPA_rows_sparse(A, k);
  uvec selected_columns = res.selected_columns;

  vec cols(k);
  for (int i = 0; i < k; i++) {
    cols[i] = selected_columns[i] + 1;
  }

  py::dict out;
  out["selected_columns"] = cols;
  out["norms"] = res.column_norms;

  return out;
}

// Runs Archetypal Analysis (AA)
//
// @param A Input matrix (dense)
// @param W0 Initial estimate of archetypes
// @param max_it, min_delta Define stopping conditions
//
// @return A named list with entries 'selected_columns' and 'norms'
py::dict run_AA(arma::Mat<npdouble> &A, arma::Mat<npdouble> &W0,
                int max_it = 50, double min_delta = 1e-16) {
  field<arma::Mat<npdouble>> AA_res =
      ACTIONet::run_AA(A, W0, max_it, min_delta);

  py::dict out;

  arma::Mat<npdouble> W = A * AA_res(0);
  out["C"] = AA_res(0);
  out["H"] = AA_res(1);
  out["W"] = W;

  return out;
}

// Runs multi-level ACTION decomposition method
//
// @param S_r Reduced kernel matrix
// @param k_min Minimum number of archetypes to consider (default=2)
// @param k_max Maximum number of archetypes to consider, or "depth" of
// decomposition (default=30)
// @param thread_no Number of parallel threads (default=0)
// @param max_it,min_delta Convergence parameters for archetypal analysis
//
// @return A named list with entries 'C' and 'H', each a list for different
// values of k
py::dict run_ACTION(arma::Mat<npdouble> &S_r, int k_min = 2, int k_max = 30,
                    int thread_no = 0, int max_it = 50,
                    double min_delta = 1e-16, int normalization = 0) {
  ACTIONet::ACTION_results trace = ACTIONet::run_ACTION(
      S_r, k_min, k_max, thread_no, max_it, min_delta, normalization);

  py::dict res;

  py::list C(k_max);
  for (int i = 0; i < k_max; i++) {
    arma::Mat<npdouble> curr_C = trace.C[i + 1];
    // curr_C = clamp(curr_C, 1e-5, 1);
    // curr_C = normalise(curr_C, 1);
    C[i] = curr_C;
  }
  res["C"] = C;

  py::list H(k_max);
  for (int i = 0; i < k_max; i++) {
    arma::Mat<npdouble> curr_H = trace.H[i + 1];
    // curr_H = clamp(curr_H, 1e-5, 1);
    // curr_H = normalise(curr_H, 1);
    H[i] = curr_H;
  }
  res["H"] = H;

  return res;
}

// Filters multi-level archetypes and concatenate filtered archetypes.
// (Pre-ACTIONet archetype processing)
//
// @param C_trace,H_trace Output of ACTION
// @param min_specificity_z_threshold Defines the stringency of pruning
// nonspecific archetypes. The larger the value, the more archetypes will be
// filtered out (default=-1)
//
// @return A named list: \itemize{
// \item selected_archs: py::dict of final archetypes that passed the
// filtering/pruning step. \item C_stacked,H_stacked: Horizontal/Vertical
// concatenation of filtered C and H matrices, respectively.
// }
py::dict prune_archetypes(vector<arma::Mat<npdouble>> &C_trace,
                          vector<arma::Mat<npdouble>> &H_trace,
                          double min_specificity_z_threshold = -3,
                          int min_cells = 3) {
  int n_list = H_trace.size();
  field<arma::Mat<npdouble>> C_trace_vec(n_list + 1);
  field<arma::Mat<npdouble>> H_trace_vec(n_list + 1);
  for (int i = 0; i < n_list; i++) {
    if (H_trace[i].is_empty()) {
      continue;
    }

    C_trace_vec[i] =
        C_trace[i];  // aw::conv_to<arma::Mat<npdouble>>::from(C_trace_list[i]);
    H_trace_vec[i] =
        H_trace[i];  // aw::conv_to<arma::Mat<npdouble>>::from(H_trace_list[i]);
  }

  ACTIONet::multilevel_archetypal_decomposition results =
      ACTIONet::prune_archetypes(C_trace_vec, H_trace_vec,
                                 min_specificity_z_threshold, min_cells);

  py::dict out_list;

  for (int i = 0; i < results.selected_archs.n_elem; i++)
    results.selected_archs[i]++;
  out_list["selected_archs"] = results.selected_archs;

  out_list["C_stacked"] = results.C_stacked;
  out_list["H_stacked"] = results.H_stacked;

  return out_list;
}

// Identifies and aggregates redundant archetypes into equivalent classes
// (Post-ACTIONet archetype processing)

py::dict unify_archetypes(arma::Mat<npdouble> &S_r,
                          arma::Mat<npdouble> &C_stacked,
                          arma::Mat<npdouble> &H_stacked,
                          double backbone_density, double resolution,
                          int min_cluster_size, int thread_no,
                          int normalization) {
  ACTIONet::unification_results results = ACTIONet::unify_archetypes(
      S_r, C_stacked, H_stacked, backbone_density, resolution, min_cluster_size,
      thread_no, normalization);
  py::dict out_list;

  for (int i = 0; i < results.selected_archetypes.n_elem; i++)
    results.selected_archetypes[i]++;
  out_list["selected_archetypes"] = results.selected_archetypes;

  out_list["C_unified"] = arma::SpMat<npdouble>(results.C_unified);
  out_list["H_unified"] = arma::SpMat<npdouble>(results.H_unified);

  for (int i = 0; i < results.assigned_archetypes.n_elem; i++)
    results.assigned_archetypes[i]++;

  out_list["assigned_archetype"] = results.assigned_archetypes;
  out_list["arch_membership_weights"] = results.arch_membership_weights;

  out_list["ontology"] = results.dag_adj;
  out_list["ontology_node_attributes"] = results.dag_node_annotations;

  return out_list;
}

// Computes node centrality scores
arma::vec compute_archetype_core_centrality(arma::SpMat<npdouble> &G,
                                            uvec assignments) {
  vec node_centrality =
      ACTIONet::compute_archetype_core_centrality(G, assignments);

  return node_centrality;
}

// Compute node coreness
arma::vec compute_core_number(arma::SpMat<npdouble> &G) {
  uvec core_num = ACTIONet::compute_core_number(G);

  return (conv_to<vec>::from(core_num));
}

// Builds an interaction network from the multi-level archetypal decompositions
//
// @param H_stacked Output of the prune_archetypes() function.
// @param algorithm Nearest neighbor alogirthm: k*nn, knn. (default="k*nn")
// @param distance_metric Distance metric to use: jsd, l2, ip. (default="jsd")
// @param density Overall density of constructed graph. The higher the density,
// the more edges are retained (default=1.0).
// @param thread_no Number of parallel threads (default=0).
// @param mutual_edges_only Symmetrization strategy for nearest-neighbor edges.
// If it is true, only mutual-nearest-neighbors are returned (default=true).
// @param k Optional parameter specifying k for knn algorithm (default=10).
// @param M 'M' parameter to pass to UMAP (default=16).
// @param ef_construction 'ef_construction' parameter to pass to UMAP
// (default=200).
// @param ef 'ef' parameter to pass to UMAP (default=50).
//
// @return G Adjacency matrix of the ACTIONet graph.
arma::SpMat<npdouble> buildNetwork(arma::Mat<npdouble> &H,
                                   string algorithm = "k*nn",
                                   string distance_metric = "jsd",
                                   double density = 1.0, int thread_no = 0,
                                   bool mutual_edges_only = true, int k = 10) {
  double M = 16, ef_construction = 200, ef = 50;

  arma::SpMat<npdouble> G =
      ACTIONet::buildNetwork(H, algorithm, distance_metric, density, thread_no,
                             M, ef_construction, ef, mutual_edges_only, k);
  return G;
}

// Performs stochastic force-directed layout on the input graph
//
// @param G Adjacency matrix of the ACTIONet graph
// @param initial_position Reduced kernel matrix (is used for reproducible
// initialization).
// @param algorithm Algorithm to use for visualization layout (default="TUMAP").
// @param compactness_level A value between 0-100, indicating the compactness of
// ACTIONet layout (default=50)
// @param n_epochs Number of epochs for SGD algorithm (default=500).
// @param thread_no Number of threads (default=0).
// @param seed Seed for random initialization (default=0).
//
// @return A named list \itemize{
// \item coordinates 2D coordinates of vertices.
// \item coordinates_3D 3D coordinates of vertices.
// \item colors De novo color of nodes inferred from their 3D embedding.
// }
py::dict layoutNetwork(arma::SpMat<npdouble> &G,
                       arma::Mat<npdouble> &initial_position,
                       const std::string &method = "umap",
                       bool presmooth_network = false, double min_dist = 1,
                       double spread = 1, double gamma = 1.0,
                       unsigned int n_epochs = 500, int thread_no = 0,
                       int seed = 0, double learning_rate = 1.0,
                       int sim2dist = 2) {
  field<arma::Mat<npdouble>> res = ACTIONet::layoutNetwork_xmap(
      G, initial_position, presmooth_network, method, min_dist, spread, gamma,
      n_epochs, thread_no, seed, learning_rate, sim2dist);

  py::dict out_list;
  out_list["coordinates"] = res(0);
  out_list["coordinates_3D"] = res(1);
  out_list["colors"] = res(2);

  return out_list;
}

// Computes graph clustering using Leiden algorith over signed graphs
//
// @param G Adjacency matrix of the input graph
// @param resolution_parameter Granularity of clustering. Larger values result
// in more clusters (default = 1.0)
// @param initial_clusters_ Initialization vector for clusters (if available)
// @param seed Random seed
//
// @return clusters Assignment vector of samples to clusters
arma::vec signed_cluster(arma::SpMat<npdouble> &A,
                         double resolution_parameter = 1.0,
                         uvec initial_clusters_ = uvec(), int seed = 0) {
  uvec initial_clusters_uvec(A.n_rows);
  if (initial_clusters_.n_elem == A.n_rows) {
    for (int i = 0; i < A.n_rows; i++)
      initial_clusters_uvec(i) = initial_clusters_(i);
  } else {
    for (int i = 0; i < A.n_rows; i++) initial_clusters_uvec(i) = i;
  }

  vec clusters = ACTIONet::signed_cluster(A, resolution_parameter,
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
arma::vec unsigned_cluster(arma::SpMat<npdouble> &A,
                           double resolution_parameter = 1.0,
                           arma::uvec initial_clusters_ = arma::uvec(),
                           int seed = 0) {
  uvec initial_clusters_uvec(A.n_rows);
  if (initial_clusters_.n_elem == A.n_rows) {
    for (int i = 0; i < A.n_rows; i++)
      initial_clusters_uvec(i) = initial_clusters_(i);
  } else {
    for (int i = 0; i < A.n_rows; i++) initial_clusters_uvec(i) = i;
  }

  vec clusters = ACTIONet::unsigned_cluster(A, resolution_parameter,
                                            initial_clusters_uvec, seed);

  return clusters;
}

arma::SpMat<npdouble> normalize_adj(arma::SpMat<npdouble> &G,
                                    int norm_type = 0) {
  arma::SpMat<npdouble> P = ACTIONet::normalize_adj(G, norm_type);

  return (P);
}

arma::Mat<npdouble> compute_network_diffusion_fast(arma::SpMat<npdouble> &G,
                                                   arma::SpMat<npdouble> &X0,
                                                   int thread_no = 0,
                                                   double alpha = 0.85,
                                                   int max_it = 5) {
  arma::Mat<npdouble> X = ACTIONet::compute_network_diffusion_fast(
      G = G, X0 = X0, thread_no = thread_no, alpha = alpha, max_it = max_it);

  return (X);
}

arma::Mat<npdouble> compute_network_diffusion_approx(
    arma::SpMat<npdouble> &G, arma::Mat<npdouble> X0, int thread_no = 0,
    double alpha = 0.85, int max_it = 5, double res_threshold = 1e-8,
    int norm_type = 0) {
  if (G.n_rows != X0.n_rows) {
    stderr_printf("Dimension mismatch: G (%dx%d) and X0 (%dx%d)\n", G.n_rows,
                  G.n_cols, X0.n_rows, X0.n_cols);
    return (arma::Mat<npdouble>());
  }

  arma::SpMat<npdouble> P = ACTIONet::normalize_adj(G, norm_type);

  arma::Mat<npdouble> X = ACTIONet::compute_network_diffusion_Chebyshev(
      P, X0, thread_no, alpha, max_it, res_threshold);

  return (X);
}

arma::vec run_LPA(arma::SpMat<npdouble> &G, arma::vec labels,
                  arma::vec fixed_labels, double lambda_param, int iters,
                  double sig_threshold, int thread_no = 0) {
  /*
arma::uvec fixed_labels_vec;

if (!fixed_labels.is_empty()) {
arma::uvec fixed_labels_vec(fixed_labels_.size());
for (int i = 0; i < fixed_labels_.size(); i++) {
fixed_labels_vec(i) = fixed_labels_(i);
}
}
*/
  arma::uvec fixed_labels_uvec = conv_to<uvec>::from(fixed_labels);

  arma::vec new_labels =
      ACTIONet::LPA(G, labels, lambda_param, iters, sig_threshold,
                    fixed_labels_uvec, thread_no);
  return (new_labels);
}

py::dict compute_archetype_feature_specificity_bin(arma::SpMat<npdouble> &S,
                                                   arma::Mat<npdouble> &H,
                                                   int thread_no = 0) {
  field<arma::Mat<npdouble>> res =
      ACTIONet::compute_feature_specificity_bin(S, H, thread_no);

  py::dict out_list;
  out_list["archetypes"] = res(0);
  out_list["upper_significance"] = res(1);
  out_list["lower_significance"] = res(2);

  return (out_list);
}

py::dict compute_archetype_feature_specificity(arma::SpMat<npdouble> &S,
                                               arma::Mat<npdouble> &H,
                                               int thread_no = 0) {
  field<arma::Mat<npdouble>> res =
      ACTIONet::compute_feature_specificity(S, H, thread_no);

  py::dict out_list;
  out_list["archetypes"] = res(0);
  out_list["upper_significance"] = res(1);
  out_list["lower_significance"] = res(2);

  return (out_list);
}

py::dict compute_archetype_feature_specificity_full(arma::Mat<npdouble> &S,
                                                    arma::Mat<npdouble> &H,
                                                    int thread_no = 0) {
  field<arma::Mat<npdouble>> res =
      ACTIONet::compute_feature_specificity(S, H, thread_no);

  py::dict out_list;
  out_list["archetypes"] = res(0);
  out_list["upper_significance"] = res(1);
  out_list["lower_significance"] = res(2);

  return (out_list);
}

py::dict compute_cluster_feature_specificity(arma::SpMat<npdouble> &S,
                                             arma::uvec sample_assignments,
                                             int thread_no = 0) {
  field<arma::Mat<npdouble>> res =
      ACTIONet::compute_feature_specificity(S, sample_assignments, thread_no);

  py::dict out_list;
  out_list["average_profile"] = res(0);
  out_list["upper_significance"] = res(1);
  out_list["lower_significance"] = res(2);

  return (out_list);
}

py::dict compute_cluster_feature_specificity_full(arma::Mat<npdouble> &S,
                                                  arma::uvec sample_assignments,
                                                  int thread_no = 0) {
  field<arma::Mat<npdouble>> res =
      ACTIONet::compute_feature_specificity(S, sample_assignments, thread_no);

  py::dict out_list;
  out_list["average_profile"] = res(0);
  out_list["upper_significance"] = res(1);
  out_list["lower_significance"] = res(2);

  return (out_list);
}

py::dict autocorrelation_Geary(arma::SpMat<npdouble> &G,
                               arma::Mat<npdouble> scores,
                               int normalization_method = 1, int perm_no = 30,
                               int thread_no = 0) {
  arma::field<arma::vec> out = ACTIONet::autocorrelation_Geary(
      G, scores, normalization_method, perm_no, thread_no);

  py::dict res;
  res["Geary_C"] = out[0];
  res["zscore"] = out[1];
  res["mu"] = out[2];
  res["sigma"] = out[3];

  return (res);
}

py::dict autocorrelation_Geary_full(arma::Mat<npdouble> &G,
                                    arma::Mat<npdouble> scores,
                                    int normalization_method = 1,
                                    int perm_no = 30, int thread_no = 0) {
  arma::field<arma::vec> out = ACTIONet::autocorrelation_Geary(
      G, scores, normalization_method, perm_no, thread_no);

  py::dict res;
  res["Geary_C"] = out[0];
  res["zscore"] = out[1];
  res["mu"] = out[2];
  res["sigma"] = out[3];

  return (res);
}

py::dict autocorrelation_Moran(arma::SpMat<npdouble> &G,
                               arma::Mat<npdouble> scores,
                               int normalization_method = 1, int perm_no = 30,
                               int thread_no = 0) {
  arma::field<arma::vec> out = ACTIONet::autocorrelation_Moran(
      G, scores, normalization_method, perm_no, thread_no);

  py::dict res;
  res["Moran_I"] = out[0];
  res["zscore"] = out[1];
  res["mu"] = out[2];
  res["sigma"] = out[3];

  return (res);
}

py::dict autocorrelation_Moran_full(arma::Mat<npdouble> &G,
                                    arma::Mat<npdouble> scores,
                                    int normalization_method = 1,
                                    int perm_no = 30, int thread_no = 0) {
  arma::field<arma::vec> out = ACTIONet::autocorrelation_Moran(
      G, scores, normalization_method, perm_no, thread_no);

  py::dict res;
  res["Moran_I"] = out[0];
  res["zscore"] = out[1];
  res["mu"] = out[2];
  res["sigma"] = out[3];

  return (res);
}

arma::Mat<npdouble> normalize_mat(arma::Mat<npdouble> &X, int normalization = 0,
                                  int dim = 0) {
  arma::Mat<npdouble> X_norm = ACTIONet::normalize_mat(X, normalization, dim);

  return (X_norm);
}

arma::SpMat<npdouble> normalize_spmat(arma::SpMat<npdouble> &X,
                                      int normalization = 0, int dim = 0) {
  arma::SpMat<npdouble> X_norm = ACTIONet::normalize_mat(X, normalization, dim);

  return (X_norm);
}

arma::Mat<npdouble> aggregate_genesets_mahalanobis_2archs(
    arma::SpMat<npdouble> &G, arma::SpMat<npdouble> &S,
    arma::SpMat<npdouble> &marker_mat, int network_normalization_method,
    int expression_normalization_method, int gene_scaling_method,
    double pre_alpha, double post_alpha, int thread_no) {
  arma::Mat<npdouble> stats = ACTIONet::aggregate_genesets_mahalanobis_2archs(
      G, S, marker_mat, network_normalization_method,
      expression_normalization_method, gene_scaling_method, pre_alpha,
      post_alpha, thread_no);

  return (stats);
}

arma::Mat<npdouble> aggregate_genesets_mahalanobis_2gmm(
    arma::SpMat<npdouble> &G, arma::SpMat<npdouble> &S,
    arma::SpMat<npdouble> &marker_mat, int network_normalization_method,
    int expression_normalization_method, int gene_scaling_method,
    double pre_alpha, double post_alpha, int thread_no) {
  arma::Mat<npdouble> stats = ACTIONet::aggregate_genesets_mahalanobis_2gmm(
      G, S, marker_mat, network_normalization_method,
      expression_normalization_method, gene_scaling_method, pre_alpha,
      post_alpha, thread_no);

  return (stats);
}

arma::Mat<npdouble> aggregate_genesets_weighted_enrichment(
    arma::SpMat<npdouble> &G, arma::SpMat<npdouble> &S,
    arma::SpMat<npdouble> &marker_mat, int network_normalization_method,
    int expression_normalization_method, int gene_scaling_method,
    double pre_alpha, double post_alpha, int thread_no) {
  arma::Mat<npdouble> stats = ACTIONet::aggregate_genesets_weighted_enrichment(
      G, S, marker_mat, network_normalization_method,
      expression_normalization_method, gene_scaling_method, pre_alpha,
      post_alpha, thread_no);

  return (stats);
}

arma::Mat<npdouble> aggregate_genesets_weighted_enrichment_permutation(
    arma::SpMat<npdouble> &G, arma::SpMat<npdouble> &S,
    arma::SpMat<npdouble> &marker_mat, int network_normalization_method,
    int expression_normalization_method, int gene_scaling_method,
    double pre_alpha, double post_alpha, int thread_no, int perm_no) {
  arma::Mat<npdouble> stats =
      ACTIONet::aggregate_genesets_weighted_enrichment_permutation(
          G, S, marker_mat, network_normalization_method,
          expression_normalization_method, gene_scaling_method, pre_alpha,
          post_alpha, thread_no, perm_no);

  return (stats);
}

py::dict assess_enrichment(arma::Mat<npdouble> &scores,
                           arma::SpMat<npdouble> &associations,
                           int thread_no = 0) {
  field<mat> res = ACTIONet::assess_enrichment(scores, associations, thread_no);

  py::dict out_list;
  out_list["logPvals"] = res[0];
  out_list["thresholds"] = res[1];

  return (out_list);
}

arma::vec xicor(arma::vec xvec, arma::vec yvec, bool compute_pval, int seed) {
  arma::vec res = ACTIONet::xicor(xvec, yvec, compute_pval, seed);

  return (res);
}

py::dict XICOR(arma::Mat<npdouble> &X, arma::Mat<npdouble> &Y,
               bool compute_pval, int seed, int thread_no) {
  field<arma::Mat<npdouble>> out =
      ACTIONet::XICOR(X, Y, compute_pval, seed, thread_no);

  py::dict res;
  res["XI"] = out[0];
  res["Z"] = out[1];

  return (res);
}

arma::Mat<npdouble> run_harmony(arma::Mat<npdouble> &X, arma::Mat<npdouble> &W0,
                                arma::vec batch, int clustering_algorithm,
                                double sigma_val, double theta_val,
                                int max_iter_cluster, int max_iter_harmony,
                                double eps_cluster, double eps_harmony,
                                double tau, double block_size,
                                double lambda_val, bool verbose, int seed) {
  mat Z_corr = ACTIONet::run_harmony(
      X, W0, batch, clustering_algorithm, sigma_val, theta_val,
      max_iter_cluster, max_iter_harmony, eps_cluster, eps_harmony, tau,
      block_size, lambda_val, verbose, seed);

  return (Z_corr);
}

py::dict aggregate_genesets(arma::sp_mat &G, arma::sp_mat &S,
                            arma::sp_mat &marker_mat,
                            int network_normalization_method, double alpha,
                            int thread_no) {
  arma::field<arma::mat> stats = ACTIONet::aggregate_genesets_vision(
      G, S, marker_mat, network_normalization_method, alpha, thread_no);

  py::dict res;

  res["stats_norm_smoothed"] = stats[0];
  res["stats_norm"] = stats[1];
  res["stats"] = stats[2];

  return (res);
}

arma::mat assess_label_enrichment(arma::sp_mat &G, arma::mat &M,
                                  int thread_no = 0) {
  arma::mat logPvals = ACTIONet::assess_label_enrichment(G, M, thread_no);

  return (logPvals);
}

PYBIND11_MODULE(_ACTIONet, m) {
  m.doc() = R"pbdoc(
        ACTIONet package
        -----------------------

        .. currentmodule:: ACTIONet

        .. autosummary::
           :toctree: _generate

    )pbdoc";
  // normalization
  m.def("normalize_mat", &normalize_mat, "Normalize input matrix.",
        py::arg("X"), py::arg("normalization") = 0, py::arg("dim") = 0);

  // normalization
  m.def("normalize_spmat", &normalize_spmat, "Normalize input matrix.",
        py::arg("X"), py::arg("normalization") = 0, py::arg("dim") = 0);

  // SVD
  m.def("IRLB_SVD", &IRLB_SVD, "Computes SVD using IRLB algorithm.",
        py::arg("A"), py::arg("dim"), py::arg("iters") = 1000,
        py::arg("seed") = 0, py::arg("verbose") = 1);

  m.def("IRLB_SVD_full", &IRLB_SVD_full, "Computes SVD using IRLB algorithm.",
        py::arg("A"), py::arg("dim"), py::arg("iters") = 1000,
        py::arg("seed") = 0, py::arg("verbose") = 1);

  m.def("FengSVD", &FengSVD, "Computes SVD using Feng et al. algorithm.",
        py::arg("A"), py::arg("dim"), py::arg("iters") = 5, py::arg("seed") = 0,
        py::arg("verbose") = 1);

  m.def("FengSVD_full", &FengSVD_full,
        "Computes SVD using Feng et al. algorithm.", py::arg("A"),
        py::arg("dim"), py::arg("iters") = 5, py::arg("seed") = 0,
        py::arg("verbose") = 1);

  m.def("HalkoSVD", &HalkoSVD, "Computes SVD using Halko et al. algorithm.",
        py::arg("A"), py::arg("dim"), py::arg("iters") = 5, py::arg("seed") = 0,
        py::arg("verbose") = 1);

  m.def("HalkoSVD_full", &HalkoSVD_full,
        "Computes SVD using Halko et al. algorithm.", py::arg("A"),
        py::arg("dim"), py::arg("iters") = 5, py::arg("seed") = 0,
        py::arg("verbose") = 1);

  // Kernel reduction
  m.def("reduce_kernel", &reduce_kernel,
        "Computes reduced kernel matrix for a given profile", py::arg("S"),
        py::arg("reduced_dim") = 50, py::arg("iters") = 1000,
        py::arg("seed") = 0, py::arg("SVD_algorithm") = 0,
        py::arg("prenormalize") = false, py::arg("verbose") = 1);

  m.def("reduce_kernel_full", &reduce_kernel_full,
        "Computes reduced kernel matrix for a given profile", py::arg("S"),
        py::arg("reduced_dim") = 50, py::arg("iters") = 1000,
        py::arg("seed") = 0, py::arg("SVD_algorithm") = 0,
        py::arg("prenormalize") = false, py::arg("verbose") = 1);

  // Lower-level functions
  m.def("run_simplex_regression", &run_simplex_regression,
        "Solves min_{X} (|| AX - B ||) s.t. simplex constraint", py::arg("A"),
        py::arg("B"), py::arg("computeXtX") = false);

  m.def("run_AA", &run_AA, "Runs Archetypal Analysis (AA) Algorithm",
        py::arg("A"), py::arg("W0"), py::arg("max_it") = 50,
        py::arg("min_delta") = 0.01);

  m.def("run_SPA", &run_SPA,
        "Runs Successive Projection Algorithm (SPA) to solve separable NMF",
        py::arg("A"), py::arg("k"));

  m.def("run_SPA_rows_sparse", &run_SPA_rows_sparse,
        "Runs Successive Projection Algorithm (SPA) to solve separable NMF",
        py::arg("A"), py::arg("k"));

  // ACTION decomposition
  m.def("run_ACTION", &run_ACTION,
        "Runs multi-level ACTION decomposition method", py::arg("S_r"),
        py::arg("k_min") = 2, py::arg("k_max") = 30, py::arg("thread_no") = 0,
        py::arg("max_it") = 50, py::arg("min_delta") = 0.01,
        py::arg("normalization") = 1);

  // Archetypes
  m.def("prune_archetypes", &prune_archetypes,
        "Filters multi-level archetypes and concatenate filtered archetypes",
        py::arg("C_trace"), py::arg("H_trace"),
        py::arg("min_specificity_z_threshold") = -3, py::arg("min_cells") = 3);

  m.def(
      "unify_archetypes", &unify_archetypes,
      "Identifies and aggregates redundant archetypes into equivalent classes",
      py::arg("S_r"), py::arg("C_stacked"), py::arg("H_stacked"),
      py::arg("backbone_density") = 0.5, py::arg("resolution") = 1.0,
      py::arg("min_cluster_size") = 3, py::arg("thread_no") = 0,
      py::arg("normalization") = 0);

  m.def("compute_archetype_core_centrality", &compute_archetype_core_centrality,
        "Computes node centrality scores based on localized coreness",
        py::arg("G"), py::arg("sample_assignments"));

  m.def("compute_core_number", &compute_core_number,
        "Computes node centrality scores based on coreness", py::arg("G"));

  // Network
  m.def("buildNetwork", &buildNetwork,
        "Builds an interaction network from a decomposition factor (H)"
        "decompositions",
        py::arg("H"), py::arg("algorithm") = "k*nn",
        py::arg("distance_metric") = "jsd", py::arg("density") = 1.0,
        py::arg("thread_no") = 0, py::arg("mutual_edges_only") = true,
        py::arg("k") = 10);

  m.def("layoutNetwork", &layoutNetwork,
        "Performs stochastic force-directed layout on the input graph",
        py::arg("G"), py::arg("initial_position"), py::arg("method") = "umap",
        py::arg("presmooth_network") = false, py::arg("min_dist") = 1.0,
        py::arg("spread") = 1.0, py::arg("min_dist") = 1.0,
        py ::arg("n_epochs") = 1000, py::arg("thread_no") = 0,
        py::arg("seed") = 0, py::arg("learning_rate") = 1.0,
        py::arg("sim2dist") = 2);

  m.def("signed_cluster", &unsigned_cluster,
        "Computes graph clustering using Leiden "
        "algorith over signed graphs",
        py::arg("A"), py::arg("resolution_parameter") = 1.0,
        py::arg("initial_clusters") = uvec(), py::arg("seed") = 0);

  m.def("unsigned_cluster", &unsigned_cluster,
        "Computes graph clustering using Leiden algorith over unsigned graphs",
        py::arg("A"), py::arg("resolution_parameter") = 1.0,
        py::arg("initial_clusters") = uvec(), py::arg("seed") = 0);

  m.def("normalize_adj", &normalize_adj,
        "Normalizes adjacency matrix using different strategies", py::arg("G"),
        py::arg("norm_type") = 0);

  m.def("compute_network_diffusion_fast", &compute_network_diffusion_fast,
        "Computes network diffusion using a given adjacency matrix",
        py::arg("G"), py::arg("X0"), py::arg("thread_no") = 0,
        py::arg("alpha") = 0.85, py::arg("max_it") = 5);

  m.def("compute_network_diffusion_approx", &compute_network_diffusion_approx,
        "Computes network diffusion using a given adjacency matrix",
        py::arg("G"), py::arg("X0"), py::arg("thread_no") = 0,
        py::arg("alpha") = 0.85, py::arg("max_it") = 5,
        py::arg("res_threshold") = 1e-8, py::arg("norm_type") = 1);

  m.def("run_LPA", &run_LPA,
        "Run label prepagation on a given set of known labels", py::arg("G"),
        py::arg("labels"), py::arg("fixed_labels") = arma::vec(),
        py::arg("lambda_param") = 1, py::arg("iters") = 3,
        py::arg("sig_threshold") = 3, py::arg("thread_no") = 0);

  m.def("compute_archetype_feature_specificity_bin",
        &compute_archetype_feature_specificity_bin,
        "Computes feature specificity of genes for each archetype",
        py::arg("S"), py::arg("H"), py::arg("thread_no") = 0);

  m.def("compute_archetype_feature_specificity",
        &compute_archetype_feature_specificity,
        "Computes feature specificity of genes for each archetype",
        py::arg("S"), py::arg("H"), py::arg("thread_no") = 0);

  m.def("compute_archetype_feature_specificity_full",
        &compute_archetype_feature_specificity_full,
        "Computes feature specificity of genes for each archetype",
        py::arg("S"), py::arg("H"), py::arg("thread_no") = 0);

  m.def("compute_cluster_feature_specificity",
        &compute_cluster_feature_specificity,
        "Computes feature specificity of genes for each archetype",
        py::arg("S"), py::arg("sample_assignments"), py::arg("thread_no") = 0);

  m.def("compute_cluster_feature_specificity_full",
        &compute_cluster_feature_specificity_full,
        "Computes feature specificity of genes for each archetype",
        py::arg("S"), py::arg("sample_assignments"), py::arg("thread_no") = 0);

  m.def("autocorrelation_Geary", &autocorrelation_Geary,
        "Computes spatial (network) autocorrelation (Geary)", py::arg("G"),
        py::arg("scores"), py::arg("normalization_method") = 1,
        py::arg("perm_no") = 30, py::arg("thread_no") = 0);

  m.def("autocorrelation_Moran", &autocorrelation_Moran,
        "Computes spatial (network) autocorrelation (Moran)", py::arg("G"),
        py::arg("scores"), py::arg("normalization_method") = 1,
        py::arg("perm_no") = 30, py::arg("thread_no") = 0);

  m.def("autocorrelation_Geary_full", &autocorrelation_Geary_full,
        "Computes spatial (network) autocorrelation (Geary)", py::arg("G"),
        py::arg("scores"), py::arg("normalization_method") = 1,
        py::arg("perm_no") = 30, py::arg("thread_no") = 0);

  m.def("autocorrelation_Moran_full", &autocorrelation_Moran_full,
        "Computes spatial (network) autocorrelation (Moran)", py::arg("G"),
        py::arg("scores"), py::arg("normalization_method") = 1,
        py::arg("perm_no") = 30, py::arg("thread_no") = 0);

  m.def("aggregate_genesets_mahalanobis_2archs",
        &aggregate_genesets_mahalanobis_2archs,
        "Scoring genesets using archetype mixture model approach", py::arg("G"),
        py::arg("S"), py::arg("marker_mat"),
        py::arg("network_normalization_method") = 2,
        py::arg("expression_normalization_method") = 1,
        py::arg("gene_scaling_method") = 0, py::arg("pre_alpha") = 0.15,
        py::arg("post_alpha") = 0.9, py::arg("thread_no") = 0);

  m.def("aggregate_genesets_mahalanobis_2gmm",
        &aggregate_genesets_mahalanobis_2gmm,
        "Scoring genesets using gaussian mixture model approach", py::arg("G"),
        py::arg("S"), py::arg("marker_mat"),
        py::arg("network_normalization_method") = 2,
        py::arg("expression_normalization_method") = 1,
        py::arg("gene_scaling_method") = 0, py::arg("pre_alpha") = 0.15,
        py::arg("post_alpha") = 0.9, py::arg("thread_no") = 0);

  m.def("aggregate_genesets_weighted_enrichment",
        &aggregate_genesets_weighted_enrichment,
        "Scoring genesets a parametric enrichment method", py::arg("G"),
        py::arg("S"), py::arg("marker_mat"),
        py::arg("network_normalization_method") = 2,
        py::arg("expression_normalization_method") = 1,
        py::arg("gene_scaling_method") = 0, py::arg("pre_alpha") = 0.15,
        py::arg("post_alpha") = 0.9, py::arg("thread_no") = 0);

  m.def("aggregate_genesets_weighted_enrichment_permutation",
        &aggregate_genesets_weighted_enrichment_permutation,
        "Scoring genesets a nonparametric enrichment method (permutation test)",
        py::arg("G"), py::arg("S"), py::arg("marker_mat"),
        py::arg("network_normalization_method") = 2,
        py::arg("expression_normalization_method") = 1,
        py::arg("gene_scaling_method") = 0, py::arg("pre_alpha") = 0.15,
        py::arg("post_alpha") = 0.9, py::arg("thread_no") = 0,
        py::arg("perm_no") = 30);

  m.def("assess_enrichment", &assess_enrichment,
        "Assess enrichment using a fast GSEA approach with sparsity "
        "considerations",
        py::arg("scores"), py::arg("associations"), py::arg("thread_no") = 0);

  m.def("xicor", &xicor, "XI correlation measure", py::arg("x"), py::arg("y"),
        py::arg("compute_pval") = true, py::arg("seed") = 0);

  m.def("XICOR", &XICOR, "XI correlation measure", py::arg("X"), py::arg("Y"),
        py::arg("compute_pval") = true, py::arg("seed") = 0,
        py::arg("thread_no") = 0);

  m.def("run_harmony", &run_harmony, "Run Harmony", py::arg("X"), py::arg("W0"),
        py::arg("batch"), py::arg("clustering_algorithm") = 1,
        py::arg("sigma_val") = 0.1, py::arg("theta_val") = 2.0,
        py::arg("max_iter_cluster") = 200, py::arg("max_iter_harmony") = 10,
        py::arg("eps_cluster") = 1e-5, py::arg("eps_harmony") = 1e-4,
        py::arg("tau") = 0, py::arg("block_size") = 0.05,
        py::arg("lambda_val") = 1.0, py::arg("verbose") = true,
        py::arg("seed") = 0);

  m.def("aggregate_genesets", &aggregate_genesets,
        "Aggregates genesets per cell", py::arg("G"), py::arg("S"),
        py::arg("marker_mat"), py::arg("network_normalization_method") = 0,
        py::arg("alpha") = 0.85, py::arg("thread_no") = 0);

  m.def("assess_label_enrichment", &assess_label_enrichment,
        "Assess enrichment of labels around each node in a network",
        py::arg("G"), py::arg("M"), py::arg("thread_no") = 0);

#ifdef VERSION_INFO
  m.attr("__version__") = VERSION_INFO;
#else
  m.attr("__version__") = "dev";
#endif
}
