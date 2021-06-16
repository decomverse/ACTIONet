#include <ACTIONet.h>
#include <arma_wrapper.h>
#include <numpy/npy_common.h>
#include <pybind11/stl.h>
#include <stdint.h>
#include <string.h>
#include <string>
#include <functional>
#include <utility>

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
  field<mat> SVD_out = ACTIONet::IRLB_SVD(A, dim, iters, seed, verbose);
  py::dict res;
  res["u"] = SVD_out(0);
  res["d"] = SVD_out(1);
  res["v"] = SVD_out(2);
  return res;
}
py::dict IRLB_SVD_full(arma::Mat<npdouble> &A, int dim, int iters = 1000,
                       int seed = 0, int verbose = 1) {
  field<mat> SVD_out = ACTIONet::IRLB_SVD(A, dim, iters, seed, verbose);
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
py::dict FengSVD(arma::SpMat<npdouble> &A, int dim, int iters = 5,
                 int seed = 0, int verbose = 1) {
  field<mat> SVD_out = ACTIONet::FengSVD(A, dim, iters, seed, verbose);
  py::dict res;
  res["u"] = SVD_out(0);
  res["d"] = SVD_out(1);
  res["v"] = SVD_out(2);
  return res;
}
py::dict FengSVD_full(arma::Mat<npdouble> &A, int dim, int iters = 5,
                      int seed = 0, int verbose = 1) {
  field<mat> SVD_out = ACTIONet::FengSVD(A, dim, iters, seed, verbose);
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
  field<mat> SVD_out = ACTIONet::HalkoSVD(A, dim, iters, seed, verbose);
  py::dict res;
  res["u"] = SVD_out(0);
  res["d"] = SVD_out(1);
  res["v"] = SVD_out(2);
  return res;
}
py::dict HalkoSVD_full(arma::Mat<npdouble> &A, int dim, int iters = 5,
                       int seed = 0, int verbose = 1) {
  field<mat> SVD_out = ACTIONet::HalkoSVD(A, dim, iters, seed, verbose);
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
  field<mat> reduction = ACTIONet::reduce_kernel(S, reduced_dim, iters, seed,
                                                 SVD_algorithm, prenormalize, verbose);

  py::dict res;

  res["V"] = reduction(0);

  vec sigma = reduction(1).col(0);
  res["sigma"] = sigma;

  mat V = reduction(2);
  // printf("%d x %d\n", V.n_rows, V.n_cols);
  for (int i = 0; i < V.n_cols; i++) {
	  vec v = V.col(i) * sigma(i);
	  v = round(v*1e5)/1e5;
    double cs = sum(v);
    if( cs < 0)
		v = -v;
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
  field<mat> reduction = ACTIONet::reduce_kernel(S, reduced_dim, iters, seed,
                                                 SVD_algorithm, prenormalize, verbose);

  py::dict res;

  res["V"] = reduction(0);

  vec sigma = reduction(1).col(0);
  res["sigma"] = sigma;

  mat V = reduction(2);
  // printf("%d x %d\n", V.n_rows, V.n_cols);
  for (int i = 0; i < V.n_cols; i++) {
	  vec v = V.col(i) * sigma(i);
	  v = round(v*1e5)/1e5;
    double cs = sum(v);
    if( cs < 0)
		v = -v;
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
py::dict run_AA(arma::Mat<npdouble> &A, arma::Mat<npdouble> W0, int max_it = 50,
                double min_delta = 1e-16) {
  field<mat> AA_res = ACTIONet::run_AA(A, W0, max_it, min_delta);

  py::dict out;
  out["C"] = AA_res(0);
  out["W"] = A * AA_res(0);
  out["H"] = AA_res(1);

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
                    double min_delta = 1e-16) {
  ACTIONet::ACTION_results trace =
      ACTIONet::run_ACTION(S_r, k_min, k_max, thread_no, max_it, min_delta);

  py::dict res;

  py::list C(k_max);
  for (int i = 0; i < k_max; i++) {
	mat curr_C = trace.C[i + 1];
    //curr_C = clamp(curr_C, 1e-5, 1);
	//curr_C = normalise(curr_C, 1);
    C[i] = curr_C;
  }
  res["C"] = C;

  py::list H(k_max);
  for (int i = 0; i < k_max; i++) {
	mat curr_H = trace.H[i + 1];
    //curr_H = clamp(curr_H, 1e-5, 1);
	//curr_H = normalise(curr_H, 1);
    H[i] = curr_H;
  }
  res["H"] = H;

  return res;
}

// Runs multi-level ACTION decomposition method
//
// @param S_r Reduced kernel matrix
// @param k_min Minimum number of archetypes to consider (default=2)
// @param k_max Maximum number of archetypes to consider, or "depth" of
// decomposition (default=30)
// @param max_it,min_delta Convergence parameters for archetypal analysis
// @param max_trial Maximum number of trials before termination
//
// @return A named list with entries 'C' and 'H', each a list for different
// values of k
py::dict run_ACTION_plus(arma::Mat<npdouble> &S_r, int k_min = 2,
                         int k_max = 30, int max_it = 50,
                         double min_delta = 1e-16, int max_trial = 3) {
  ACTIONet::ACTION_results trace = ACTIONet::run_ACTION_plus(
      S_r, k_min, k_max, max_it, min_delta, max_trial);

  py::dict res;

  py::list C(trace.H.n_elem - 1);
  for (int i = k_min; i < trace.H.n_elem; i++) {
    C[i - 1] = trace.C[i];
  }
  res["C"] = C;

  py::list H(k_max);
  for (int i = 0; i < k_max; i++) {
    H[i - 1] = trace.H[i];
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
py::dict unify_archetypes(mat &S_r, mat &C_stacked, mat &H_stacked,
                                     double violation_threshold = 0.0,
                                     int thread_no = 0) {
  ACTIONet::unification_results results = ACTIONet::unify_archetypes(S_r, C_stacked, H_stacked, violation_threshold, thread_no);

  py::dict out_list;

  for (int i = 0; i < results.selected_archetypes.n_elem; i++)
    results.selected_archetypes[i]++;
  out_list["selected_archetypes"] = results.selected_archetypes;

  out_list["C_unified"] = sp_mat(results.C_unified);
  out_list["H_unified"] = sp_mat(results.H_unified);

  for (int i = 0; i < results.assigned_archetypes.n_elem; i++)
    results.assigned_archetypes[i]++;
  out_list["assigned_archetype"] = results.assigned_archetypes;

  out_list["ontology"] = results.dag_adj;
  out_list["ontology_node_attributes"] = results.dag_node_annotations;

  return out_list;
}

// Builds an interaction network from the multi-level archetypal decompositions
//
// @param H_stacked Output of the prune_archetypes() function.
// @param density Overall density of constructed graph. The higher the density,
// the more edges are retained (default = 1.0).
// @param thread_no Number of parallel threads (default = 0).
// @param mutual_edges_only Symmetrization strategy for nearest-neighbor edges.
// If it is true, only mutual-nearest-neighbors are returned (default=TRUE).
// @param distance_metric Distance metric to use: jsd, l2, ip
// @param nn_approach Nearest neighbor alogirthm: k*nn, knn
// @param k Optional parameter specifying k for knn algorithm 
//
// @return G Adjacency matrix of the ACTIONet graph.
arma::SpMat<npdouble> build_ACTIONet(arma::Mat<npdouble> &H_stacked,
                                     double density = 1.0, int thread_no = 0,
                                     bool mutual_edges_only = true,
				     string distance_metric="jsd",
				     string nn_approach="k*nn",
				     int k=10) {
  double M = 16, ef_construction = 200, ef = 50;
						 
  arma::SpMat<npdouble> G = ACTIONet::build_ACTIONet(H_stacked,density, thread_no, M, ef_construction, ef, mutual_edges_only, distance_metric, nn_approach, k);
  return G;
}

// Performs stochastic force-directed layout on the input graph (ACTIONet)
//
// @param G Adjacency matrix of the ACTIONet graph
// @param S_r Reduced kernel matrix (is used for reproducible initialization).
// @param compactness_level A value between 0-100, indicating the compactness of
// ACTIONet layout (default=50)
// @param layout_alg Algorithm to use for visualization layout (default=0).
// @param n_epochs Number of epochs for SGD algorithm (default=100).
// @param thread_no Number of threads.
//
// @return A named list \itemize{
// \item coordinates 2D coordinates of vertices.
// \item coordinates_3D 3D coordinates of vertices.
// \item colors De novo color of nodes inferred from their 3D embedding.
// }
py::dict layout_ACTIONet(arma::SpMat<npdouble> &G, arma::Mat<npdouble> S_r,
                         int compactness_level = 50, unsigned int n_epochs = 500, int layout_alg = 0, int thread_no = 0, int seed = 0) {
  field<arma::Mat<npdouble>> res =
      ACTIONet::layout_ACTIONet(G, S_r, compactness_level, n_epochs, layout_alg, thread_no, seed);

  py::dict out_list;
  out_list["coordinates"] = res(0);
  out_list["coordinates_3D"] = res(1);
  out_list["colors"] = res(2);

  return out_list;
}

// Encrypts a set of given input ids
//
// @param ids py::dict of input string ids
// @param pass Pass phrase to use for encryption
//
// @return A string array of encoded ids
vector<string> encode_ids(vector<string> ids, string pass) {
  vector<string> encoded_ids(ids.size());

  cryptor::set_key(pass);
  for (int i = 0; i < (int)ids.size(); i++) {
    auto enc = cryptor::encrypt(ids[i]);
    encoded_ids[i] = enc;
  }

  return encoded_ids;
}

// Decrypts a set of given encrypted ids
//
// @param encoded_ids py::dict of encrypted string ids
// @param pass Pass phrase to use for decryption
//
// @return A string array of decrypted ids
vector<string> decode_ids(vector<string> encoded_ids, string pass) {
  vector<string> decoded_ids(encoded_ids.size());

  cryptor::set_key(pass);
  for (int i = 0; i < (int)encoded_ids.size(); i++) {
    auto dec = cryptor::decrypt(encoded_ids[i]);
    decoded_ids[i] = dec;
  }

  return decoded_ids;
}

// Computes pseudobulk profiles per archetype
//
// @param S Input matrix (dense)
// @param H A soft membership matrix - Typically H_unified from the
// unify_architypes() function
//
// @return S matrix aggregated within each class of sample_assignments
arma::Mat<npdouble> compute_pseudo_bulk_per_archetype_full(
    arma::Mat<npdouble> &S, arma::Mat<npdouble> &H) {
  arma::Mat<npdouble> pb = ACTIONet::compute_pseudo_bulk_per_archetype(S, H);

  return pb;
}
arma::Mat<npdouble> compute_pseudo_bulk_per_archetype(arma::SpMat<npdouble> &S,
                                                      arma::Mat<npdouble> &H) {
  arma::Mat<npdouble> pb = ACTIONet::compute_pseudo_bulk_per_archetype(S, H);

  return pb;
}

// Computes pseudobulk profiles per cluster
//
// @param S Input matrix (dense)
// @param sample_assignments Any sample clustering/annotation (it has to be in
// {1, ..., max_class_num})
//
// @return S matrix aggregated within each class of sample_assignments
arma::Mat<npdouble> compute_pseudo_bulk_per_cluster_full(
    arma::Mat<npdouble> &S, uvec sample_assignments) {
  arma::Mat<npdouble> pb =
      ACTIONet::compute_pseudo_bulk_per_cluster(S, sample_assignments);

  return pb;
}
arma::Mat<npdouble> compute_pseudo_bulk_per_cluster(arma::SpMat<npdouble> &S,
                                                    uvec sample_assignments) {
  arma::Mat<npdouble> pb =
      ACTIONet::compute_pseudo_bulk_per_cluster(S, sample_assignments);

  return pb;
}

// Computes pseudobulk profiles (groups[k1] x individuals[k2])
//
// @param S Input matrix (dense)
// @param sample_assignments Any primary grouping - typically based on cell
// type/state (it has to be in {1, ..., k1})
// @param individuals Any Secondary grouping - typically corresponds to
// individuals (it has to be in {1, ..., k2})
//
// @return A list of pseudobulk profile, where each entry is matrix
// corresponding to one cell type/state
field<arma::Mat<npdouble>> compute_pseudo_bulk_per_cluster_and_ind_full(
    arma::Mat<npdouble> &S, uvec sample_assignments, uvec individuals) {
  field<arma::Mat<npdouble>> pbs_list =
      ACTIONet::compute_pseudo_bulk_per_cluster_and_ind(S, sample_assignments, individuals);

  return pbs_list;
}
field<arma::Mat<npdouble>> compute_pseudo_bulk_per_cluster_and_ind(arma::SpMat<npdouble> &S,
                                                       uvec sample_assignments,
                                                       uvec individuals) {
  field<arma::Mat<npdouble>> pbs_list =
      ACTIONet::compute_pseudo_bulk_per_cluster_and_ind(S, sample_assignments, individuals);

  return pbs_list;
}

// Renormalized input matrix to minimize differences in means
//
// @param S Input matrix (dense)
// @param sample_assignments Any primary grouping - typically based on cell
// type/state (it has to be in {1, ..., k1})
//
// @return A list with the first entry being the renormalized input matrix
arma::Mat<npdouble> renormalize_input_matrix_full(arma::Mat<npdouble> &S,
                                                  uvec sample_assignments) {
  arma::Mat<npdouble> S_norm =
      ACTIONet::renormalize_input_matrix(S, sample_assignments);

  return S_norm;
}
arma::SpMat<npdouble> renormalize_input_matrix(arma::SpMat<npdouble> &S,
                                               uvec sample_assignments) {
  arma::SpMat<npdouble> S_norm =
      ACTIONet::renormalize_input_matrix(S, sample_assignments);

  return S_norm;
}

// Compute feature specificity (from archetype footprints and binary input)
//
// @param S Input matrix (sparseMatrix - binary)
// @param H A soft membership matrix - Typically H_unified from the
// unify_archetypes() function.
//
// @return A list with the over/under-logPvals
py::dict compute_archetype_feature_specificity_bin(arma::SpMat<npdouble> &S,
                                                   arma::Mat<npdouble> &H) {
  field<arma::Mat<npdouble>> res =
      ACTIONet::compute_feature_specificity_bin(S, H);

  py::dict out_list;
  out_list["archetypes"] = res(0);
  out_list["upper_significance"] = res(1);
  out_list["lower_significance"] = res(2);

  return out_list;
}

// Compute feature specificity (discriminative scores)
//
// @param S Input matrix (dense)
// @param H A soft membership matrix - Typically H_unified from the
// unify_archetypes() function.
//
// @return A list with the over/under-logPvals
py::dict compute_archetype_feature_specificity_full(arma::Mat<npdouble> &S,
                                                    arma::Mat<npdouble> &H) {
  field<arma::Mat<npdouble>> res = ACTIONet::compute_feature_specificity(S, H);

  py::dict out_list;
  out_list["archetypes"] = res(0);
  out_list["upper_significance"] = res(1);
  out_list["lower_significance"] = res(2);

  return out_list;
}
py::dict compute_archetype_feature_specificity(arma::SpMat<npdouble> &S,
                                               arma::Mat<npdouble> &H) {
  field<arma::Mat<npdouble>> res = ACTIONet::compute_feature_specificity(S, H);

  py::dict out_list;
  out_list["archetypes"] = res(0);
  out_list["upper_significance"] = res(1);
  out_list["lower_significance"] = res(2);

  return out_list;
}

// Compute feature specificity (from cluster assignments)
//
// @param S Input matrix ("sparseMatrix")
// @param sample_assignments Vector of cluster assignments
//
// @return A list with the over/under-logPvals
py::dict compute_cluster_feature_specificity(arma::SpMat<npdouble> &S,
                                             arma::uvec &sample_assignments) {
  field<arma::Mat<npdouble>> res =
      ACTIONet::compute_feature_specificity(S, sample_assignments);

  py::dict out_list;
  out_list["archetypes"] = res(0);
  out_list["upper_significance"] = res(1);
  out_list["lower_significance"] = res(2);

  return out_list;
}
py::dict compute_cluster_feature_specificity_full(
    arma::Mat<npdouble> &S, arma::uvec &sample_assignments) {
  field<arma::Mat<npdouble>> res =
      ACTIONet::compute_feature_specificity(S, sample_assignments);

  py::dict out_list;
  out_list["archetypes"] = res(0);
  out_list["upper_significance"] = res(1);
  out_list["lower_significance"] = res(2);

  return out_list;
}

// Compute coreness of graph vertices
//
// @param G Input graph
//
// @return cn core-number of each graph node
uvec compute_core_number(sp_mat &G) {
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
arma::vec compute_archetype_core_centrality(arma::SpMat<npdouble> &G,
                                            arma::uvec sample_assignments) {
  vec conn = ACTIONet::compute_archetype_core_centrality(G, sample_assignments);

  return (conn);
}

//' OBSOLETE
//' Computes network diffusion over a given network, starting with an arbitrarty
// set of initial scores
//'
//' @param G Input graph
//' @param X0 Matrix of initial values per diffusion (ncol(G) == nrow(G) ==
// ncol(X0)) ' @param thread_no Number of parallel threads ' @param alpha
// Random-walk depth ( between [0, 1] ) ' @param max_it PageRank iterations
//'
//' @return Matrix of diffusion scores
arma::Mat<npdouble> compute_network_diffusion(arma::SpMat<npdouble> &G,
                                              arma::SpMat<npdouble> &X0,
                                              int thread_no = 0,
                                              double alpha = 0.85,
                                              int max_it = 3) {
  mat Diff = ACTIONet::compute_network_diffusion(G, X0, thread_no, alpha, max_it);

  return (Diff);
}

//' Quickly computes network diffusion over a given network, starting with an arbitrarty
//' return Matrix of diffusion scores
arma::Mat<npdouble> compute_network_diffusion_fast(arma::SpMat<npdouble> &G,
                                              arma::SpMat<npdouble> &X0,
                                              int thread_no = 0,
                                              double alpha = 0.85,
                                              int max_it = 5) {
  mat Diff = ACTIONet::compute_network_diffusion_fast(G, X0, thread_no, alpha, max_it);

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
arma::SpMat<npdouble> compute_sparse_network_diffusion(
    arma::SpMat<npdouble> &G, arma::SpMat<npdouble> &X0, double alpha = 0.85,
    double rho = 1e-4, double epsilon = 0.001, int max_iter = 20) {
  sp_mat scores = ACTIONet::compute_sparse_network_diffusion(G, X0, alpha, rho,
                                                             epsilon, max_iter);

  return (scores);
}

//' Computes feature enrichment wrt a given annotation
//'
//' @param scores Specificity scores of features
//' @param associations Binary matrix of annotations
//' @param L Length of the top-ranked scores to scan
//'
//' @return Matrix of log-pvalues
py::dict assess_enrichment(arma::Mat<npdouble> &scores,
                           arma::SpMat<npdouble> &associations,
                           int thread_no = 0) {
  field<mat> res = ACTIONet::assess_enrichment(scores, associations, thread_no);

  py::dict out_list;

  out_list["logPvals"] = res(0);
  out_list["thresholds"] = res(1);

  return (out_list);
}

// Clusters data points using the hierarchical DBSCAN algorithm.
//
// @param X Input data matrix with each row being a data point
//
// @return A list with \itemize{
// \item labels
// \item membershipProbabilities
// \item outlierScores
//}
py::dict run_HDBSCAN(arma::Mat<npdouble> &X, int min_points = 5,
                     int min_cluster_size = 5) {
  field<vec> res = ACTIONet::run_HDBSCAN(X, min_points, min_cluster_size);

  py::dict out_list;
  out_list["labels"] = res(0);
  out_list["membershipProbabilities"] = res(1);
  out_list["outlierScores"] = res(2);

  return (out_list);
}

// Computes the maximum-weight bipartite graph matching
//
// @param G Adjacency matrix of the input graph
//
// @return G_matched An adjacency matrix with a maximum of one nonzero entry on
// rows/columns
arma::Mat<npdouble> MWM_hungarian(mat &G) {
  mat G_matched = ACTIONet::MWM_hungarian(G);

  return G_matched;
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
vec signed_cluster(arma::SpMat<npdouble> &A, double resolution_parameter = 1.0,
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

arma::Mat<npdouble> unsigned_cluster_batch(arma::SpMat<npdouble> &A,
                                           vec resolutions,
                                           uvec initial_clusters_ = uvec(),
                                           int seed = 0) {
  uvec initial_clusters_uvec(A.n_rows);
  if (initial_clusters_.n_elem == A.n_rows) {
    for (int i = 0; i < A.n_rows; i++)
      initial_clusters_uvec(i) = initial_clusters_(i);
  } else {
    for (int i = 0; i < A.n_rows; i++) initial_clusters_uvec(i) = i;
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
vec unsigned_cluster(arma::SpMat<npdouble> &A,
                     double resolution_parameter = 1.0,
                     uvec initial_clusters_ = uvec(), int seed = 0) {
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

arma::Mat<npdouble> Prune_PageRank(arma::Mat<npdouble> &U,
                                   double density = 1.0) {
  mat G_matched = ACTIONet::Prune_PageRank(U, density);

  return G_matched;
}

//' Project a new data into current embedding
//'
//' @param W bipartite graph weights
//' @param coor2D, coor3D, colRGB Outpput of layout_ACTIONet()
//' @param Compactness compactness of the layout
//' @param n_epochs SGD iterts
//' @param thread_no # of threads
//'
//' @return Embedding(s)& colors
py::dict transform_layout(arma::SpMat<npdouble> &W, arma::Mat<npdouble> &coor2D,
                          Mat<npdouble> &coor3D, arma::Mat<npdouble> &colRGB,
                          int compactness_level = 50,
                          unsigned int n_epochs = 500, int thread_no = 0, int seed = 0) {
  field<mat> res = ACTIONet::transform_layout(
      W, coor2D, coor3D, colRGB, compactness_level, n_epochs, thread_no, seed);

  py::dict out_list;
  out_list["coordinates"] = res(0);
  out_list["coordinates_3D"] = res(1);
  out_list["colors"] = res(2);

  return out_list;
}

arma::Mat<npdouble> compute_full_sim(arma::Mat<npdouble> &H,
                                     int thread_no = 0) {
  arma::Mat<npdouble> G = ACTIONet::computeFullSim(H, thread_no);

  return G;
}

arma::vec run_LPA(sp_mat &G, arma::vec labels, double lambda = 1, int iters = 3, double sig_threshold = 3, arma::vec fixed_labels_ = arma::vec()) {
  arma::uvec fixed_labels_vec;
  if (!fixed_labels_.is_empty()) {
    arma::uvec fixed_labels_vec(fixed_labels_.size());
    for(int i = 0; i < fixed_labels_.size(); i++) {
		    fixed_labels_vec(i) = fixed_labels_(i);
    }
  }
  arma::vec new_labels = ACTIONet::LPA(G, labels, lambda, iters, sig_threshold, fixed_labels_vec);
  return(new_labels);
}

PYBIND11_MODULE(_ACTIONet, m) {
  m.doc() = R"pbdoc(
        ACTIONet package
        -----------------------

        .. currentmodule:: ACTIONet

        .. autosummary::
           :toctree: _generate

    )pbdoc";
  // SVD
  m.def("IRLB_SVD", &IRLB_SVD,
    "Computes SVD using IRLB algorithm.",
    py::arg("A"), py::arg("dim"), py::arg("iters") = 1000, py::arg("seed") = 0, py::arg("verbose") = 1);

  m.def("IRLB_SVD_full", &IRLB_SVD_full,
    "Computes SVD using IRLB algorithm.",
    py::arg("A"), py::arg("dim"), py::arg("iters") = 1000, py::arg("seed") = 0, py::arg("verbose") = 1);

  m.def("FengSVD", &FengSVD,
    "Computes SVD using Feng et al. algorithm.",
    py::arg("A"), py::arg("dim"), py::arg("iters") = 5, py::arg("seed") = 0, py::arg("verbose") = 1);

  m.def("FengSVD_full", &FengSVD_full,
    "Computes SVD using Feng et al. algorithm.",
    py::arg("A"), py::arg("dim"), py::arg("iters") = 5, py::arg("seed") = 0, py::arg("verbose") = 1);

  m.def("HalkoSVD", &HalkoSVD,
    "Computes SVD using Halko et al. algorithm.",
    py::arg("A"), py::arg("dim"), py::arg("iters") = 5, py::arg("seed") = 0, py::arg("verbose") = 1);

  m.def("HalkoSVD_full", &HalkoSVD_full,
    "Computes SVD using Halko et al. algorithm.",
    py::arg("A"), py::arg("dim"), py::arg("iters") = 5, py::arg("seed") = 0, py::arg("verbose") = 1);

  // Kernel reduction
  m.def("reduce_kernel", &reduce_kernel,
    "Computes reduced kernel matrix for a given profile",
    py::arg("S"), py::arg("reduced_dim") = 50, py::arg("iters") = 5,
    py::arg("seed") = 0, py::arg("SVD_algorithm") = 1, py::arg("prenormalize") = false, py::arg("verbose") = 1);

  m.def("reduce_kernel_full", &reduce_kernel_full,
    "Computes reduced kernel matrix for a given profile",
    py::arg("S"), py::arg("reduced_dim") = 50, py::arg("iters") = 5,
    py::arg("seed") = 0, py::arg("SVD_algorithm") = 1, py::arg("prenormalize") = false, py::arg("verbose") = 1);

  // Lower-level functions
  m.def("run_simplex_regression", &run_simplex_regression,
    "Solves min_{X} (|| AX - B ||) s.t. simplex constraint",
    py::arg("A"), py::arg("B"), py::arg("computeXtX") = false);

  m.def("run_AA", &run_AA,
    "Runs Archetypal Analysis (AA) Algorithm",
    py::arg("A"), py::arg("W0"), py::arg("max_it") = 50, py::arg("min_delta") = 0.01);

  m.def("run_SPA", &run_SPA,
    "Runs Successive Projection Algorithm (SPA) to solve separable NMF",
    py::arg("A"), py::arg("k"));

  m.def("run_SPA_rows_sparse", &run_SPA_rows_sparse,
    "Runs Successive Projection Algorithm (SPA) to solve separable NMF",
    py::arg("A"), py::arg("k"));

  m.def("renormalize_input_matrix", &renormalize_input_matrix,
    "Renormalized input matrix to minimize differences in means",
    py::arg("S"), py::arg("sample_assignments"));

  m.def("renormalize_input_matrix_full", &renormalize_input_matrix_full,
    "Renormalized input matrix to minimize differences in means",
    py::arg("S"), py::arg("sample_assignments"));

  // ACTION decomposition
  m.def("run_ACTION", &run_ACTION,
    "Runs multi-level ACTION decomposition method",
    py::arg("S_r"), py::arg("k_min") = 2, py::arg("k_max") = 30, py::arg("thread_no") = 0,
    py::arg("max_it") = 50, py::arg("min_delta") = 0.01);

  m.def("run_ACTION_plus", &run_ACTION_plus,
    "Runs multi-level ACTION decomposition method",
    py::arg("S_r"), py::arg("k_min") = 2, py::arg("k_max") = 30,
    py::arg("max_it") = 50, py::arg("min_delta") = 0.01, py::arg("max_trial") = 3);

  // Archetypes
  m.def("prune_archetypes", &prune_archetypes,
    "Filters multi-level archetypes and concatenate filtered archetypes",
    py::arg("C_trace"), py::arg("H_trace"), py::arg("min_specificity_z_threshold") = -3,
    py::arg("min_cells") = 3);

  m.def("unify_archetypes", &unify_archetypes,
    "Identifies and aggregates redundant archetypes into equivalent classes",
    py::arg("S_r"), py::arg("C_stacked"), py::arg("H_stacked"),
    py::arg("violation_threshold") = 0.0, py::arg("thread_no") = 0);

  // Network
  m.def("build_ACTIONet", &build_ACTIONet,
        "Builds an interaction network from the multi-level archetypal decompositions",
        py::arg("H_stacked"), py::arg("density") = 1.0,
        py::arg("thread_no") = 0, py::arg("mutual_edges_only") = true,
        py::arg("distance_metric")="jsd",py::arg("nn_approach")="k*nn",py::arg("k"));

  
  m.def("layout_ACTIONet", &layout_ACTIONet,
    "Performs stochastic force-directed layout on the input graph (ACTIONet)",
    py::arg("G"), py::arg("S_r"), py::arg("compactness_level") = 50, py::arg("n_epochs") = 500, py::arg("layout_alg") = 0, py::arg("thread_no") = 0, py::arg("seed") = 0);

  // Pseudobulk
  m.def("compute_pseudo_bulk_per_archetype", &compute_pseudo_bulk_per_archetype,
    "Computes pseudobulk profiles",
    py::arg("S"), py::arg("sample_assignments"));

  m.def("compute_pseudo_bulk_per_archetype_full", &compute_pseudo_bulk_per_archetype_full,
    "Computes pseudobulk profiles",
    py::arg("S"), py::arg("sample_assignments"));

  m.def("compute_pseudo_bulk_per_cluster", &compute_pseudo_bulk_per_cluster,
    "Computes pseudobulk profiles",
    py::arg("S"), py::arg("sample_assignments"));

  m.def("compute_pseudo_bulk_per_cluster_full", &compute_pseudo_bulk_per_cluster_full,
    "Computes pseudobulk profiles",
    py::arg("S"), py::arg("sample_assignments"));

  m.def("compute_pseudo_bulk_per_cluster_and_ind", &compute_pseudo_bulk_per_cluster_and_ind,
    "Computes pseudobulk profiles (groups[k1] x individuals[k2])",
    py::arg("S"), py::arg("sample_assignments"), py::arg("individuals"));

  m.def("compute_pseudo_bulk_per_cluster_ind_full", &compute_pseudo_bulk_per_cluster_and_ind_full,
    "Computes pseudobulk profiles (groups[k1] x individuals[k2])",
    py::arg("S"), py::arg("sample_assignments"), py::arg("individuals"));

  m.def("renormalize_input_matrix", &renormalize_input_matrix,
    "Renormalized input matrix to minimize differences in means",
    py::arg("S"), py::arg("sample_assignments"));

  m.def("renormalize_input_matrix_full", &renormalize_input_matrix_full,
    "Renormalized input matrix to minimize differences in means",
    py::arg("S"), py::arg("sample_assignments"));

  // Feature specificity
  m.def("compute_archetype_feature_specificity_bin", &compute_archetype_feature_specificity_bin,
    "Compute feature specificity (from archetype footprints and binary input)",
    py::arg("S"), py::arg("H"));

  m.def("compute_archetype_feature_specificity", &compute_archetype_feature_specificity,
    "Compute feature specificity (discriminative scores)",
    py::arg("S"), py::arg("sample_assignments"));

  m.def("compute_archetype_feature_specificity_full", &compute_archetype_feature_specificity_full,
    "Compute feature specificity (discriminative scores)",
    py::arg("S"), py::arg("sample_assignments"));

  m.def("compute_cluster_feature_specificity", &compute_cluster_feature_specificity,
    "Compute feature specificity (discriminative scores)",
    py::arg("S"), py::arg("sample_assignments"));

  m.def("compute_cluster_feature_specificity_full", &compute_cluster_feature_specificity_full,
    "Compute feature specificity (discriminative scores)",
    py::arg("S"), py::arg("sample_assignments"));

  m.def("compute_core_number", &compute_core_number,
    "Compute coreness of graph vertices",
    py::arg("G"));

  m.def("compute_archetype_core_centrality", &compute_archetype_core_centrality,
    "Compute the overall connectivity of each node",
    py::arg("G"), py::arg("sample_assignments"));

  m.def("compute_network_diffusion", &compute_network_diffusion,
    "Computes PageRank for a selected set of nodes",
    py::arg("G"), py::arg("X0"), py::arg("thread_no") = 0,
    py::arg("alpha") = 0.85, py::arg("max_it") = 3);

  m.def("compute_network_diffusion_fast", &compute_network_diffusion_fast,
    "Computes PageRank for a selected set of nodes",
    py::arg("G"), py::arg("X0"), py::arg("thread_no") = 0,
    py::arg("alpha") = 0.85, py::arg("max_it") = 5);

  m.def("compute_sparse_network_diffusion", &compute_sparse_network_diffusion,
    "Computes L1-regularized PageRank for a selected set of nodes",
    py::arg("G"), py::arg("X0"), py::arg("alpha") = 0.85,
    py::arg("rho") = 1e-4, py::arg("epsilon") = 0.001,
    py::arg("max_iter") = 20);

  m.def("assess_enrichment", &assess_enrichment,
    "Performs enrichment analysis",
    py::arg("scores"), py::arg("associations"), py::arg("thread_no") = 0);

  m.def("encode_ids", &encode_ids,
    "Encrypts a set of given input ids",
    py::arg("ids"), py::arg("pass"));

  m.def("decode_ids", &decode_ids,
    "Decrypts a set of given encrypted ids",
    py::arg("encoded_ids"), py::arg("pass"));

  m.def("run_HDBSCAN", &run_HDBSCAN,
    "Clusters data points using the hierarchical DBSCAN algorithm",
    py::arg("X"), py::arg("min_points"), py::arg("min_cluster_size"));

  m.def("MWM_hungarian", &MWM_hungarian,
    "Computes the maximum-weight bipartite graph matching",
    py::arg("G"));

  m.def("signed_cluster", &unsigned_cluster,
    "Computes graph clustering using Leiden algorith over signed graphs",
    py::arg("A"), py::arg("resolution_parameter") = 1.0,
    py::arg("initial_clusters") = uvec(), py::arg("seed") = 0);

  m.def("unsigned_cluster_batch", &unsigned_cluster_batch,
    "Computes graph clustering using Leiden algorith over signed graphs",
    py::arg("A"), py::arg("resolution_parameter") = 1.0,
    py::arg("initial_clusters") = uvec(), py::arg("seed") = 0);

  m.def("unsigned_cluster", &unsigned_cluster,
    "Computes graph clustering using Leiden algorith over unsigned graphs",
    py::arg("A"), py::arg("resolution_parameter") = 1.0,
    py::arg("initial_clusters") = uvec(), py::arg("seed") = 0);

  m.def("Prune_PageRank", &Prune_PageRank,
    "",
    py::arg("U"), py::arg("density"));

  m.def("transform_layout", &transform_layout,
    "Project a new data into current embedding",
    py::arg("W"), py::arg("coor2D"), py::arg("coor3D"), py::arg("colRGB"),
    py::arg("compactness_level") = 50, py::arg("n_epochs") = 500,
    py::arg("thread_no") = 0, py::arg("seed") = 0);

  m.def("compute_full_sim", &compute_full_sim,
    "",
    py::arg("H"), py::arg("thread_no") = 0);

  m.def("run_LPA", &run_LPA,
    "Run label prepagation on a given set of known labels",
    py::arg("G"), py::arg("labels"), py::arg("lambda") = 1, py::arg("iters") = 3,
    py::arg("sig_threshold") = 3, py::arg("fixed_labels_") = arma::vec());

#ifdef VERSION_INFO
  m.attr("__version__") = VERSION_INFO;
#else
  m.attr("__version__") = "dev";
#endif
}
