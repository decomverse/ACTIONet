#include <stdint.h>
#include <utility>
#include <numpy/npy_common.h>
#include <functional>
#include <arma_wrapper.h>
#include <pybind11/stl.h>
#include <ACTIONet.h>

using aw::npint;
using aw::npdouble;
using aw::intmat;
using aw::dmat;
using aw::dcube; 
using aw::dvec; 

namespace py = pybind11;
using namespace py::literals;

// Computes reduced kernel matrix for a given (single-cell) profile
//
// @param S Input matrix (dense)
// @param reduced_dim Dimension of the reduced kernel matrix (default=50)
// @param iters Number of SVD iterations (default=5)
// @param seed Random seed (default=0)
// @param reduction_algorithm Kernel reduction algorithm. Currently only ACTION method (1) is implemented (default=1)
// @param SVD_algorithm SVD algorithm to use. Currently supported methods are Halko (1) and Feng (2) (default=1)
// 
// @return A named list with S_r, V, lambda, and exp_var. \itemize{
// \item S_r: reduced kernel matrix of size reduced_dim x #samples.
// \item V: Associated left singular-vectors (useful for reconstructing discriminative scores for features, such as genes).
// \item lambda, exp_var: Summary statistics of the sigular-values.
// }
py::dict reduce_kernel_full(arma::Mat<npdouble> &S, int reduced_dim = 50, int iters = 5, int seed = 0, int SVD_algorithm = 0, bool prenormalize = false) {
	
	field<mat> reduction = ACTIONet::reduce_kernel(S, reduced_dim, iters, seed, SVD_algorithm, prenormalize);				


    py::dict res;

	res["V"] = reduction(0);
	
	vec sigma = reduction(1).col(0);
	res["sigma"] = sigma;
	
	mat V = reduction(2);
	for(int i = 0; i < V.n_cols; i++) {
		V.col(i) *= sigma(i);
	}	
	res["S_r"] = trans(V);
		

	res["A"] = reduction(3);
	res["B"] = reduction(4);
	
		
	return res;
}
py::dict reduce_kernel(arma::SpMat<npdouble> &S, int reduced_dim = 50, int iters = 5, int seed = 0, int SVD_algorithm = 0, bool prenormalize = false) {
	
	field<mat> reduction = ACTIONet::reduce_kernel(S, reduced_dim, iters, seed, SVD_algorithm, prenormalize);				


    py::dict res;

	res["V"] = reduction(0);
	
	vec sigma = reduction(1).col(0);
	res["sigma"] = sigma;
	
	mat V = reduction(2);
	for(int i = 0; i < V.n_cols; i++) {
		V.col(i) *= sigma(i);
	}	
	res["S_r"] = trans(V);
		

	res["A"] = reduction(3);
	res["B"] = reduction(4);
	
		
	return res;
}




// Solves min_{X} (|| AX - B ||) s.t. simplex constraint
//
// @param A Input matrix (dense)
// @param B Input matrix (dense)
// 
// @return X Solution
arma::Mat<npdouble> run_simplex_regression(arma::Mat<npdouble>& A, arma::Mat<npdouble>& B, bool computeXtX = false) {	
	arma::Mat<npdouble> X = ACTIONet::run_simplex_regression(A, B, computeXtX);
	
	return X;
}


// Runs Successive Projection Algorithm (SPA) to solve separable NMF
//
// @param A Input matrix (dense)
// @param k Number of columns to select
// 
// @return A named list with entries 'selected_columns' and 'norms'
py::dict run_SPA(arma::Mat<npdouble>& A, int k) {	
	
	ACTIONet::SPA_results res = ACTIONet::run_SPA(A, k);
	uvec selected_columns = res.selected_columns;
	
	vec cols(k);
	for(int i = 0; i < k; i++) {
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
py::dict run_AA(arma::Mat<npdouble>&A, arma::Mat<npdouble> W0, int max_it = 50, double min_delta = 0.01) {	

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
// @param k_max Maximum number of archetypes to consider, or "depth" of decomposition (default=30)
// @param thread_no Number of parallel threads (default=4)
// 
// @return A named list with entries 'C' and 'H', each a list for different values of k
py::dict run_ACTION(arma::Mat<npdouble>& S_r, int k_min = 2, int k_max=30, int thread_no = 4, int max_it = 50, double min_delta = 0.01) {	

	ACTIONet::ACTION_results trace = ACTIONet::run_ACTION(S_r, k_min, k_max, thread_no, max_it, min_delta);

	py::dict res;
	
    py::list C(k_max);
	for (int i = 0; i < k_max; i++) {
		C[i] = trace.C[i+1];
	}
	res["C"] = C;	

    py::list H(k_max);
	for (int i = 0; i < k_max; i++) {
		H[i] = trace.H[i+1];
	}
	res["H"] = H;
	
		
	return res;
}


// Filters multi-level archetypes and concatenate filtered archetypes.
// (Pre-ACTIONet archetype processing)
//
// @param C_trace,H_trace Output of ACTION
// @param min_specificity_z_threshold Defines the stringency of pruning nonspecific archetypes. 
// The larger the value, the more archetypes will be filtered out (default=-1)
// 
// @return A named list: \itemize{
// \item selected_archs: py::dict of final archetypes that passed the filtering/pruning step.
// \item C_stacked,H_stacked: Horizontal/Vertical concatenation of filtered C and H matrices, respectively.
// }
py::dict prune_archetypes(vector<arma::Mat<npdouble>> &C_trace_list, vector<arma::Mat<npdouble>> & H_trace_list, double min_specificity_z_threshold = -1, int min_cells = 3) {	

	field<arma::Mat<npdouble>> C_trace(C_trace_list.size());
	field<arma::Mat<npdouble>> H_trace(H_trace_list.size());	
	for(int i = 0; i < C_trace_list.size(); i++) {
		C_trace[i] = C_trace_list[i]; //aw::conv_to<arma::Mat<npdouble>>::from(C_trace_list[i]);
		H_trace[i] = H_trace_list[i]; //aw::conv_to<arma::Mat<npdouble>>::from(H_trace_list[i]);
	}
	
	ACTIONet::multilevel_archetypal_decomposition results = ACTIONet::prune_archetypes(C_trace, H_trace, min_specificity_z_threshold, min_cells);
	
		
	py::dict out_list;		

	// It is 0-indexed, but 1-index makes more sense for assignments
	for(int i = 0; i < (int)results.selected_archs.n_elem; i++) results.selected_archs[i]++;
		out_list["selected_archs"] = results.selected_archs;
	
	out_list["C_stacked"] = results.C_stacked;
	out_list["H_stacked"] = results.H_stacked;


    return out_list;
}


// Identifies and aggregates redundant archetypes into equivalent classes
// (Post-ACTIONet archetype processing)
//
// @param G Adjacency matrix of the ACTIONet graph
// @param S_r Reduced kernel profile
// @param C_stacked,H_stacked Output of reconstruct_archetypes()
// 
// @return A named list: \itemize{
// \item archetype_groups: Equivalent classes of archetypes (non-redundant)
// \item C_unified,H_unified: C and H matrices of unified archetypes
// \item sample_assignments: Assignment of samples/cells to unified archetypes
// }
py::dict unify_archetypes(mat &S_r, mat &C_stacked, mat &H_stacked, double min_edge_weight = 0.5, int min_coreness = 2, double resolution = 1.0, int min_repeat = 2, int thread_no = 0, double alpha = 0.05, double beta = 0.5, double outlier_threshold = 0.5, int minPoints = 5, int minClusterSize = 5, double cond_threshold = 10, int normalization_type = 0, bool preprocess_adj = true, bool reduce_G = false, int method_type = 3) {
	
	ACTIONet::unification_results results = ACTIONet::unify_archetypes(S_r, C_stacked, H_stacked, min_edge_weight, min_coreness, resolution, min_repeat, thread_no, alpha, beta, outlier_threshold, minPoints, minClusterSize, cond_threshold, normalization_type, preprocess_adj, reduce_G, method_type);
	
		
	py::dict out_list;		

	for(int i = 0; i < results.selected_archetypes.n_elem; i++) results.selected_archetypes[i]++;
	out_list["selected_archetypes"] = results.selected_archetypes;


	out_list["C_unified"] = results.C_unified;
	out_list["H_unified"] = results.H_unified;
	
	for(int i = 0; i < results.assigned_archetypes.n_elem; i++) results.assigned_archetypes[i]++;
	out_list["assigned_archetypes"] = results.assigned_archetypes;

    return out_list;
}



// Builds an interaction network from the multi-level archetypal decompositions
//
// @param H_stacked Output of the prune_archetypes() function.
// @param density Overall density of constructed graph. The higher the density, the more edges are retained (default = 1.0).
// @param thread_no Number of parallel threads (default = 4).
// @param mutual_edges_only Symmetrization strategy for nearest-neighbor edges. 
// If it is true, only mutual-nearest-neighbors are returned (default=TRUE).
// 
// @return G Adjacency matrix of the ACTIONet graph.
arma::SpMat<npdouble> build_ACTIONet(arma::Mat<npdouble>& H_stacked, double density = 1.0, int thread_no=8, bool mutual_edges_only = true) {
	double M = 16, ef_construction = 200, ef = 50;
	
	const arma::SpMat<npdouble> G = ACTIONet::build_ACTIONet(H_stacked, density, thread_no, M, ef_construction, ef, mutual_edges_only);
	
    return G;
}


// Performs stochastic force-directed layout on the input graph (ACTIONet)
//
// @param G Adjacency matrix of the ACTIONet graph
// @param S_r Reduced kernel matrix (is used for reproducible initialization).
// @param compactness_level A value between 0-100, indicating the compactness of ACTIONet layout (default=50)
// @param n_epochs Number of epochs for SGD algorithm (default=100).
// @param thread_no Number of threads.
// 
// @return A named list \itemize{
// \item coordinates 2D coordinates of vertices.
// \item coordinates_3D 3D coordinates of vertices.
// \item colors De novo color of nodes inferred from their 3D embedding.
// }
py::dict layout_ACTIONet(arma::SpMat<npdouble>& G, arma::Mat<npdouble>& S_r, int compactness_level, unsigned int n_epochs, int thread_no) {

	field<arma::Mat<npdouble>> res = ACTIONet::layout_ACTIONet(G, S_r, compactness_level, n_epochs, thread_no);
    
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
    for(int i = 0; i < (int)ids.size(); i++) {
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
    for(int i = 0; i < (int)encoded_ids.size(); i++) {
		auto dec = cryptor::decrypt(encoded_ids[i]);
		decoded_ids[i] = dec;
	}
	
    return decoded_ids;
}



// Computes pseudobulk profiles
//
// @param S Input matrix (dense)
// @param sample_assignments Any sample clustering/annotation (it has to be in {1, ..., max_class_num})
// 
// @return S matrix aggregated within each class of sample_assignments
arma::Mat<npdouble> compute_pseudo_bulk_full(arma::Mat<npdouble> &S, uvec sample_assignments) {
	
	arma::Mat<npdouble> pb = ACTIONet::compute_pseudo_bulk_per_cluster(S, sample_assignments);
	
    return pb;
}
arma::Mat<npdouble> compute_pseudo_bulk(arma::SpMat<npdouble> &S, uvec sample_assignments) {
	
	arma::Mat<npdouble> pb = ACTIONet::compute_pseudo_bulk_per_cluster(S, sample_assignments);
	
    return pb;
}


// Computes pseudobulk profiles (groups[k1] x individuals[k2])
//
// @param S Input matrix (dense)
// @param sample_assignments Any primary grouping - typically based on cell type/state (it has to be in {1, ..., k1})
// @param individuals Any Secondary grouping - typically corresponds to individuals (it has to be in {1, ..., k2})
// 
// @return A list of pseudobulk profile, where each entry is matrix corresponding to one cell type/state
field<arma::Mat<npdouble>> compute_pseudo_bulk_per_ind_full(arma::Mat<npdouble>& S, uvec sample_assignments, uvec individuals) {

	field<arma::Mat<npdouble>> pbs_list = ACTIONet::compute_pseudo_bulk_per_ind(S, sample_assignments, individuals);
	
    return pbs_list;
}
field<arma::Mat<npdouble>> compute_pseudo_bulk_per_ind(arma::SpMat<npdouble>& S, uvec sample_assignments, uvec individuals) {

	field<arma::Mat<npdouble>> pbs_list = ACTIONet::compute_pseudo_bulk_per_ind(S, sample_assignments, individuals);
	
    return pbs_list;
}


// Renormalized input matrix to minimize differences in means
//
// @param S Input matrix (dense)
// @param sample_assignments Any primary grouping - typically based on cell type/state (it has to be in {1, ..., k1})
// 
// @return A list with the first entry being the renormalized input matrix
arma::Mat<npdouble> renormalize_input_matrix_full(arma::Mat<npdouble>& S, uvec sample_assignments) {
	arma::Mat<npdouble> S_norm = ACTIONet::renormalize_input_matrix(S, sample_assignments);

	return(S_norm);
}
arma::SpMat<npdouble> renormalize_input_matrix(arma::SpMat<npdouble>& S, uvec sample_assignments) {
	arma::SpMat<npdouble> S_norm = ACTIONet::renormalize_input_matrix(S, sample_assignments);

	return(S_norm);
}


// Compute feature specificity (discriminative scores)
//
// @param S Input matrix (dense)
// @param H A soft membership matrix - Typically H_unified from the unify_archetypes() function.
// 
// @return A list with the over/under-logPvals
py::dict compute_archetype_feature_specificity_full(arma::Mat<npdouble> &S, arma::Mat<npdouble>& H) {
    
	field<arma::Mat<npdouble>> res = ACTIONet::compute_feature_specificity(S, H);

    py::dict out_list;
	out_list["archetypes"] = res(0);
	out_list["upper_significance"] = res(1);
	out_list["lower_significance"] = res(2);


	return(out_list);
}
py::dict compute_archetype_feature_specificity(arma::SpMat<npdouble> &S, arma::Mat<npdouble>& H) {
    
	field<arma::Mat<npdouble>> res = ACTIONet::compute_feature_specificity(S, H);

    py::dict out_list;
	out_list["archetypes"] = res(0);
	out_list["upper_significance"] = res(1);
	out_list["lower_significance"] = res(2);


	return(out_list);
}


py::dict compute_archetype_feature_specificity_full(arma::Mat<npdouble> &S, arma::uvec sample_assignments) {
    
	field<arma::Mat<npdouble>> res = ACTIONet::compute_feature_specificity(S, sample_assignments);

    py::dict out_list;
	out_list["archetypes"] = res(0);
	out_list["upper_significance"] = res(1);
	out_list["lower_significance"] = res(2);


	return(out_list);
}
py::dict compute_cluster_feature_specificity(arma::SpMat<npdouble> &S, arma::uvec sample_assignments) {
    
	field<arma::Mat<npdouble>> res = ACTIONet::compute_feature_specificity(S, sample_assignments);

    py::dict out_list;
	out_list["archetypes"] = res(0);
	out_list["upper_significance"] = res(1);
	out_list["lower_significance"] = res(2);


	return(out_list);
}


//' Compute coreness of subgraph vertices induced by each archetype
//'
//' @param G Input graph
//' @param sample_assignments Archetype discretization (output of unify_archetypes())
//' 
//' @return cn core-number of each graph node
arma::vec compute_archetype_core_centrality(arma::SpMat<npdouble> &G, arma::uvec sample_assignments) {

	vec conn = ACTIONet::compute_archetype_core_centrality(G, sample_assignments);

	return(conn);	
}	


//' Computes network diffusion over a given network, starting with an arbitrarty set of initial scores
//'
//' @param G Input graph
//' @param X0 Matrix of initial values per diffusion (ncol(G) == nrow(G) == ncol(X0))
//' @param thread_no Number of parallel threads
//' @param alpha Random-walk depth ( between [0, 1] )
//' @param max_it PageRank iterations
//' 
//' @return Matrix of diffusion scores
arma::Mat<npdouble> compute_network_diffusion(arma::SpMat<npdouble> &G, arma::SpMat<npdouble> &X0, int thread_no = 4, double alpha = 0.85, int max_it = 3) {

	mat Diff = ACTIONet::compute_network_diffusion(G, X0, thread_no, alpha, max_it);

	return(Diff);
}




//' Computes sparse network diffusion over a given network, starting with an arbitrarty set of initial scores
//'
//' @param G Input graph
//' @param X0 Matrix of initial values per diffusion (ncol(G) == nrow(G) == ncol(X0))
//' @param alpha Random-walk depth ( between [0, 1] )
//' @param rho Sparsity controling parameter
//' @param epsilon,max_it Conditions on the length of diffusion 
//' 
//' @return Matrix of sparse diffusion scores
arma::SpMat<npdouble> compute_sparse_network_diffusion(arma::SpMat<npdouble> &G, arma::SpMat<npdouble> &X0, double alpha = 0.85, double rho = 1e-4, double epsilon = 0.001, int max_iter = 20) {		

	sp_mat scores = ACTIONet::compute_sparse_network_diffusion(G, X0, alpha, rho, epsilon, max_iter);

	return(scores);
}



//' Computes feature enrichment wrt a given annotation
//'
//' @param scores Specificity scores of features
//' @param associations Binary matrix of annotations
//' @param L Length of the top-ranked scores to scan
//' 
//' @return Matrix of log-pvalues
py::dict assess_enrichment(arma::Mat<npdouble> &scores, arma::SpMat<npdouble> &associations, int thread_no = 0) {
	
	field<mat> res = ACTIONet::assess_enrichment(scores, associations, thread_no);

	py::dict out_list;	
		
	out_list["logPvals"] = res(0);
	out_list["thresholds"] = res(1);

	return(out_list);
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
py::dict transform_layout(arma::SpMat<npdouble> &W, arma::Mat<npdouble>& coor2D, Mat<npdouble>& coor3D, arma::Mat<npdouble>& colRGB, int compactness_level = 50, unsigned int n_epochs = 500, int thread_no = 0) {

	field<mat> res = ACTIONet::transform_layout(W, coor2D, coor3D, colRGB, compactness_level, n_epochs, thread_no);
    
	py::dict out_list;		
	out_list["coordinates"] = res(0);
	out_list["coordinates_3D"] = res(1);
	out_list["colors"] = res(2);

    return out_list;
    	
}

//' Computes graph clustering using Leiden algorith over signed graphs
//'
//' @param G Adjacency matrix of the input graph
//' @param resolution_parameter Granularity of clustering. Larger values result in more clusters (default = 1.0)
//' @param initial_clusters_ Initialization vector for clusters (if available)
//' @param seed Random seed
//' 
//' @return clusters Assignment vector of samples to clusters
vec signed_cluster(arma::SpMat<npdouble> &A, double resolution_parameter = 1.0, uvec initial_clusters = uvec(), int seed = 0) {
	vec clusters = ACTIONet::signed_cluster(A, resolution_parameter, initial_clusters, seed);

    return clusters;	
}


//' Computes graph clustering using Leiden algorith over unsigned graphs
//'
//' @param G Adjacency matrix of the input graph
//' @param resolution_parameter Granularity of clustering. Larger values result in more clusters (default = 1.0)
//' @param initial_clusters_ Initialization vector for clusters (if available)
//' @param seed Random seed
//' 
//' @return clusters Assignment vector of samples to clusters
vec unsigned_cluster(arma::SpMat<npdouble> &A, double resolution_parameter = 1.0, uvec initial_clusters = uvec(), int seed = 0) {
	vec clusters = ACTIONet::unsigned_cluster(A, resolution_parameter, initial_clusters, seed);

    return clusters;	
}

PYBIND11_MODULE(_ACTIONet, m) {
    m.doc() = R"pbdoc(
        ACTIONet package
        -----------------------

        .. currentmodule:: ACTIONet

        .. autosummary::
           :toctree: _generate

    )pbdoc";
    	
		m.def("reduce_kernel", py::overload_cast<arma::SpMat<npdouble> &, int, int, int, int, bool>(&reduce_kernel), "Computes reduced kernel matrix for a given (single-cell) profile", py::arg("S"), py::arg("reduced_dim")=50, py::arg("iters")=5, py::arg("seed")=0, py::arg("SVD_algorithm")=1, py::arg("prenormalize") = false);		
		m.def("reduce_kernel_full", py::overload_cast<arma::Mat<npdouble> &, int, int, int, int, bool>(&reduce_kernel_full), "Computes reduced kernel matrix for a given (single-cell) profile", py::arg("S"), py::arg("reduced_dim")=50, py::arg("iters")=5, py::arg("seed")=0, py::arg("SVD_algorithm")=1, py::arg("prenormalize") = false);		



		m.def("run_simplex_regression", py::overload_cast<arma::Mat<npdouble>&, arma::Mat<npdouble>&, bool>(&run_simplex_regression), "Solves min_{X} (|| AX - B ||) s.t. simplex constraint", py::arg("A"), py::arg("B"), py::arg("computeXtX")=false);
		
		
		
		m.def("run_SPA", py::overload_cast<arma::Mat<npdouble>&, int>(&run_SPA), "Runs Successive Projection Algorithm (SPA) to solve separable NMF", py::arg("A"), py::arg("k"));
		m.def("run_AA", py::overload_cast<arma::Mat<npdouble>&, arma::Mat<npdouble>, int, double>(&run_AA), "Runs Archetypal Analysis (AA) Algorithm", py::arg("A"), py::arg("W0"), py::arg("max_it") = 50, py::arg("min_delta") = 0.01);
		m.def("run_ACTION", py::overload_cast<arma::Mat<npdouble>&, int, int, int, int, double>(&run_ACTION), "Runs multi-level ACTION decomposition method", py::arg("S_r"), py::arg("k_min")=2, py::arg("k_max")=30, py::arg("thread_no")=4, py::arg("max_it")=50, py::arg("min_delta")=0.01);


		m.def("prune_archetypes", py::overload_cast<vector<arma::Mat<npdouble>>&, vector<arma::Mat<npdouble>> &, double, int>(&prune_archetypes), "Filters multi-level archetypes and concatenate filtered archetypes", py::arg("C_trace"), py::arg("H_trace"), py::arg("min_specificity_z_threshold")=-1, py::arg("int min_cells") = 3);		
		m.def("unify_archetypes", py::overload_cast<arma::Mat<npdouble>&, arma::Mat<npdouble>&, arma::Mat<npdouble>&, double, int, double, int, int, double, double, double, int, int, double, int, bool, bool, int>(&unify_archetypes), "Identifies and aggregates redundant archetypes into equivalent classes", py::arg("S_r"), py::arg("C_stacked"), py::arg("H_stacked"), py::arg("min_edge_weight") = 0.5, py::arg("min_coreness") = 2, py::arg("resolution") = 1.0, py::arg("min_repeat") = 2, py::arg("thread_no") = 0, py::arg("alpha") = 0.05, py::arg("beta") = 0.5, py::arg("outlier_threshold") = 0.5, py::arg("minPoints") = 5, py::arg("minClusterSize") = 5, py::arg("cond_threshold") = 10, py::arg("normalization_type") = 0, py::arg("preprocess_adj") = true, py::arg("reduce_G") = false, py::arg("method_type") = 3);
	
		m.def("build_ACTIONet", py::overload_cast<arma::Mat<npdouble>&, double, int, bool>(&build_ACTIONet), "Builds an interaction network from the multi-level archetypal decompositions", py::arg("H_stacked"), py::arg("density") = 1.0, py::arg("thread_no")=8, py::arg("mutual_edges_only")=true);		
		m.def("layout_ACTIONet", py::overload_cast<arma::SpMat<npdouble>&, arma::Mat<npdouble>&, int, unsigned int, int>(&layout_ACTIONet), "Performs stochastic force-directed layout on the input graph (ACTIONet)", py::arg("G"), py::arg("S_r"), py::arg("compactness_level")= 50, py::arg("n_epochs")=100, py::arg("thread_no")= 4);

		m.def("compute_pseudo_bulk_full", py::overload_cast<arma::Mat<npdouble>&, uvec>(&compute_pseudo_bulk_full), "Computes pseudobulk profiles", py::arg("S"), py::arg("sample_assignments"));
		m.def("compute_pseudo_bulk_per_ind_full", py::overload_cast<arma::Mat<npdouble>&, uvec, uvec>(&compute_pseudo_bulk_per_ind_full), "Computes pseudobulk profiles (groups[k1] x individuals[k2])", py::arg("S"), py::arg("sample_assignments"), py::arg("individuals"));
		
		m.def("renormalize_input_matrix_full", py::overload_cast<arma::Mat<npdouble>&, uvec>(&renormalize_input_matrix_full), "Renormalized input matrix to minimize differences in means", py::arg("S"), py::arg("sample_assignments"));
		m.def("compute_archetype_feature_specificity_full", py::overload_cast<arma::Mat<npdouble> &, arma::Mat<npdouble>&>(&compute_archetype_feature_specificity_full), "Compute feature specificity (discriminative scores)", py::arg("S"), py::arg("H"));
		m.def("compute_archetype_feature_specificity_full", py::overload_cast<arma::Mat<npdouble> &, uvec>(&compute_archetype_feature_specificity_full), "Compute feature specificity (discriminative scores)", py::arg("S"), py::arg("sample_assignments"));

		m.def("encode_ids", py::overload_cast<vector<string>, string>(&encode_ids), "Encrypts a set of given input ids", py::arg("ids"), py::arg("passphrase"));
		m.def("decode_ids", py::overload_cast<vector<string>, string>(&decode_ids), "Decrypts a set of given encrypted ids", py::arg("ids"), py::arg("passphrase"));		

		m.def("compute_archetype_core_centrality", py::overload_cast<arma::SpMat<npdouble> &, arma::uvec>(&compute_archetype_core_centrality), "Compute the overall connectivity of each node", py::arg("G"), py::arg("sample_assignments"));		
		m.def("compute_network_diffusion", py::overload_cast<arma::SpMat<npdouble> &, arma::SpMat<npdouble> &, int, double, int>(&compute_network_diffusion), "Computes PageRank for a selected set of nodes", py::arg("G"), py::arg("X0"), py::arg("thread_no") = 4, py::arg("alpha") = 0.85, py::arg("max_it") = 3);					
		m.def("compute_sparse_network_diffusion", py::overload_cast<arma::SpMat<npdouble> &, arma::SpMat<npdouble> &, double, double, double, int>(&compute_sparse_network_diffusion), "Computes L1-regularized PageRank for a selected set of nodes", py::arg("G"), py::arg("X0"), py::arg("alpha") = 0.85, py::arg("rho") = 1e-4, py::arg("epsilon") = 0.001, py::arg("max_iter") = 20);		
		
		m.def("assess_enrichment", py::overload_cast<arma::Mat<npdouble> &, arma::SpMat<npdouble> &, int>(&assess_enrichment), "Performs enrichment analysis", py::arg("scores"), py::arg("associations"), py::arg("thread_no") = 0);		



		// Sparse overloads
		m.def("compute_pseudo_bulk", py::overload_cast<arma::SpMat<npdouble>&, uvec>(&compute_pseudo_bulk), "Computes pseudobulk profiles", py::arg("S"), py::arg("sample_assignments"));
		m.def("compute_pseudo_bulk_per_ind", py::overload_cast<arma::SpMat<npdouble>&, uvec, uvec>(&compute_pseudo_bulk_per_ind), "Computes pseudobulk profiles (groups[k1] x individuals[k2])", py::arg("S"), py::arg("sample_assignments"), py::arg("individuals"));
		m.def("renormalize_input_matrix", py::overload_cast<arma::SpMat<npdouble>&, uvec>(&renormalize_input_matrix), "Renormalized input matrix to minimize differences in means", py::arg("S"), py::arg("sample_assignments"));
		m.def("compute_archetype_feature_specificity", py::overload_cast<arma::SpMat<npdouble> &, arma::Mat<npdouble>&>(&compute_archetype_feature_specificity), "Compute feature specificity (discriminative scores)", py::arg("S"), py::arg("H"));
		m.def("compute_cluster_feature_specificity", py::overload_cast<arma::SpMat<npdouble> &, uvec>(&compute_cluster_feature_specificity), "Compute feature specificity (discriminative scores)", py::arg("S"), py::arg("sample_assignments"));


		m.def("transform_layout", py::overload_cast<arma::SpMat<npdouble> &, Mat<npdouble>&, Mat<npdouble>&, Mat<npdouble>&, int, unsigned int, int>(&transform_layout), "Project a new data into current embedding", py::arg("W"), py::arg("coor2D"), py::arg("coor3D"), py::arg("colRGB"), py::arg("compactness_level") = 50, py::arg("n_epochs") = 500, py::arg("thread_no") = 0);


		m.def("unsigned_cluster", py::overload_cast<arma::SpMat<npdouble> &, double, uvec, int>(&unsigned_cluster), "Computes graph clustering using Leiden algorith over unsigned graphs", py::arg("A"), py::arg("resolution_parameter") = 1.0, py::arg("initial_clusters") = uvec(), py::arg("seed") = 0);
		m.def("signed_cluster", py::overload_cast<arma::SpMat<npdouble> &, double, uvec, int>(&unsigned_cluster), "Computes graph clustering using Leiden algorith over signed graphs", py::arg("A"), py::arg("resolution_parameter") = 1.0, py::arg("initial_clusters") = uvec(), py::arg("seed") = 0);







#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif        
}

