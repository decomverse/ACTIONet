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
py::dict reduce_kernel(arma::Mat<npdouble> &S, int reduced_dim = 50, int iters = 5, int seed = 0, int reduction_algorithm = 1, int SVD_algorithm = 1) {
	
	ACTIONet::ReducedKernel reduction;	
	reduction = ACTIONet::reduce_kernel(S, reduced_dim, iters, seed, reduction_algorithm, SVD_algorithm);				
	

    py::dict res;
    
	res["S_r"] = reduction.S_r;		
	res["V"] = reduction.V;
	res["lambda"] = reduction.lambda;
	res["explained_var"] = reduction.exp_var;	
		
	return res;
}
py::dict reduce_kernel_sparse(arma::SpMat<npdouble> &S, int reduced_dim = 50, int iters = 5, int seed = 0, int reduction_algorithm = 1, int SVD_algorithm = 1) {
	
	ACTIONet::ReducedKernel reduction;	
	reduction = ACTIONet::reduce_kernel(S, reduced_dim, iters, seed, reduction_algorithm, SVD_algorithm);				
	

    py::dict res;
    
	res["S_r"] = reduction.S_r;		
	res["V"] = reduction.V;
	res["lambda"] = reduction.lambda;
	res["explained_var"] = reduction.exp_var;	
		
	return res;
}




// Solves min_{X} (|| AX - B ||) s.t. simplex constraint
//
// @param A Input matrix (dense)
// @param B Input matrix (dense)
// 
// @return X Solution
arma::Mat<npdouble> run_simplex_regression(arma::Mat<npdouble>& A, arma::Mat<npdouble>& B) {	
	arma::Mat<npdouble> X = ACTIONet::run_simplex_regression(A, B);
	
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


	py::dict out;	

	// It is 0-indexed, and that makes sense for column-indices
	out["selected_columns"] = selected_columns;		
	
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

	arma::field<arma::Mat<npdouble>> AA_res = ACTIONet::run_AA(A, W0, max_it, min_delta);
	
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
py::dict prune_archetypes(vector<arma::Mat<npdouble>> &C_trace_list, vector<arma::Mat<npdouble>> & H_trace_list, double min_specificity_z_threshold = -1) {	

	field<arma::Mat<npdouble>> C_trace(C_trace_list.size());
	field<arma::Mat<npdouble>> H_trace(H_trace_list.size());	
	for(int i = 0; i < C_trace_list.size(); i++) {
		C_trace[i] = C_trace_list[i]; //aw::conv_to<arma::Mat<npdouble>>::from(C_trace_list[i]);
		H_trace[i] = H_trace_list[i]; //aw::conv_to<arma::Mat<npdouble>>::from(H_trace_list[i]);
	}
	
	ACTIONet::multilevel_archetypal_decomposition results = ACTIONet::prune_archetypes(C_trace, H_trace, min_specificity_z_threshold);
	
		
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
py::dict unify_archetypes(arma::SpMat<npdouble>& G, arma::Mat<npdouble>& S_r, arma::Mat<npdouble> &C_stacked, arma::Mat<npdouble>&H_stacked, int minPoints, int minClusterSize, double outlier_threshold) {

	ACTIONet::unification_results results = ACTIONet::unify_archetypes(G, S_r, C_stacked, H_stacked, minPoints, minClusterSize, outlier_threshold);	
		
	py::dict out_list;		
	
	// It is 0-indexed, but 1-index makes more sense for assignments
	for(int i = 0; i < (int)results.archetype_groups.n_elem; i++) results.archetype_groups[i]++;
	out_list["archetype_groups"] = results.archetype_groups;
	
	out_list["C_unified"] = results.C_unified;
	out_list["H_unified"] = results.H_unified;

	// It is 0-indexed, but 1-index makes more sense for assignments
	for(int i = 0; i < (int)results.sample_assignments.n_elem; i++) results.sample_assignments[i]++;
	out_list["sample_assignments"] = results.sample_assignments;

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
arma::Mat<npdouble> compute_pseudo_bulk(arma::Mat<npdouble> &S, uvec sample_assignments) {
	
	arma::Mat<npdouble> pb = ACTIONet::compute_pseudo_bulk_per_cluster(S, sample_assignments);
	
    return pb;
}
arma::Mat<npdouble> compute_pseudo_bulk_sparse(arma::SpMat<npdouble> &S, uvec sample_assignments) {
	
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
field<arma::Mat<npdouble>> compute_pseudo_bulk_per_ind(arma::Mat<npdouble>& S, uvec sample_assignments, uvec individuals) {

	field<arma::Mat<npdouble>> pbs_list = ACTIONet::compute_pseudo_bulk_per_ind(S, sample_assignments, individuals);
	
    return pbs_list;
}
field<arma::Mat<npdouble>> compute_pseudo_bulk_per_ind_sparse(arma::SpMat<npdouble>& S, uvec sample_assignments, uvec individuals) {

	field<arma::Mat<npdouble>> pbs_list = ACTIONet::compute_pseudo_bulk_per_ind(S, sample_assignments, individuals);
	
    return pbs_list;
}


// Renormalized input matrix to minimize differences in means
//
// @param S Input matrix (dense)
// @param sample_assignments Any primary grouping - typically based on cell type/state (it has to be in {1, ..., k1})
// 
// @return A list with the first entry being the renormalized input matrix
arma::Mat<npdouble> renormalize_input_matrix(arma::Mat<npdouble>& S, uvec sample_assignments) {
	arma::Mat<npdouble> S_norm = ACTIONet::renormalize_input_matrix(S, sample_assignments);

	return(S_norm);
}
arma::SpMat<npdouble> renormalize_input_matrix_sparse(arma::SpMat<npdouble>& S, uvec sample_assignments) {
	arma::SpMat<npdouble> S_norm = ACTIONet::renormalize_input_matrix(S, sample_assignments);

	return(S_norm);
}


// Compute feature specificity (discriminative scores)
//
// @param S Input matrix (dense)
// @param H A soft membership matrix - Typically H_unified from the unify_archetypes() function.
// 
// @return A list with the over/under-logPvals
py::dict compute_archetype_feature_specificity(arma::Mat<npdouble> &S, arma::Mat<npdouble>& H) {
    
	field<arma::Mat<npdouble>> res = ACTIONet::compute_feature_specificity(S, H);

    py::dict out_list;
	out_list["archetypes"] = res(0);
	out_list["upper_significance"] = res(1);
	out_list["lower_significance"] = res(2);


	return(out_list);
}
py::dict compute_archetype_feature_specificity_sparse(arma::SpMat<npdouble> &S, arma::Mat<npdouble>& H) {
    
	field<arma::Mat<npdouble>> res = ACTIONet::compute_feature_specificity(S, H);

    py::dict out_list;
	out_list["archetypes"] = res(0);
	out_list["upper_significance"] = res(1);
	out_list["lower_significance"] = res(2);


	return(out_list);
}


py::dict compute_cluster_feature_specificity(arma::Mat<npdouble> &S, arma::uvec sample_assignments) {
    
	field<arma::Mat<npdouble>> res = ACTIONet::compute_feature_specificity(S, sample_assignments);

    py::dict out_list;
	out_list["archetypes"] = res(0);
	out_list["upper_significance"] = res(1);
	out_list["lower_significance"] = res(2);


	return(out_list);
}
py::dict compute_cluster_feature_specificity_sparse(arma::SpMat<npdouble> &S, arma::uvec sample_assignments) {
    
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
arma::Mat<npdouble> assess_enrichment(arma::Mat<npdouble> &scores, arma::Mat<npdouble> &associations, int L) {
	
	mat logPvals = ACTIONet::assess_enrichment(scores, associations, L);

	return(logPvals);
}



PYBIND11_MODULE(_ACTIONet, m) {
    m.doc() = R"pbdoc(
        ACTIONet package
        -----------------------

        .. currentmodule:: ACTIONet

        .. autosummary::
           :toctree: _generate

    )pbdoc";
    	
		m.def("reduce_kernel", py::overload_cast<arma::Mat<npdouble> &, int, int, int, int, int>(&reduce_kernel), "Computes reduced kernel matrix for a given (single-cell) profile", py::arg("S"), py::arg("reduced_dim")=50, py::arg("iters")=5, py::arg("seed")=0, py::arg("reduction_algorithm")=1, py::arg("SVD_algorithm")=1);		

		m.def("run_simplex_regression", py::overload_cast<arma::Mat<npdouble>&, arma::Mat<npdouble>&>(&run_simplex_regression), "Solves min_{X} (|| AX - B ||) s.t. simplex constraint", py::arg("A"), py::arg("B"));
		m.def("run_SPA", py::overload_cast<arma::Mat<npdouble>&, int>(&run_SPA), "Runs Successive Projection Algorithm (SPA) to solve separable NMF", py::arg("A"), py::arg("k"));
		m.def("run_AA", py::overload_cast<arma::Mat<npdouble>&, arma::Mat<npdouble>, int, double>(&run_AA), "Runs Archetypal Analysis (AA) Algorithm", py::arg("A"), py::arg("W0"), py::arg("max_it") = 50, py::arg("min_delta") = 0.01);
		m.def("run_ACTION", py::overload_cast<arma::Mat<npdouble>&, int, int, int, int, double>(&run_ACTION), "Runs multi-level ACTION decomposition method", py::arg("S_r"), py::arg("k_min")=2, py::arg("k_max")=30, py::arg("thread_no")=4, py::arg("max_it")=50, py::arg("min_delta")=0.01);


		m.def("prune_archetypes", py::overload_cast<vector<arma::Mat<npdouble>>&, vector<arma::Mat<npdouble>> &, double>(&prune_archetypes), "Filters multi-level archetypes and concatenate filtered archetypes", py::arg("C_trace"), py::arg("H_trace"), py::arg("min_specificity_z_threshold")=-1);		
		m.def("unify_archetypes", py::overload_cast<arma::SpMat<npdouble>&, arma::Mat<npdouble>&, arma::Mat<npdouble>&, arma::Mat<npdouble>&, int, int, double>(&unify_archetypes), "Identifies and aggregates redundant archetypes into equivalent classes", py::arg("G"), py::arg("S_r"), py::arg("C_stacked"), py::arg("H_stacked"), py::arg("minPoints") = 5, py::arg("minClusterSize") = 5, py::arg("outlier_threshold") = 0.0);

	
		m.def("build_ACTIONet", py::overload_cast<arma::Mat<npdouble>&, double, int, bool>(&build_ACTIONet), "Builds an interaction network from the multi-level archetypal decompositions", py::arg("H_stacked"), py::arg("density") = 1.0, py::arg("thread_no")=8, py::arg("mutual_edges_only")=true);		
		m.def("layout_ACTIONet", py::overload_cast<arma::SpMat<npdouble>&, arma::Mat<npdouble>&, int, unsigned int, int>(&layout_ACTIONet), "Performs stochastic force-directed layout on the input graph (ACTIONet)", py::arg("G"), py::arg("S_r"), py::arg("compactness_level")= 50, py::arg("n_epochs")=100, py::arg("thread_no")= 4);

		m.def("compute_pseudo_bulk", py::overload_cast<arma::Mat<npdouble>&, uvec>(&compute_pseudo_bulk), "Computes pseudobulk profiles", py::arg("S"), py::arg("sample_assignments"));
		m.def("compute_pseudo_bulk_per_ind", py::overload_cast<arma::Mat<npdouble>&, uvec, uvec>(&compute_pseudo_bulk_per_ind), "Computes pseudobulk profiles (groups[k1] x individuals[k2])", py::arg("S"), py::arg("sample_assignments"), py::arg("individuals"));
		
		m.def("renormalize_input_matrix", py::overload_cast<arma::Mat<npdouble>&, uvec>(&renormalize_input_matrix), "Renormalized input matrix to minimize differences in means", py::arg("S"), py::arg("sample_assignments"));
		m.def("compute_archetype_feature_specificity", py::overload_cast<arma::Mat<npdouble> &, arma::Mat<npdouble>&>(&compute_archetype_feature_specificity), "Compute feature specificity (discriminative scores)", py::arg("S"), py::arg("H"));
		m.def("compute_cluster_feature_specificity", py::overload_cast<arma::Mat<npdouble> &, uvec>(&compute_cluster_feature_specificity), "Compute feature specificity (discriminative scores)", py::arg("S"), py::arg("sample_assignments"));

		m.def("encode_ids", py::overload_cast<vector<string>, string>(&encode_ids), "Encrypts a set of given input ids", py::arg("ids"), py::arg("passphrase"));
		m.def("decode_ids", py::overload_cast<vector<string>, string>(&decode_ids), "Decrypts a set of given encrypted ids", py::arg("ids"), py::arg("passphrase"));		

		m.def("compute_archetype_core_centrality", py::overload_cast<arma::SpMat<npdouble> &, arma::uvec>(&compute_archetype_core_centrality), "Compute the overall connectivity of each node", py::arg("G"), py::arg("sample_assignments"));		
		m.def("compute_network_diffusion", py::overload_cast<arma::SpMat<npdouble> &, arma::SpMat<npdouble> &, int, double, int>(&compute_network_diffusion), "Computes PageRank for a selected set of nodes", py::arg("G"), py::arg("X0"), py::arg("thread_no") = 4, py::arg("alpha") = 0.85, py::arg("max_it") = 3);					
		m.def("compute_sparse_network_diffusion", py::overload_cast<arma::SpMat<npdouble> &, arma::SpMat<npdouble> &, double, double, double, int>(&compute_sparse_network_diffusion), "Computes L1-regularized PageRank for a selected set of nodes", py::arg("G"), py::arg("X0"), py::arg("alpha") = 0.85, py::arg("rho") = 1e-4, py::arg("epsilon") = 0.001, py::arg("max_iter") = 20);		
		
		m.def("assess_enrichment", py::overload_cast<arma::Mat<npdouble> &, arma::Mat<npdouble> &, int>(&assess_enrichment), "Performs enrichment analysis", py::arg("scores"), py::arg("associations"), py::arg("L") = 1000);		



		// Sparse overloads
		m.def("reduce_kernel_sparse", py::overload_cast<arma::SpMat<npdouble> &, int, int, int, int, int>(&reduce_kernel_sparse), "Computes reduced kernel matrix for a given (single-cell) profile", py::arg("S"), py::arg("reduced_dim")=50, py::arg("iters")=5, py::arg("seed")=0, py::arg("reduction_algorithm")=1, py::arg("SVD_algorithm")=1);		
		m.def("compute_pseudo_bulk_sparse", py::overload_cast<arma::SpMat<npdouble>&, uvec>(&compute_pseudo_bulk_sparse), "Computes pseudobulk profiles", py::arg("S"), py::arg("sample_assignments"));
		m.def("compute_pseudo_bulk_per_ind_sparse", py::overload_cast<arma::SpMat<npdouble>&, uvec, uvec>(&compute_pseudo_bulk_per_ind_sparse), "Computes pseudobulk profiles (groups[k1] x individuals[k2])", py::arg("S"), py::arg("sample_assignments"), py::arg("individuals"));
		m.def("renormalize_input_matrix_sparse", py::overload_cast<arma::SpMat<npdouble>&, uvec>(&renormalize_input_matrix_sparse), "Renormalized input matrix to minimize differences in means", py::arg("S"), py::arg("sample_assignments"));
		m.def("compute_archetype_feature_specificity_sparse", py::overload_cast<arma::SpMat<npdouble> &, arma::Mat<npdouble>&>(&compute_archetype_feature_specificity_sparse), "Compute feature specificity (discriminative scores)", py::arg("S"), py::arg("H"));
		m.def("compute_cluster_feature_specificity_sparse", py::overload_cast<arma::SpMat<npdouble> &, uvec>(&compute_cluster_feature_specificity_sparse), "Compute feature specificity (discriminative scores)", py::arg("S"), py::arg("sample_assignments"));




#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif        
}

