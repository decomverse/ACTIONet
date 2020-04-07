// Computes reduced kernel matrix for a given (single-cell) profile
//
// @param S Input matrix (sparse)
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
// 
// @examples
// S = logcounts(sce)
// reduction.out = reduce(S, reduced_dim = 50)
// S_r = reduction.out$S_r
py::dict reduce_kernel(const arma::SpMat<npdouble>&S, int reduced_dim = 50, int iter = 5, int seed = 0, int reduction_algorithm = 1, int SVD_algorithm = 1) {
	
	ACTIONet::ReducedKernel reduction;	
	reduction = ACTIONet::reduce_kernel(S, reduced_dim, iter, seed, reduction_algorithm, SVD_algorithm);				
	

    py::dict res;
    
	res["S_r"] = reduction.S_r;		
	res["V"] = reduction.V;
	res["lambda"] = reduction.lambda;
	res["explained_var"] = reduction.exp_var;	
		
	return res;
}


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
// 
// @examples
// S = logcounts(sce)
// reduction.out = reduce(S, reduced_dim = 50)
// S_r = reduction.out$S_r
py::dict reduce_kernel(const dmat&S, int reduced_dim = 50, int iter = 5, int seed = 0, int reduction_algorithm = 1, int SVD_algorithm = 1) {
	
	ACTIONet::ReducedKernel reduction;	
	reduction = ACTIONet::reduce_kernel(S, reduced_dim, iter, seed, reduction_algorithm, SVD_algorithm);				
	

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
// 
// @examples
// C = ACTION.out$C[[10]]
// A = S_r * C
// B = S_r
// H = run_simplex_regression(A, B)
dmat run_simplex_regression(const dmat& A, const dmat& B) {	
	dmat X = ACTIONet::run_simplex_regression(A, B);
	
	return X;
}


// Runs Successive Projection Algorithm (SPA) to solve separable NMF
//
// @param A Input matrix (dense)
// @param k Number of columns to select
// 
// @return A named list with entries 'selected_columns' and 'norms'
// @examples
// H = run_SPA(S_r, 10)
py::dict run_SPA(const dmat& A, int k) {	

	ACTIONet::SPA_results res = ACTIONet::run_SPA(A, k);
	uvec selected_columns = res.selected_columns;
	
	vec cols = conv_to<vec>::from(selected_columns);

	py::dict out;	
	out["selected_columns"] = cols;		
	out["norms"] = res.column_norms;
		
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
// @examples
// ACTION.out = run_ACTION(S_r, k_max = 10)
// H8 = ACTION.out$H[[8]]
// cell.assignments = apply(H8, 2, which.max)
py::dict run_ACTION(const dmat& S_r, int k_min = 2, int k_max=30, int thread_no = 4) {	

	ACTIONet::ACTION_results trace = ACTIONet::run_ACTION(S_r, k_min, k_max, thread_no);

	py::dict res;
	
    py::list C(k_max+1);
	for (int i = k_min; i <= k_max; i++) {
		C[i] = trace.C[i];
	}
	res["C"] = C;	

    py::list H(k_max);
	for (int i = k_min; i <= k_max; i++) {
		H.append(trace.H[i]);
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
// @examples
// S = logcounts(sce)
// reduction.out = reduce(S, reduced_dim = 50)
// S_r = reduction.out$S_r
// ACTION.out = run_ACTION(S_r, k_max = 10)
// reconstruction.out = reconstruct_archetypes(S, ACTION.out$C, ACTION.out$H)
py::dict prune_archetypes(const py::list& C_trace, const py::list& H_trace, double min_specificity_z_threshold = -1) {	
	int n_list = H_trace.size();
	vector<const dmat&> C_trace_vec(n_list+1);
	vector<const dmat&> H_trace_vec(n_list+1);
	for (int i = 0; i < n_list; i++) {
		if(Rf_isNull(H_trace[i])) {
			continue;
		}
		C_trace_vec[i+1] = (as<const dmat&>(C_trace[i]));
		H_trace_vec[i+1] = (as<const dmat&>(H_trace[i]));
	}
		
	ACTIONet::multilevel_archetypal_decomposition results = ACTIONet::prune_archetypes(C_trace_vec, H_trace_vec, min_specificity_z_threshold);
	
		
	py::dict out_list;		

	for(int i = 0; i < results.selected_archs.n_elem; i++) results.selected_archs[i]++;
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
// @examples
// prune.out = prune_archetypes(ACTION.out$C, ACTION.out$H)
//	G = build_ACTIONet(prune.out$H_stacked)
// unification.out = unify_archetypes(G, S_r, prune.out$C_stacked, prune.out$H_stacked)
// cell.clusters = unification.out$sample_assignments
py::dict unify_archetypes(const arma::SpMat<npdouble>& &G, const dmat& &S_r, const dmat& &C_stacked, const dmat& &H_stacked) {
	ACTIONet::unification_results results = ACTIONet::unify_archetypes(G, S_r, C_stacked, H_stacked);
	
		
	py::dict out_list;		
	
	for(int i = 0; i < results.archetype_groups.n_elem; i++) results.archetype_groups[i]++;
	out_list["archetype_groups"] = results.archetype_groups;

	out_list["C_unified"] = results.C_unified;
	out_list["H_unified"] = results.H_unified;
	
	for(int i = 0; i < results.sample_assignments.n_elem; i++) results.sample_assignments[i]++;
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
// 
// @examples
// prune.out = prune_archetypes(ACTION.out$C, ACTION.out$H)
//	G = build_ACTIONet(prune.out$H_stacked)
const arma::SpMat<npdouble>& build_ACTIONet(const dmat& &H_stacked, double density = 1.0, int thread_no=8, bool mutual_edges_only = true) {
	double M = 16, ef_construction = 200, ef = 50;
	const arma::SpMat<npdouble>& G = ACTIONet::build_ACTIONet(H_stacked, density, thread_no, M, ef_construction, ef, mutual_edges_only);
	
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
// 
// @examples
//	G = build_ACTIONet(prune.out$H_stacked)
//	vis.out = layout_ACTIONet(G, S_r)
py::dict layout_ACTIONet(const arma::SpMat<npdouble>& &G, const dmat& S_r, int compactness_level= 50, unsigned int n_epochs = 100, int thread_no = 4) {

	field<const dmat&> res = ACTIONet::layout_ACTIONet(G, S_r, compactness_level, n_epochs, thread_no);
    
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
// 
// @examples
//	encoded.ids = encode_ids(colnames(sce))
vector<string> encode_ids(vector<string> ids, string pass) {
	
	vector<string> encoded_ids(ids.size());
	
    cryptor::set_key(pass);
    for(int i = 0; i < ids.size(); i++) {
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
// 
// @examples
//	ids = decode_ids(encoded.ids)
vector<string> decode_ids(vector<string> encoded_ids, string pass) {
	
	vector<string> decoded_ids(encoded_ids.size());
	
    cryptor::set_key(pass);
    for(int i = 0; i < encoded_ids.size(); i++) {
		auto dec = cryptor::decrypt(encoded_ids[i]);
		decoded_ids[i] = dec;
	}
	
    return decoded_ids;
}


// Computes pseudobulkt profiles
//
// @param S Input matrix
// @param sample_assignments Any sample clustering/annotation (it has to be in {1, ..., max_class_num})
// 
// @return S matrix aggregated within each class of sample_assignments
// 
// @examples
// prune.out = prune_archetypes(ACTION.out$C, ACTION.out$H)
//	G = build_ACTIONet(prune.out$H_stacked)
// unification.out = unify_archetypes(G, S_r, prune.out$C_stacked, prune.out$H_stacked)
// cell.clusters = unification.out$sample_assignments
//	pbs = compute_pseudo_bulk(S, cell.clusters)
const dmat& compute_pseudo_bulk(SEXP S, uvec sample_assignments) {
	
	const dmat& pb;
    if (Rf_isS4(S)) {
		const arma::SpMat<npdouble>& tmp = as<arma::const arma::SpMat<npdouble>&>(S);
		pb = ACTIONet::compute_pseudo_bulk(tmp, sample_assignments);
    } else {
		const dmat& tmp = as<arma::const dmat&>(S);
		pb = ACTIONet::compute_pseudo_bulk(tmp, sample_assignments);
    } 	
	
    return pb;
}



// Computes pseudobulk profiles (groups[k1] x individuals[k2])
//
// @param S Input matrix
// @param sample_assignments Any primary grouping - typically based on cell type/state (it has to be in {1, ..., k1})
// @param individuals Any Secondary grouping - typically corresponds to individuals (it has to be in {1, ..., k2})
// 
// @return A list of pseudobulk profile, where each entry is matrix corresponding to one cell type/state
// 
// @examples
// prune.out = prune_archetypes(ACTION.out$C, ACTION.out$H)
//	G = build_ACTIONet(prune.out$H_stacked)
// unification.out = unify_archetypes(G, S_r, prune.out$C_stacked, prune.out$H_stacked)
// cell.clusters = unification.out$sample_assignments
//	pbs.list = compute_pseudo_bulk(S, cell.clusters, sce$individuals)
field<const dmat&> compute_pseudo_bulk_per_ind(SEXP S, uvec sample_assignments, uvec individuals) {
	
	field<const dmat&> pbs_list;
    if (Rf_isS4(S)) {
		const arma::SpMat<npdouble>& tmp = as<arma::const arma::SpMat<npdouble>&>(S);
		pbs_list = ACTIONet::compute_pseudo_bulk_per_ind(tmp, sample_assignments, individuals);
    } else {
		const dmat& tmp = as<arma::const dmat&>(S);
		pbs_list = ACTIONet::compute_pseudo_bulk_per_ind(tmp, sample_assignments, individuals);
    } 		
	
    return pbs_list;
}




// Renormalized input matrix to minimize differences in means
//
// @param S Input matrix
// @param sample_assignments Any primary grouping - typically based on cell type/state (it has to be in {1, ..., k1})
// 
// @return A list with the first entry being the renormalized input matrix
// 
// @examples
// prune.out = prune_archetypes(ACTION.out$C, ACTION.out$H)
//	G = build_ACTIONet(prune.out$H_stacked)
// unification.out = unify_archetypes(G, S_r, prune.out$C_stacked, prune.out$H_stacked)
// cell.clusters = unification.out$sample_assignments
//	S.norm = renormalize_input_matrix(S, cell.clusters)
py::dict renormalize_input_matrix(SEXP S, uvec sample_assignments) {
    py::dict res;
    if (Rf_isS4(S)) {
		const arma::SpMat<npdouble>& tmp = as<arma::const arma::SpMat<npdouble>&>(S);
		const arma::SpMat<npdouble>& S_norm = ACTIONet::renormalize_input_matrix(tmp, sample_assignments);
		res["S_norm"] = S_norm;
    } else {
		const dmat& tmp = as<arma::const dmat&>(S);
		const dmat& S_norm = ACTIONet::renormalize_input_matrix(tmp, sample_assignments);
		res["S_norm"] = S_norm;
    } 		
    

	return(res);
}


// Compute feature specificity (discriminative scores)
//
// @param S Input matrix
// @param H A soft membership matrix - Typically H_unified from the unify_archetypes() function.
// 
// @return A list with the over/under-logPvals
// 
// @examples
// prune.out = prune_archetypes(ACTION.out$C, ACTION.out$H)
//	G = build_ACTIONet(prune.out$H_stacked)
// unification.out = unify_archetypes(G, S_r, prune.out$C_stacked, prune.out$H_stacked)
//	logPvals.list = renormalize_input_matrix(S, unification.out$H_unified)
py::dict compute_feature_specificity(SEXP S, const dmat& H) {
    
    field<const dmat&> res;
    
    if (Rf_isS4(S)) {
		const arma::SpMat<npdouble>& tmp = as<arma::const arma::SpMat<npdouble>&>(S);
		res = ACTIONet::compute_feature_specificity(tmp, H);
    } else {
		const dmat& tmp = as<arma::const dmat&>(S);
		res = ACTIONet::compute_feature_specificity(tmp, H);
    } 		

    py::dict out_list;
	out_list["archetypes"] = res(0);
	out_list["upper_significance"] = res(1);
	out_list["lower_significance"] = res(2);


	return(out_list);
}
