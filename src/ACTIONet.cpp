#include <RcppArmadillo.h>
#include <ACTIONet.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

#define ARMA_USE_CXX11_RNG

//' Computes SVD decomposition
//'
//' This is direct implementation of the randomized SVD algorithm for sparse matrices:
//' Xu Feng, Yuyang Xie, and Yaohang Li, "Fast Randomzied SVD for Sparse Data," in Proc. the 10th Asian Conference on Machine Learning (ACML), Beijing, China, Nov. 2018.
//' 
//' @param A Input matrix (either a "matrix" or "sparseMatrix")
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
List FengSVD(SEXP A, int dim, int iters = 5, int seed = 0) {	
	
	field<mat> SVD_out;
    if (Rf_isS4(A)) {
		sp_mat tmp = as<arma::sp_mat>(A);
		SVD_out = ACTIONet::FengSVD(tmp, dim, iters, seed);            
    } else {
		mat tmp = as<arma::mat>(A);
		SVD_out = ACTIONet::FengSVD(tmp, dim, iters, seed);            
    } 
        
	List res;
	res["U"] = SVD_out(0);	
	res["sigma"] = SVD_out(1);	
	res["V"] = SVD_out(2);	
		
	return res;
}



//' Computes SVD decomposition
//'
//' This is direct implementation of the randomized SVD algorithm:
//' XFrom: N Halko, P. G Martinsson, and J. A Tropp. Finding structure with randomness: Probabilistic algorithms for constructing approximate matrix decompositions. Siam Review, 53(2):217-288, 2011.
//' 
//' @param A Input matrix (either a "matrix" or "sparseMatrix")
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
List HalkoSVD(SEXP A, int dim, int iters = 5, int seed = 0) {	

	field<mat> SVD_out;	
    if (Rf_isS4(A)) {
		sp_mat tmp = as<arma::sp_mat>(A);
		SVD_out = ACTIONet::HalkoSVD(tmp, dim, iters, seed);            
    } else {
		mat tmp = as<arma::mat>(A);
		SVD_out = ACTIONet::HalkoSVD(tmp, dim, iters, seed);            
    } 
	
	List res;
	
	res["U"] = SVD_out(0);	
	res["sigma"] = SVD_out(1);	
	res["V"] = SVD_out(2);	
	
	return res;
}




//' Computes reduced kernel matrix for a given (single-cell) profile
//'
//' @param S Input matrix (either a "matrix" or "sparseMatrix")
//' @param reduced_dim Dimension of the reduced kernel matrix (default=50)
//' @param iters Number of SVD iterations (default=5)
//' @param seed Random seed (default=0)
//' @param reduction_algorithm Kernel reduction algorithm. Currently only ACTION method (1) is implemented (default=1)
//' @param SVD_algorithm SVD algorithm to use. Currently supported methods are Halko (1) and Feng (2) (default=1)
//' 
//' @return A named list with S_r, V, lambda, and exp_var. \itemize{
//' \item S_r: reduced kernel matrix of size reduced_dim x #samples.
//' \item V: Associated left singular-vectors (useful for reconstructing discriminative scores for features, such as genes).
//' \item lambda, exp_var: Summary statistics of the sigular-values.
//' }
//' 
//' @examples
//' S = logcounts(sce)
//' reduction.out = reduce(S, reduced_dim = 50)
//' S_r = reduction.out$S_r
// [[Rcpp::export]]
List reduce_kernel(SEXP S, int reduced_dim = 50, int iter = 5, int seed = 0, int reduction_algorithm = 1, int SVD_algorithm = 1) {
	
	ACTIONet::ReducedKernel reduction;	
    if (Rf_isS4(S)) {
		sp_mat tmp = as<arma::sp_mat>(S);
		reduction = ACTIONet::reduce_kernel(tmp, reduced_dim, iter, seed, reduction_algorithm, SVD_algorithm);				
    } else {
		mat tmp = as<arma::mat>(S);
		reduction = ACTIONet::reduce_kernel(tmp, reduced_dim, iter, seed, reduction_algorithm, SVD_algorithm);				
    } 	
			
	List res;	
	res["S_r"] = reduction.S_r;		
	res["V"] = reduction.V;
	res["lambda"] = reduction.lambda;
	res["explained_var"] = reduction.exp_var;	
		
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
mat run_simplex_regression(mat &A, mat &B) {	
	mat X = ACTIONet::run_simplex_regression(A, B);

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
List run_SPA(mat A, int k) {	

	ACTIONet::SPA_results res = ACTIONet::run_SPA(A, k);
	uvec selected_columns = res.selected_columns;
	
	vec cols(k);
	for(int i = 0; i < k; i++) {
		cols[i] = selected_columns[i] + 1;
	}
	

	List out;	
	out["selected_columns"] = cols;		
	out["norms"] = res.column_norms;
		
	return out;
}

//' Runs multi-level ACTION decomposition method
//'
//' @param S_r Reduced kernel matrix
//' @param k_min Minimum number of archetypes to consider (default=2)
//' @param k_max Maximum number of archetypes to consider, or "depth" of decomposition (default=30)
//' @param thread_no Number of parallel threads (default=4)
//' 
//' @return A named list with entries 'C' and 'H', each a list for different values of k
//' @examples
//' ACTION.out = run_ACTION(S_r, k_max = 10)
//' H8 = ACTION.out$H[[8]]
//' cell.assignments = apply(H8, 2, which.max)
// [[Rcpp::export]]
List run_ACTION(mat S_r, int k_min = 2, int k_max=30, int thread_no = 4, int max_it = 50, double min_delta = 0.01, int type = 1) {	

	ACTIONet::ACTION_results trace = ACTIONet::run_ACTION(S_r, k_min, k_max, thread_no, max_it, min_delta, type);

	List res;
	
	List C(k_max);
	for (int i = k_min; i <= k_max; i++) {
		C[i-1] = trace.C[i];
	}
	res["C"] = C;	

	List H(k_max);
	for (int i = k_min; i <= k_max; i++) {
		H[i-1] = trace.H[i];
	}
	res["H"] = H;
	
		
	return res;
}




//' Filters multi-level archetypes and concatenate filtered archetypes.
//' (Pre-ACTIONet archetype processing)
//'
//' @param C_trace,H_trace Output of ACTION
//' @param min_specificity_z_threshold Defines the stringency of pruning nonspecific archetypes. 
//' The larger the value, the more archetypes will be filtered out (default=-1)
//' 
//' @return A named list: \itemize{
//' \item selected_archs: List of final archetypes that passed the filtering/pruning step.
//' \item C_stacked,H_stacked: Horizontal/Vertical concatenation of filtered C and H matrices, respectively.
//' }
//' @examples
//' S = logcounts(sce)
//' reduction.out = reduce(S, reduced_dim = 50)
//' S_r = reduction.out$S_r
//' ACTION.out = run_ACTION(S_r, k_max = 10)
//' reconstruction.out = reconstruct_archetypes(S, ACTION.out$C, ACTION.out$H)
// [[Rcpp::export]]
List prune_archetypes(const List& C_trace, const List& H_trace, double min_specificity_z_threshold = -1) {	
	int n_list = H_trace.size();
	field<mat> C_trace_vec(n_list+1);
	field<mat> H_trace_vec(n_list+1);
	for (int i = 0; i < n_list; i++) {
		if(Rf_isNull(H_trace[i])) {
			continue;
		}
		C_trace_vec[i+1] = (as<mat>(C_trace[i]));
		H_trace_vec[i+1] = (as<mat>(H_trace[i]));
	}
		
	ACTIONet::multilevel_archetypal_decomposition results = ACTIONet::prune_archetypes(C_trace_vec, H_trace_vec, min_specificity_z_threshold);
	
		
	List out_list;		

	for(int i = 0; i < results.selected_archs.n_elem; i++) results.selected_archs[i]++;
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
//' @param C_stacked,H_stacked Output of reconstruct_archetypes()
//' 
//' @return A named list: \itemize{
//' \item archetype_groups: Equivalent classes of archetypes (non-redundant)
//' \item C_unified,H_unified: C and H matrices of unified archetypes
//' \item sample_assignments: Assignment of samples/cells to unified archetypes
//' }
//' @examples
//' prune.out = prune_archetypes(ACTION.out$C, ACTION.out$H)
//'	G = build_ACTIONet(prune.out$H_stacked)
//' unification.out = unify_archetypes(G, S_r, prune.out$C_stacked, prune.out$H_stacked)
//' cell.clusters = unification.out$sample_assignments
// [[Rcpp::export]]
List unify_archetypes(sp_mat &G, mat &S_r, mat &C_stacked, mat &H_stacked, int minPoints = 5, int minClusterSize = 5, double outlier_threshold = 0.0) {
	ACTIONet::unification_results results = ACTIONet::unify_archetypes(G, S_r, C_stacked, H_stacked, minPoints, minClusterSize, outlier_threshold);
	
		
	List out_list;		
	
	out_list["archetype_groups"] = results.archetype_groups;

	for(int i = 0; i < results.selected_archetypes.n_elem; i++) results.selected_archetypes[i]++;
	out_list["selected_archetypes"] = results.selected_archetypes;


	out_list["C_unified"] = results.C_unified;
	out_list["H_unified"] = results.H_unified;
	
	for(int i = 0; i < results.sample_assignments.n_elem; i++) results.sample_assignments[i]++;
	out_list["sample_assignments"] = results.sample_assignments;

    return out_list;
}



//' Builds an interaction network from the multi-level archetypal decompositions
//'
//' @param H_stacked Output of the prune_archetypes() function.
//' @param density Overall density of constructed graph. The higher the density, the more edges are retained (default = 1.0).
//' @param thread_no Number of parallel threads (default = 4).
//' @param mutual_edges_only Symmetrization strategy for nearest-neighbor edges. 
//' If it is true, only mutual-nearest-neighbors are returned (default=TRUE).
//' 
//' @return G Adjacency matrix of the ACTIONet graph.
//' 
//' @examples
//' prune.out = prune_archetypes(ACTION.out$C, ACTION.out$H)
//'	G = build_ACTIONet(prune.out$H_stacked)
// [[Rcpp::export]]
sp_mat build_ACTIONet(mat H_stacked, double density = 1.0, int thread_no=8, bool mutual_edges_only = true) {
	double M = 16, ef_construction = 200, ef = 50;
	sp_mat G = ACTIONet::build_ACTIONet(H_stacked, density, thread_no, M, ef_construction, ef, mutual_edges_only);
	
    return G;
}



//' Performs stochastic force-directed layout on the input graph (ACTIONet)
//'
//' @param G Adjacency matrix of the ACTIONet graph
//' @param S_r Reduced kernel matrix (is used for reproducible initialization).
//' @param compactness_level A value between 0-100, indicating the compactness of ACTIONet layout (default=50)
//' @param n_epochs Number of epochs for SGD algorithm (default=100).
//' @param thread_no Number of threads.
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
List layout_ACTIONet(sp_mat &G, mat S_r, int compactness_level= 50, unsigned int n_epochs = 500, int thread_no = 4) {

	field<mat> res = ACTIONet::layout_ACTIONet(G, S_r, compactness_level, n_epochs, thread_no);
    
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
vector<string> encode_ids(vector<string> ids, string pass) {
	
	vector<string> encoded_ids(ids.size());
	
    cryptor::set_key(pass);
    for(int i = 0; i < ids.size(); i++) {
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
vector<string> decode_ids(vector<string> encoded_ids, string pass) {
	
	vector<string> decoded_ids(encoded_ids.size());
	
    cryptor::set_key(pass);
    for(int i = 0; i < encoded_ids.size(); i++) {
		auto dec = cryptor::decrypt(encoded_ids[i]);
		decoded_ids[i] = dec;
	}
	
    return decoded_ids;
}


//' Computes pseudobulk profiles
//'
//' @param S Input matrix
//' @param sample_assignments Any sample clustering/annotation (it has to be in {1, ..., max_class_num})
//' 
//' @return S matrix aggregated within each class of sample_assignments
//' 
//' @examples
//' prune.out = prune_archetypes(ACTION.out$C, ACTION.out$H)
//'	G = build_ACTIONet(prune.out$H_stacked)
//' unification.out = unify_archetypes(G, S_r, prune.out$C_stacked, prune.out$H_stacked)
//' cell.clusters = unification.out$sample_assignments
//'	pbs = compute_pseudo_bulk(S, cell.clusters)
// [[Rcpp::export]]
mat compute_pseudo_bulk(SEXP S, arma::Col<unsigned long long> sample_assignments) {
	
	mat pb;
    if (Rf_isS4(S)) {
		sp_mat tmp = as<arma::sp_mat>(S);
		pb = ACTIONet::compute_pseudo_bulk_per_cluster(tmp, sample_assignments);
    } else {
		mat tmp = as<arma::mat>(S);
		pb = ACTIONet::compute_pseudo_bulk_per_cluster(tmp, sample_assignments);
    } 	
	
    return pb;
}



//' Computes pseudobulk profiles (groups[k1] x individuals[k2])
//'
//' @param S Input matrix
//' @param sample_assignments Any primary grouping - typically based on cell type/state (it has to be in {1, ..., k1})
//' @param individuals Any Secondary grouping - typically corresponds to individuals (it has to be in {1, ..., k2})
//' 
//' @return A list of pseudobulk profile, where each entry is matrix corresponding to one cell type/state
//' 
//' @examples
//' prune.out = prune_archetypes(ACTION.out$C, ACTION.out$H)
//'	G = build_ACTIONet(prune.out$H_stacked)
//' unification.out = unify_archetypes(G, S_r, prune.out$C_stacked, prune.out$H_stacked)
//' cell.clusters = unification.out$sample_assignments
//'	pbs.list = compute_pseudo_bulk(S, cell.clusters, sce$individuals)
// [[Rcpp::export]]
field<mat> compute_pseudo_bulk_per_ind(SEXP S, arma::Col<unsigned long long> sample_assignments, arma::Col<unsigned long long> individuals) {
	
	field<mat> pbs_list;
    if (Rf_isS4(S)) {
		sp_mat tmp = as<arma::sp_mat>(S);
		pbs_list = ACTIONet::compute_pseudo_bulk_per_ind(tmp, sample_assignments, individuals);
    } else {
		mat tmp = as<arma::mat>(S);
		pbs_list = ACTIONet::compute_pseudo_bulk_per_ind(tmp, sample_assignments, individuals);
    } 		
	
    return pbs_list;
}




//' Renormalized input matrix to minimize differences in means
//'
//' @param S Input matrix
//' @param sample_assignments Any primary grouping - typically based on cell type/state (it has to be in {1, ..., k1})
//' 
//' @return A list with the first entry being the renormalized input matrix
//' 
//' @examples
//' prune.out = prune_archetypes(ACTION.out$C, ACTION.out$H)
//'	G = build_ACTIONet(prune.out$H_stacked)
//' unification.out = unify_archetypes(G, S_r, prune.out$C_stacked, prune.out$H_stacked)
//' cell.clusters = unification.out$sample_assignments
//'	S.norm = renormalize_input_matrix(S, cell.clusters)
// [[Rcpp::export]]
List renormalize_input_matrix(SEXP S, arma::Col<unsigned long long> sample_assignments) {
    List res;
    if (Rf_isS4(S)) {
		sp_mat tmp = as<arma::sp_mat>(S);
		sp_mat S_norm = ACTIONet::renormalize_input_matrix(tmp, sample_assignments);
		res["S_norm"] = S_norm;
    } else {
		mat tmp = as<arma::mat>(S);
		mat S_norm = ACTIONet::renormalize_input_matrix(tmp, sample_assignments);
		res["S_norm"] = S_norm;
    } 		
    

	return(res);
}


//' Compute feature specificity (from archetype footprints)
//'
//' @param S Input matrix
//' @param H A soft membership matrix - Typically H_unified from the unify_archetypes() function.
//' 
//' @return A list with the over/under-logPvals
//' 
//' @examples
//' prune.out = prune_archetypes(ACTION.out$C, ACTION.out$H)
//'	G = build_ACTIONet(prune.out$H_stacked)
//' unification.out = unify_archetypes(G, S_r, prune.out$C_stacked, prune.out$H_stacked)
//' cell.clusters = unification.out$sample_assignments
//'	S.norm = renormalize_input_matrix(S, cell.clusters)
//'	logPvals.list = compute_archetype_feature_specificity(S.norm, unification.out$H_unified)
//' specificity.scores = logPvals.list$upper_significance
// [[Rcpp::export]]
List compute_archetype_feature_specificity(SEXP S, mat &H) {
    
    field<mat> res;
    
    if (Rf_isS4(S)) {
		sp_mat tmp = as<arma::sp_mat>(S);
		res = ACTIONet::compute_feature_specificity(tmp, H);
    } else {
		mat tmp = as<arma::mat>(S);
		res = ACTIONet::compute_feature_specificity(tmp, H);
    } 		

    List out_list;
	out_list["archetypes"] = res(0);
	out_list["upper_significance"] = res(1);
	out_list["lower_significance"] = res(2);


	return(out_list);
}

//' Compute feature specificity (from cluster assignments)
//'
//' @param S Input matrix
//' @param sample_assignments Vector of cluster assignments
//' 
//' @return A list with the over/under-logPvals
//' 
//' @examples
//' prune.out = prune_archetypes(ACTION.out$C, ACTION.out$H)
//'	G = build_ACTIONet(prune.out$H_stacked)
//' unification.out = unify_archetypes(G, S_r, prune.out$C_stacked, prune.out$H_stacked)
//' cell.clusters = unification.out$sample_assignments
//'	S.norm = renormalize_input_matrix(S, cell.clusters)
//'	logPvals.list = compute_cluster_feature_specificity(S.norm, cell.clusters)
//' specificity.scores = logPvals.list$upper_significance
// [[Rcpp::export]]
List compute_cluster_feature_specificity(SEXP S, uvec sample_assignments) {
    
    field<mat> res;
    
    if (Rf_isS4(S)) {
		sp_mat tmp = as<arma::sp_mat>(S);
		res = ACTIONet::compute_feature_specificity(tmp, sample_assignments);
    } else {
		mat tmp = as<arma::mat>(S);
		res = ACTIONet::compute_feature_specificity(tmp, sample_assignments);
    } 		

    List out_list;
	out_list["average_profile"] = res(0);
	out_list["upper_significance"] = res(1);
	out_list["lower_significance"] = res(2);


	return(out_list);
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
uvec compute_core_number(sp_mat &G) {

	uvec core_num = ACTIONet::compute_core_number(G);

	return(core_num);
}


//' Compute coreness of subgraph vertices induced by each archetype
//'
//' @param G Input graph
//' @param sample_assignments Archetype discretization (output of unify_archetypes())
//' 
//' @return cn core-number of each graph node
//' 
//' @examples
//' G = colNets(ace)$ACTIONet
//' assignments = ace$archetype.assignment
//' connectivity = compute_core_number(G, assignments)
// [[Rcpp::export]]
vec compute_archetype_core_centrality(sp_mat &G, uvec sample_assignments) {

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
//' 
//' @examples
//' G = colNets(ace)$ACTIONet
//' gene.expression = Matrix::t(logcounts(ace))[c("CD19", "CD14", "CD16"), ]
//' smoothed.expression = compute_network_diffusion(G, gene.expression)
// [[Rcpp::export]]
mat compute_network_diffusion(sp_mat &G, sp_mat &X0, int thread_no = 4, double alpha = 0.85, int max_it = 3) {

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
//' 
//' @examples
//' G = colNets(ace)$ACTIONet
//' gene.expression = Matrix::t(logcounts(ace))[c("CD19", "CD14", "CD16"), ]
//' smoothed.expression = compute_sparse_network_diffusion(G, gene.expression)
// [[Rcpp::export]]
sp_mat compute_sparse_network_diffusion(sp_mat &G, sp_mat &X0, double alpha = 0.85, double rho = 1e-4, double epsilon = 0.001, int max_iter = 20) {		

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
//' 
//' @examples
//' data("gProfilerDB_human")
//' G = colNets(ace)$ACTIONet
//' associations = gProfilerDB_human$SYMBOL$REAC
//' common.genes = intersect(rownames(ace), rownames(associations))
//' specificity_scores = rowFactors(ace)[["archetype_gene_specificity"]]
//' logPvals = compute_feature_specificity(specificity_scores[common.genes, ], annotations[common.genes, ])
//' rownames(logPvals) = colnames(specificity_scores)
//' colnames(logPvals) = colnames(annotations)
// [[Rcpp::export]]
mat assess_enrichment(mat &scores, mat &associations, int L) {
	
	mat logPvals = ACTIONet::assess_enrichment(scores, associations, L);

	return(logPvals);
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
vec NetDBSCAN(SEXP G, int minPts = 10, double eps = 0.5, double alpha = 0.85) {

	sp_mat Adj;
    if (Rf_isS4(G)) {
		Adj = as<arma::sp_mat>(G);
    } else {
		Adj = sp_mat(as<arma::mat>(G));
    }     
	
	vec clusters = ACTIONet::NetDBSCAN(Adj, minPts, eps, alpha);

	return(clusters);
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
List run_HDBSCAN(mat &X, int minPoints = 5, int minClusterSize = 5) {
    
    field<vec> res = ACTIONet::run_HDBSCAN(X, minPoints, minClusterSize);
    
    List out_list;
	out_list["labels"] = res(0);
	out_list["membershipProbabilities"] = res(1);
	out_list["outlierScores"] = res(2);


	return(out_list);
}

