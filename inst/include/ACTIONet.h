#ifndef ACTIONet_H
#define ACTIONet_H

#include <stdio.h>
#include <vector>
#include <string>
#include <cmath>
#include <sys/types.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <errno.h>
#include <sys/stat.h>
#include <getopt.h>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <unordered_map>
#include <ctime>
#include <thread>
#include <atomic>


#include <arma_base.h>
#include <my_utils.h>
#include <gradient.h>
#include <sampler.h>
#include <tauprng.h>
#include <colorspace.h>
#include <cryptor.hpp>
#include <hdbscan.hpp>

// SVD algorithms
#define HALKO_ALG 1
#define FENG_ALG 2

// Kernel reduction algorithms
#define ACTIONRED_ALG 1


// Symmetrization methods for ACTIONet network edges
#define ACTIONet_AND_SYM 1
#define ACTIONet_OR_SYM 2



    
namespace ACTIONet {
	// Main structures	
		// To store the output of reduce_kernel()
		struct ReducedKernel {
			mat S_r;
			mat V;
			vec lambda;
			vec exp_var;
		};

		// To store the output of run_SPA()
		struct SPA_results {
			uvec selected_columns;
			vec column_norms;
		};
		
		// To store the output of run_ACTION()
		struct ACTION_results {
			field<uvec> selected_cols;
			field<mat> H;
			field<mat> C;
		};

		// To store the output of reconstruct_archetypes()
		struct multilevel_archetypal_decomposition {
			uvec selected_archs; // If hub removal requested, this will hold the indices of retained archetypes
			mat C_stacked; // Stacking of C matrices, after potentially removing the hub archetypes
			mat H_stacked; // Stacking of H matrices, after potentially removing the hub archetypes
		};			

		// To store the output of unify_archetypes()
		struct unification_results {
			vec archetype_groups; 
			uvec selected_archetypes;
			mat C_unified;
			mat H_unified;
			uvec sample_assignments;
		};	
		
	// Low-level functions
	// *********************************
		// Basic (randomized) SVD algorithms
		field<mat> FengSVD(sp_mat &A, int dim, int iters, int seed);
		field<mat> FengSVD(mat &A, int dim, int iters, int seed);
		
		field<mat> HalkoSVD(mat &A, int dim, int max_it, int seed);		
		field<mat> HalkoSVD(sp_mat &A, int dim, int max_it, int seed);	

		// ACTION-based kernel reduction algorithms
		ReducedKernel ACTION_reduction(sp_mat &A, int dim, int iter, int seed, int SVD_algorithm);
		ReducedKernel ACTION_reduction(mat &A, int dim, int iter, int seed, int SVD_algorithm);
		
		// Successive Projection Algorithm (SPA) to solve separable NMF
		SPA_results run_SPA(mat M, int k);
		
		// min_{X} (|| AX - B ||) s.t. simplex constraint using ACTIVE Set Method
		mat run_simplex_regression(mat &A, mat &B);
		mat run_simplex_regression_proxdist(mat &X, mat &Y, int pmaxiter, int pincmaxiter);


		// Robust archetypal analysis method
		field<mat> run_AA (mat &S, mat &W0, int max_it, double min_delta);	
		
	// *********************************
		
		
	// Entry-points to compute a reduced kernel matrix	
		ReducedKernel reduce_kernel(sp_mat &S, int dim, int iter, int seed, int reduction_algorithm, int SVD_algorithm);
		ReducedKernel reduce_kernel(mat &S, int dim, int iter, int seed, int reduction_algorithm, int SVD_algorithm);

	// ACTION decomposition
		ACTION_results run_ACTION(mat &S_r, int k_min, int k_max, int thread_no, int max_it, double min_delta);
		ACTION_results run_ACTION_dev(mat &S_r, int k_min, int k_max, int thread_no, bool auto_stop, int max_it, double min_delta);

		ACTION_results run_weighted_ACTION(mat &S_r, vec w, int k_min, int k_max, int thread_no, int max_it, double min_delta);

	// Pre-ACTIONet archetype filtering/aggregation
	// To prune archetypes across different levels and concatenate the resulting archetypes
		multilevel_archetypal_decomposition prune_archetypes(field<mat> C_trace, field<mat> H_trace, double min_specificity_z_threshold);
		
		
	// Post-ACTIONet archetype filtering/aggregation
	// To unify redundant archetypes across different levels
		unification_results unify_archetypes(sp_mat &G, mat &S_r, mat &C_stacked, mat &H_stacked, int minPoints, int minClusterSize, double outlier_threshold);
	
	
	// Main functions to build an interaction network from multi-level archetypal decompositions
		sp_mat build_ACTIONet_JS_KstarNN(mat H_stacked, double density, int thread_no, double M, double ef_construction, double ef, bool mutual_edges_only);
		sp_mat build_ACTIONet_JS_KstarNN_v2(mat H_stacked, double density, int thread_no, double M, double ef_construction, double ef, bool mutual_edges_only);
		sp_mat build_ACTIONet_JS_KNN(mat H_stacked, int k, int thread_no, double M, double ef_construction, double ef, bool mutual_edges_only);
		
		sp_mat build_ACTIONet(mat H_stacked, double density, int thread_no, double M, double ef_construction, double ef, bool mutual_edges_only);


	// SGD-based force-directed layout (adopted and modified from the UMAP implementation)
		field<mat> layout_ACTIONet(sp_mat &G, mat &S_r, int compactness_level, unsigned int n_epochs, int thread_no);	
	
	
	// Methods for pseudo-bulk construction
		mat compute_pseudo_bulk_per_archetype(sp_mat &S, mat& H);
		mat compute_pseudo_bulk_per_archetype(mat &S, mat& H);
		mat compute_pseudo_bulk_per_cluster(sp_mat &S, arma::Col<unsigned long long> sample_assignments);
		mat compute_pseudo_bulk_per_cluster(mat &S, arma::Col<unsigned long long> sample_assignments);
		field<mat> compute_pseudo_bulk_per_ind(sp_mat &S, arma::Col<unsigned long long> sample_assignments, arma::Col<unsigned long long> individuals);
		field<mat> compute_pseudo_bulk_per_ind(mat &S, arma::Col<unsigned long long> sample_assignments, arma::Col<unsigned long long> individuals);
		
		
	// Methods for renormalizing input matrix within and between each class
		mat renormalize_input_matrix(mat &S, arma::Col<unsigned long long> sample_assignments);
		sp_mat renormalize_input_matrix(sp_mat &S, arma::Col<unsigned long long> sample_assignments);
		
		
	// Methods for computing feature specificity/discriminative-scores
		field<mat> compute_feature_specificity(sp_mat &S, mat &H);
		field<mat> compute_feature_specificity(mat &S, mat &H);
		field<mat> compute_feature_specificity(sp_mat &S, uvec sample_assignments);
		field<mat> compute_feature_specificity(mat &S, uvec sample_assignments);

	
	// Methods for feature enrichment analysis
		mat assess_enrichment(mat &scores, mat &associations, int L);


	// Network tools
		uvec compute_core_number(sp_mat &G);
		vec compute_archetype_core_centrality(sp_mat &G, uvec sample_assignments);		
		mat compute_network_diffusion(sp_mat &G, sp_mat &X0, int thread_no, double alpha, int max_it);
		sp_mat compute_sparse_network_diffusion(sp_mat &G, sp_mat &X0, double alpha, double rho, double epsilon, int max_iter);
		vec NetDBSCAN(sp_mat& G, int minPts, double eps, double alpha_val);
	
		field<vec> run_HDBSCAN(mat &X, int minPoints, int minClusterSize);

}

#endif
