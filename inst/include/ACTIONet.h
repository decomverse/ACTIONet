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


// Visualization associated parameter settings
#define TUMAP_LAYOUT 0
#define UMAP_LAYOUT 1
#define GRAPHVIS_LAYOUT 2

const double UMAP_A[101] = {1.93280839781719,1.89560586588002,1.85873666431227,1.82221007490834,1.78603612060048,1.75022496320214,1.71478579945151,1.67972997626197,1.64506544270902,1.610800661285,1.57694346052399,1.54350101780511,1.51047986323257,1.47788588612333,1.44572435168023,1.41399925414561,1.38271638006498,1.35187804260518,1.3214872860387,1.29154663185922,1.26205810311418,1.23302325071067,1.20444317424075,1.17631854866857,1.14864964274379,1.12143634262879,1.09467817152021,1.0683743100033,1.04252361298475,1.01712481754341,0.992175611624647,0.967674513244996,0.943619207179927,0.920007077834315,0.896835219021839,0.874100443595699,0.851800999392949,0.829931994792615,0.808490430178554,0.787472613514984,0.766873638278737,0.746690990400437,0.726919886947928,0.707556026044195,0.688594985599233,0.670032232635194,0.651864066568649,0.634084192553475,0.616688494561969,0.599672088669339,0.583030020204371,0.5667572718654,0.550848768322639,0.535299383967892,0.520103947257001,0.505257246260431,0.490754031684977,0.476589022213249,0.46275690208242,0.449252325341552,0.436069912245555,0.423205974605747,0.4106531652521,0.39840668039948,0.386461380891047,0.374811984314975,0.363453224264704,0.352379851902848,0.341586644916259,0.331068403184832,0.320819956874279,0.31083616902857,0.301110995958752,0.291641183389757,0.282420831386121,0.273444955588216,0.264708614833586,0.256206914916444,0.247935008593902,0.239888099677924,0.232061441819675,0.224450342118235,0.217050162160312,0.209856317524031,0.202864281204524,0.196069583611474,0.189467814398248,0.183054621446351,0.176825713015038,0.17077685928726,0.164903890637922,0.159202699934773,0.153669242222215,0.148299535941784,0.143089661250278,0.138035764053223,0.133134049958711,0.12838079222654,0.123772324007265,0.119305671122251,0.114976081494676};
const double UMAP_B[101] = {0.790494973419029,0.80063784415826,0.810876441425738,0.821199202674006,0.831595366275022,0.84205539236769,0.852571713401325,0.863135518043442,0.873741680140683,0.884384956993888,0.895060878257082,0.905765637284042,0.916495998501859,0.927249214280422,0.938022954467018,0.948815759038301,0.95962499558526,0.970449732070657,0.981288783823989,0.992141168965973,1.00300608092206,1.01388286515112,1.02477099750548,1.03567006898871,1.04657977025277,1.05749987674998,1.06843023939592,1.07937077470387,1.09032145585694,1.10128169075827,1.11225322117536,1.12323470900213,1.13422639755358,1.14522861434516,1.15624176559097,1.16726633179917,1.17830241385901,1.18934945144456,1.20040819996369,1.21147891097075,1.22256381651844,1.23366041866219,1.24477022428392,1.2558936051142,1.26703094885274,1.27818265467871,1.28934756395537,1.30052872175886,1.31172539107843,1.32293800168803,1.3341669930459,1.34541281413396,1.35667592718974,1.36795680610473,1.37925594017143,1.39057383474783,1.40191101858967,1.41326804557094,1.42464550789942,1.436044048272,1.44746436980037,1.45890393087319,1.47036701291879,1.48185337703821,1.49336326709497,1.50489726618312,1.51645596605121,1.52803997486173,1.53964990048402,1.55128637349183,1.56295003156298,1.57464152150044,1.58636409305622,1.59811350189048,1.60989278253114,1.62170263415549,1.63354377154668,1.64541692037945,1.65732282325244,1.66926223230814,1.68123591907029,1.69324466615879,1.70528927262371,1.71737055545595,1.72948934595558,1.74164649289645,1.75384285823827,1.76607932576738,1.77835679827623,1.79067619009556,1.80303844043406,1.81544450541945,1.82789536263139,1.84039200538657,1.85293545544251,1.86552674229068,1.87816693701183,1.89085711093115,1.90359837758981,1.91638829237987,1.92923479503841};

#define NEGATIVE_SAMPLE_RATE 5.0
#define LEARNING_RATE 1.0	
#define UMAP_SEED 0	
#define GAMMA 1.0 

vector<double> optimize_layout(
    const apumap_gradient& gradient,
    vector<double>& head_embedding,
    vector<double>& tail_embedding,
    const vector<unsigned int>& positive_head,
    const vector<unsigned int>& positive_tail,
    unsigned int n_epochs, 
    unsigned int n_vertices,
    const vector<double>& epochs_per_sample,
    double initial_alpha,
    double negative_sample_rate,
    unsigned int seed);
    
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
		ACTION_results run_ACTION(mat S_r, int k_min, int k_max, int thread_no, int max_it, double min_delta, int type);
		ACTION_results run_ACTION_dev(mat S_r, int k_min, int k_max, int thread_no, bool auto_stop, int max_it, double min_delta);

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
