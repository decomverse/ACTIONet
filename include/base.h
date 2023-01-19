#ifndef BASE_H
#define BASE_H

#include <errno.h>
#include <getopt.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <ctime>
#include <map>
#include <set>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include <gradient.h>
#include <sampler.h>
#include <tauprng.h>
#include <mini_thread/mini_thread.h>
#include <colorspace.h>
#include <my_utils.h>
#include <math/aarand/aarand.hpp>
#include <hdbscan.hpp>
#include <pcg_random.hpp>

using namespace mini_thread;

#define STATS_GO_INLINE
#define STATS_ENABLE_ARMA_WRAPPERS
#include <stats.hpp>

// SVD algorithms
#define FULL_SVD -1
#define IRLB_ALG 0
#define HALKO_ALG 1
#define FENG_ALG 2

// Kernel reduction algorithms
#define ACTIONRED_ALG 1

// Symmetrization methods for ACTIONet network edges
#define ACTIONet_AND_SYM 1
#define ACTIONet_OR_SYM 2

#define SYS_THREADS_DEF (std::thread::hardware_concurrency() - 2)

// Visualization associated parameter settings
#define TUMAP_LAYOUT 0
#define UMAP_LAYOUT 1
#define GRAPHVIS_LAYOUT 2

struct mvtrace_obj
{
  vector<uvec> selected_cols;

  vector<mat> H_primary;
  vector<mat> C_primary;

  vector<mat> H_secondary;
  vector<mat> C_secondary;

  vector<mat> C_consensus;
};

struct full_trace
{
  vector<mvtrace_obj> indiv_trace;
  vector<mat> H_consensus;
};

// s_gd2 visualization
void layout_unweighted(int n, double *X, int m, int *I, int *J, int t_max,
                       double eps, int seed);
void layout_weighted(int n, double *X, int m, int *I, int *J, double *V,
                     int t_max, double eps, int seed);
void layout_unweighted_convergent(int n, double *X, int m, int *I, int *J,
                                  int t_max, double eps, double delta,
                                  int t_maxmax, int seed);
void layout_weighted_convergent(int n, double *X, int m, int *I, int *J,
                                double *V, int t_max, double eps, double delta,
                                int t_maxmax, int seed);

void layout_sparse_unweighted(int n, double *X, int m, int *I, int *J, int p,
                              int t_max, double eps, int seed);
void layout_sparse_weighted(int n, double *X, int m, int *I, int *J, double *V,
                            int p, int t_max, double eps, int seed);

void mds_direct(int n, int kd, double *X, double *d, double *w, int t_max,
                double *etas, int seed);

namespace ACTIONet
{
  struct PCHAkernel_ret
  {
    arma::mat S;
    arma::mat C;
    float SSE;
    arma::mat alphaC;
  };

  // Main structures
  // To store the output of compute_AA_coreset()
  struct Coreset
  {
    mat S_coreset;
    vec w_coreset;
    uvec index;
  };

  // To store the output of run_SPA()
  struct SPA_results
  {
    uvec selected_columns;
    vec column_norms;
  };

  // To store the output of run_ACTION()
  struct ACTION_results
  {
    field<uvec> selected_cols;
    field<mat> H;
    field<mat> C;
  };

  // To store the output of run_ACTION()
  struct Online_ACTION_results
  {
    field<uvec> selected_cols;
    field<mat> A;
    field<mat> B;
    field<mat> C;
    field<mat> D;
  };

  // To store the output of reconstruct_archetypes()
  struct multilevel_archetypal_decomposition
  {
    uvec selected_archs; // If hub removal requested, this will hold the indices
                         // of retained archetypes
    mat C_stacked;       // Stacking of C matrices, after potentially removing the hub
                         // archetypes
    mat H_stacked;       // Stacking of H matrices, after potentially removing the hub
                         // archetypes
  };

  // To store the output of unify_archetypes()
  struct unification_results
  {
    mat dag_adj;
    vec dag_node_annotations;
    uvec selected_archetypes;
    mat C_unified;
    mat H_unified;
    uvec assigned_archetypes;
    vec archetype_group;
    mat arch_membership_weights;
  };

  // Low-level functions
  // *********************************
  // Basic (randomized) SVD algorithms
  field<mat> FengSVD(sp_mat &A, int dim, int iters, int seed, int verbose);
  field<mat> FengSVD(mat &A, int dim, int iters, int seed, int verbose);

  field<mat> HalkoSVD(mat &A, int dim, int max_it, int seed, int verbose);
  field<mat> HalkoSVD(sp_mat &A, int dim, int max_it, int seed, int verbose);

  field<mat> IRLB_SVD(mat &A, int dim, int max_it, int seed, int verbose);
  // field<mat> IRLB_SVD(sp_mat &A, int dim, int max_it, int seed);
  field<mat> IRLB_SVD(sp_mat &A, int dim, int iters, int seed, int verbose);

  // Successive Projection Algorithm (SPA) to solve separable NMF
  SPA_results run_SPA(mat &M, int k);
  SPA_results run_SPA_rows_sparse(sp_mat &A, int k);

  // min_{X} (|| AX - B ||) s.t. simplex constraint using ACTIVE Set Method
  // mat run_simplex_regression(mat &A, mat &B);
  mat run_simplex_regression(mat &A, mat &B, bool computeXtX);
  mat run_simplex_regression_proxdist(mat &X, mat &Y, int pmaxiter,
                                      int pincmaxiter);

  // Robust archetypal analysis method
  field<mat> run_AA(mat &A, mat &W0, int max_it = 100, double min_delta = 1e-6);

  // Online archetypal analysis method (Online Dictionary Learning for Approximate
  // Archetypal Analysis)
  field<mat> Online_update_AA(mat &Xt, mat &D, mat &A, mat &B);
  field<mat> run_online_AA(mat &X, mat &D0, field<uvec> samples);
  Online_ACTION_results run_online_ACTION(mat &S_r, field<uvec> samples,
                                          int k_min, int k_max, int thread_no);

  // *********************************

  // Entry-points to compute a reduced kernel matrix
  field<mat> perturbedSVD(field<mat> SVD_results, mat &A, mat &B);

  field<mat> SVD2ACTIONred(sp_mat &S, field<mat> SVD_results);
  field<mat> SVD2ACTIONred(mat &S, field<mat> SVD_results);

  field<mat> PCA2ACTIONred(sp_mat &S, field<mat> PCA_results);
  field<mat> PCA2ACTIONred(mat &S, field<mat> PCA_results);

  field<mat> SVD2PCA(sp_mat &S, field<mat> SVD_results);
  field<mat> SVD2PCA(mat &S, field<mat> SVD_results);
  field<mat> PCA2SVD(sp_mat &S, field<mat> PCA_results);
  field<mat> PCA2SVD(mat &S, field<mat> PCA_results);

  field<mat> reduce_kernel(sp_mat &S, int dim, int iter, int seed,
                           int SVD_algorithm, bool prenormalize, int verbose);
  field<mat> reduce_kernel(mat &S, int dim, int iter, int seed, int SVD_algorithm,
                           bool prenormalize, int verbose);

  field<mat> ACTIONred2SVD(field<mat> SVD_results);

  field<mat> deflate_reduction(field<mat> SVD_results, mat &A, mat &B);

  field<mat> orthogonalize_batch_effect(sp_mat &S, field<mat> SVD_results,
                                        mat &design);
  field<mat> orthogonalize_batch_effect(mat &S, field<mat> SVD_results,
                                        mat &design);

  field<mat> orthogonalize_basal(sp_mat &S, field<mat> SVD_results, mat &basal);
  field<mat> orthogonalize_basal(mat &S, field<mat> SVD_results, mat &basal);

  // ACTION decomposition
  ACTION_results run_ACTION(mat &S_r, int k_min, int k_max, int thread_no,
                            int max_it, double min_delta, int normalization);
  ACTION_results run_subACTION(mat &S_r, mat &W_parent, mat &H_parent, int kk,
                               int k_min, int k_max, int thread_no, int max_it,
                               double min_delta);

  ACTION_results run_weighted_ACTION(mat &S_r, vec w, int k_min, int k_max,
                                     int thread_no, int max_it, double min_delta);

  ACTION_results run_ACTION_plus(mat &S_r, int k_min, int k_max, int max_it,
                                 double min_delta, int max_trial);

  // Pre-ACTIONet archetype filtering/aggregation
  // To prune archetypes across different levels and concatenate the resulting
  // archetypes
  multilevel_archetypal_decomposition prune_archetypes(
      field<mat> C_trace, field<mat> H_trace, double min_specificity_z_threshold,
      int min_cells);

  // Post-ACTIONet archetype filtering/aggregation
  // To unify redundant archetypes across different levels
  // unification_results unify_archetypes(sp_mat &G, mat &S_r, mat &archetypes,
  // mat &C_stacked, mat &H_stacked, int minPoints, int minClusterSize, double
  // outlier_threshold, int reduced_dim);
  unification_results unify_archetypes(mat &S_r, mat &C_stacked, mat &H_stacked,
                                       double backbone_density = 0.5,
                                       double resolution = 1.0,
                                       int min_cluster_size = 3,
                                       int thread_no = 0, int normalization = 0);

  // Main functions to build an interaction network from multi-level archetypal
  // decompositions
  sp_mat buildNetwork_KstarNN(mat H_stacked, double density, int thread_no,
                              double M, double ef_construction, double ef,
                              bool mutual_edges_only, string distance_metric);
  sp_mat buildNetwork_KstarNN_v2(mat H_stacked, double density, int thread_no,
                                 double M, double ef_construction, double ef,
                                 bool mutual_edges_only, string distance_metric);
  sp_mat buildNetwork_KNN(mat H_stacked, int k, int thread_no, double M,
                          double ef_construction, double ef,
                          bool mutual_edges_only, string distance_metric);

  sp_mat buildNetwork(mat H, string algorithm = "k*nn",
                      string distance_metric = "jsd", double density = 1.0,
                      int thread_no = 0, double M = 16,
                      double ef_construction = 200, double ef = 200,
                      bool mutual_edges_only = true, int k = 10);

  mat computeFullSim(mat &H, int thread_no);

  // SGD-based force-directed layout (adopted and modified from the UMAP (uwot)
  field<mat> layoutNetwork_xmap(sp_mat &G, mat &initial_position,
                                bool presmooth_network = false,
                                const std::string &method = "umap",
                                double min_dist = 1, double spread = 1,
                                double gamma = 1.0, unsigned int n_epochs = 500,
                                int thread_no = 0, int seed = 0,
                                double learning_rate = 1.0, int sim2dist = 2);

  mat transform_layout(sp_mat &G, mat &reference_layout, bool presmooth_network,
                       const std::string &method, double min_dist, double spread,
                       double gamma, unsigned int n_epochs, int thread_no,
                       int seed, double learning_rate, int sim2dist);

  // Methods for pseudo-bulk construction
  mat compute_pseudo_bulk_per_archetype(sp_mat &S, mat &H);
  mat compute_pseudo_bulk_per_archetype(mat &S, mat &H);
  field<mat> compute_pseudo_bulk_per_archetype_and_ind(
      sp_mat &S, mat &H, arma::Col<unsigned long long> sample_assignments);
  field<mat> compute_pseudo_bulk_per_archetype_and_ind(
      mat &S, mat &H, arma::Col<unsigned long long> sample_assignments);

  mat compute_grouped_rowsums(sp_mat &S,
                              arma::Col<unsigned long long> sample_assignments);
  mat compute_grouped_rowsums(mat &S,
                              arma::Col<unsigned long long> sample_assignments);

  mat compute_grouped_rowmeans(sp_mat &S,
                               arma::Col<unsigned long long> sample_assignments);
  mat compute_grouped_rowmeans(mat &S,
                               arma::Col<unsigned long long> sample_assignments);

  mat compute_grouped_rowvars(sp_mat &S,
                              arma::Col<unsigned long long> sample_assignments);
  mat compute_grouped_rowvars(mat &S,
                              arma::Col<unsigned long long> sample_assignments);

  // Methods for renormalizing input matrix within and between each class
  mat renormalize_input_matrix(mat &S,
                               arma::Col<unsigned long long> sample_assignments);
  sp_mat renormalize_input_matrix(
      sp_mat &S, arma::Col<unsigned long long> sample_assignments);

  // Methods for computing feature specificity/discriminative-scores
  field<mat> compute_feature_specificity_bin(sp_mat &Sb, mat &H, int thread_no);
  field<mat> compute_feature_specificity(sp_mat &S, mat &H, int thread_no);
  field<mat> compute_feature_specificity(mat &S, mat &H, int thread_no);
  field<mat> compute_feature_specificity(sp_mat &S, uvec sample_assignments,
                                         int thread_no);
  field<mat> compute_feature_specificity(mat &S, uvec sample_assignments,
                                         int thread_no);

  // Methods for feature enrichment analysis
  field<mat> assess_enrichment(mat &scores, sp_mat &associations, int thread_no);

  // Network tools
  uvec compute_core_number(sp_mat &G);
  vec compute_archetype_core_centrality(sp_mat &G, uvec sample_assignments);
  mat compute_network_diffusion(sp_mat &G, sp_mat &X0, int thread_no,
                                double alpha, int max_it);
  mat compute_network_diffusion_fast(sp_mat &G, sp_mat &X0, int thread_no,
                                     double alpha, int max_it);
  mat compute_network_diffusion_direct(sp_mat &G, sp_mat &X0, int thread_no,
                                       double alpha);
  mat compute_network_diffusion_SFMULT(sp_mat &G, sp_mat &X0, double alpha,
                                       int max_it);

  vec NetDBSCAN(sp_mat &G, int minPts, double eps, double alpha_val);

  field<vec> run_HDBSCAN(mat &X, int minPoints, int minClusterSize);

  mat MWM_hungarian(mat &G);
  umat MWM_rank1(vec u, vec v, double u_threshold, double v_threshold);

  mat Prune_PageRank(mat &U, double density);

  vec unsigned_cluster(sp_mat A, double resolution_parameter,
                       uvec initial_clusters, int seed);
  vec signed_cluster(sp_mat A, double resolution_parameter, uvec initial_clusters,
                     int seed);

  Coreset compute_AA_coreset(sp_mat &S, int m);

  mat NetEnh(mat Adj);

  mat unsigned_cluster_batch(sp_mat A, vec resolutions, uvec initial_clusters,
                             int seed);

  vec LPA(sp_mat &G, vec labels, double lambda = 0, int iters = 3,
          double sig_threshold = 3, uvec fixed_labels = uvec(),
          int thread_no = 0);

  mat compute_marker_aggregate_stats(sp_mat &G, sp_mat &S, sp_mat &marker_mat,
                                     double alpha, int max_it, int thread_no,
                                     bool ignore_baseline_expression);

  field<mat> run_AA_with_batch_correction(mat &Z, mat &W0, vec batch, int max_it,
                                          int max_correction_rounds,
                                          double lambda, double min_delta);

  ACTION_results run_ACTION_with_batch_correction(
      mat &S_r, vec batch, int k_min, int k_max, int thread_no, int max_it,
      int max_correction_rounds, double lambda, double min_delta);

  mat compute_marker_aggregate_stats_basic_sum(sp_mat &S, sp_mat &marker_mat);
  mat compute_marker_aggregate_stats_basic_sum_perm(sp_mat &S, sp_mat &marker_mat,
                                                    int perm_no, int thread_no);
  mat compute_marker_aggregate_stats_basic_sum_perm_smoothed(
      sp_mat &G, sp_mat &S, sp_mat &marker_mat, double alpha, int max_it,
      int perm_no, int thread_no);
  mat compute_marker_aggregate_stats_basic_sum_smoothed(sp_mat &G, sp_mat &S,
                                                        sp_mat &marker_mat,
                                                        double alpha, int max_it,
                                                        int perm_no,
                                                        int thread_no);
  mat compute_marker_aggregate_stats_basic_sum_smoothed_normalized(
      sp_mat &G, sp_mat &S, sp_mat &marker_mat, double, int max_it, int perm_no,
      int thread_no);
  mat compute_marker_aggregate_stats_basic_sum_perm_smoothed_v2(
      sp_mat &G, sp_mat &S, sp_mat &marker_mat, double alpha, int max_it,
      int perm_no, int thread_no);

  mat compute_marker_aggregate_stats_TFIDF_sum_smoothed(
      sp_mat &G, sp_mat &S, sp_mat &marker_mat, double alpha, int max_it,
      int perm_no, int thread_no, int normalization);

  sp_mat LSI(sp_mat &X, double size_factor);

  field<vec> autocorrelation_Moran(mat G, mat scores,
                                   int normalization_method = 1, int perm_no = 30,
                                   int thread_no = 0);
  field<vec> autocorrelation_Moran(sp_mat G, mat scores,
                                   int normalization_method = 1, int perm_no = 30,
                                   int thread_no = 0);

  field<vec> autocorrelation_Geary(mat G, mat scores,
                                   int normalization_method = 1, int perm_no = 30,
                                   int thread_no = 0);
  field<vec> autocorrelation_Geary(sp_mat G, mat scores,
                                   int normalization_method = 1, int perm_no = 30,
                                   int thread_no = 0);

  sp_mat normalize_adj(sp_mat &G, int norm_type = 1);
  mat compute_network_diffusion_Chebyshev(sp_mat &P, mat &X, int thread_no = 0,
                                          double alpha = 0.85, int max_it = 5,
                                          double res_threshold = 1e-8);
  mat compute_marker_aggregate_stats_nonparametric(mat &S, sp_mat &marker_mat,
                                                   int thread_no = 0);

  full_trace runACTION_muV(vector<mat> cell_signatures, int k_min, int k_max,
                           vec alpha, double lambda = 1, int AA_iters = 50,
                           int Opt_iters = 0, int thread_no = 0);
  mat compute_markers_eigengene(mat &S, sp_mat &marker_mat, int normalization = 0,
                                int thread_no = 0);

  vec sweepcut(sp_mat &A, vec s, int min_size = 5, int max_size = -1);

  mat normalize_scores(mat scores, int method = 1, int thread_no = 0);

  mat aggregate_genesets(sp_mat &G, sp_mat &S, sp_mat &marker_mat,
                         int network_normalization_method = 0,
                         int expression_normalization_method = 0,
                         int gene_scaling_method = 0, double post_alpha = 0.85,
                         int thread_no = 0);

  mat aggregate_genesets_mahalanobis_2archs(
      sp_mat &G, sp_mat &S, sp_mat &marker_mat,
      int network_normalization_method = 0,
      int expression_normalization_method = 0, int gene_scaling_method = 3,
      double pre_alpha = 0.15, double post_alpha = 0.85, int thread_no = 0);
  mat aggregate_genesets_mahalanobis_2gmm(
      sp_mat &G, sp_mat &S, sp_mat &marker_mat,
      int network_normalization_method = 0,
      int expression_normalization_method = 0, int gene_scaling_method = 3,
      double pre_alpha = 0.15, double post_alpha = 0.85, int thread_no = 0);
  mat aggregate_genesets_weighted_enrichment(
      sp_mat &G, sp_mat &S, sp_mat &marker_mat,
      int network_normalization_method = 0,
      int expression_normalization_method = 0, int gene_scaling_method = 3,
      double pre_alpha = 0.15, double post_alpha = 0.85, int thread_no = 0);
  mat aggregate_genesets_weighted_enrichment_permutation(
      sp_mat &G, sp_mat &S, sp_mat &marker_mat,
      int network_normalization_method = 0,
      int expression_normalization_method = 0, int gene_scaling_method = 3,
      double pre_alpha = 0.15, double post_alpha = 0.85, int thread_no = 0,
      int perm_no = 30);

  mat run_simplex_regression_FW(mat &A, mat &B, int max_iter = -1,
                                double min_diff = 0.01);
  field<mat> recursiveNMU_mine(mat M, int dim, int max_SVD_iter,
                               int max_iter_inner);
  field<mat> recursiveNMU(mat M, int dim, int max_SVD_iter, int max_iter_inner);

  mat normalize_mat(mat &X, int normalization = 0, int dim = 0);
  sp_mat normalize_mat(sp_mat &X, int normalization = 0, int dim = 0);
  vec rank_vec(vec x, int method = 0);
  vec xicor(vec xvec, vec yvec, bool compute_pval = true, int seed = 0);
  field<mat> XICOR(mat &X, mat &Y, bool compute_pval = true, int seed = 0, int thread_no = 0);

  sp_mat buildNetwork_bipartite(mat H1, mat H2, double density = 1.0,
                                int thread_no = 0, double M = 16,
                                double ef_construction = 200,
                                double ef = 200,
                                string distance_metric = "jsd");

  //  mat aggregate_genesets_weighted_enrichment_permutation_sparse(sp_mat &G, sp_mat &S, sp_mat &marker_mat, int perm_no = 100, int network_normalization_method = 0, double post_alpha = 0.85, int thread_no = 0);

  field<mat> aggregate_genesets_vision(sp_mat &G, sp_mat &S, sp_mat &marker_mat, int network_normalization_method = 0, double alpha = 0.85, int thread_no = 0);

  mat oneHot_encoding(vec batches);
  mat run_harmony(mat &X, mat &W0, vec batch, int clustering_algorithm = 1, double sigma_val = 0.1, double theta_val = 2.0, int max_iter_cluster = 200, int max_iter_harmony = 10, double eps_cluster = 1e-5, double eps_harmony = 1e-4, double tau = 0, double block_size = 0.05, double lambda_val = 1.0, bool verbose = true, int seed = 0);
  mat assess_label_enrichment(sp_mat &H, mat &M, int thread_no = 0);

} // namespace ACTIONet

#endif
