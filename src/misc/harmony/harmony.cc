#include "ACTIONet.h"

typedef arma::mat MATTYPE;
typedef arma::vec VECTYPE;
typedef arma::rowvec ROWVECTYPE;
typedef arma::cube CUBETYPE;

namespace ACTIONet
{
  void cosine_normalize(MATTYPE &X, int margin, bool do_safe)
  {
    // to avoid Inf values, first divide by max
    if (margin == 1)
    {
      for (unsigned r = 0; r < X.n_rows; r++)
      {
        if (do_safe)
          X.row(r) = X.row(r) / X.row(r).max();
        X.row(r) = X.row(r) / norm(X.row(r), 2);
      }
    }
    else
    {
      for (unsigned c = 0; c < X.n_cols; c++)
      {
        if (do_safe)
          X.col(c) = X.col(c) / X.col(c).max();
        X.col(c) = X.col(c) / norm(X.col(c), 2);
      }
    }
  }

  MATTYPE safe_entropy(const MATTYPE &X)
  {
    MATTYPE A = X % log(X);
    A.elem(find_nonfinite(A)).zeros();
    return (A);
  }

  // Overload pow to work on a MATTYPErix and vector
  MATTYPE harmony_pow(MATTYPE A, const VECTYPE &T)
  {
    for (unsigned c = 0; c < A.n_cols; c++)
    {
      A.col(c) = pow(A.col(c), as_scalar(T.row(c)));
    }
    return (A);
  }

  MATTYPE merge_R(const MATTYPE &R, float thresh = 0.8)
  {
    MATTYPE cor_res = cor(R.t());
    int K = R.n_rows;

    // define equivalence classes
    uvec equiv_classes = linspace<uvec>(0, K - 1, K);
    int new_idx;
    for (int i = 0; i < K - 1; i++)
    {
      for (int j = i + 1; j < K; j++)
      {
        if (cor_res(i, j) > thresh)
        {
          new_idx = min(equiv_classes(i), equiv_classes(j));
          equiv_classes(i) = new_idx;
          equiv_classes(j) = new_idx;
        }
      }
    }

    // sum over equivalence classes
    uvec uclasses = unique(equiv_classes);
    MATTYPE R_new = zeros<MATTYPE>(uclasses.n_elem, R.n_cols);
    for (unsigned i = 0; i < R_new.n_rows; i++)
    {
      uvec idx = find(equiv_classes == uclasses(i));
      R_new.row(i) = sum(R.rows(idx), 0);
    }
    return R_new;
  }

  MATTYPE compute_Y(const MATTYPE &Z_cos, const MATTYPE &R)
  {
    return arma::normalise(Z_cos * R.t(), 2, 0);
  }

  MATTYPE scaleRows_dgc(const VECTYPE &x, const VECTYPE &p, const VECTYPE &i,
                        int ncol, int nrow, float thresh)
  {
    // (0) fill in non-zero elements
    MATTYPE res = arma::zeros<MATTYPE>(nrow, ncol);
    for (int c = 0; c < ncol; c++)
    {
      for (int j = p[c]; j < p[c + 1]; j++)
      {
        res(i[j], c) = x(j);
      }
    }

    // (1) compute means
    VECTYPE mean_vec = arma::zeros<VECTYPE>(nrow);
    for (int c = 0; c < ncol; c++)
    {
      for (int j = p[c]; j < p[c + 1]; j++)
      {
        mean_vec(i[j]) += x[j];
      }
    }
    mean_vec /= ncol;

    // (2) compute SDs
    VECTYPE sd_vec = arma::zeros<VECTYPE>(nrow);
    arma::uvec nz = arma::zeros<arma::uvec>(nrow);
    nz.fill(ncol);
    for (int c = 0; c < ncol; c++)
    {
      for (int j = p[c]; j < p[c + 1]; j++)
      {
        sd_vec(i[j]) += (x[j] - mean_vec(i[j])) * (x[j] - mean_vec(i[j])); // (x - mu)^2
        nz(i[j])--;
      }
    }

    // count for the zeros
    for (int r = 0; r < nrow; r++)
    {
      sd_vec(r) += nz(r) * mean_vec(r) * mean_vec(r);
    }

    sd_vec = arma::sqrt(sd_vec / (ncol - 1));

    // (3) scale values
    res.each_col() -= mean_vec;
    res.each_col() /= sd_vec;
    res.elem(find(res > thresh)).fill(thresh);
    res.elem(find(res < -thresh)).fill(-thresh);
    return res;
  }

  class harmony
  {
  public:
    /* CONSTRUCTORS etc */
    harmony(int __K);
    void setup(MATTYPE &__Z, MATTYPE &__Phi, MATTYPE &__Phi_moe, VECTYPE __Pr_b,
               VECTYPE __sigma, VECTYPE __theta, int __max_iter_kmeans,
               float __epsilon_kmeans, float __epsilon_harmony,
               int __K, float tau, float __block_size,
               MATTYPE __lambda, bool __verbose, int seed);

    /* METHODS */
    void moe_correct_ridge_cpp();
    CUBETYPE moe_ridge_get_betas_cpp();
    void init_cluster_cpp(unsigned C);
    int cluster_harmony();
    int cluster_SPA();
    int cluster_AA();

    void allocate_buffers();
    void compute_objective();
    int update_R();
    bool check_convergence(int type);
    std::mt19937_64 engine;

    /* FIELDS */
    MATTYPE R, Z_orig, Z_corr, Z_cos, Y, Y_unnormed, Phi, Phi_moe;
    VECTYPE Pr_b, theta, N_b, sigma, sigma_prior;
    MATTYPE lambda; // diagonal MATTYPErix of ridge regression penalties
    vector<float> objective_harmony;
    vector<float> objective_kmeans, objective_kmeans_dist, objective_kmeans_entropy, objective_kmeans_cross;
    vector<int> kmeans_rounds; // OLD: Kb

    //    vector<uvec> phi_map;
    float block_size, epsilon_kmeans, epsilon_harmony, merge_thresh_global;
    int N, K, B, d, max_iter_kmeans, window_size;

    // buffers
    MATTYPE _scale_dist, dist_mat, O, E, dir_prior, Phi_Rk; // N_k, N_kb, N_b, numerator, denominator, C;
    uvec update_order, cells_update;
    MATTYPE W;

    // flags
    bool ran_setup, ran_init, verbose; // do_merge_R;
  };

  harmony::harmony(int __K) : K(__K) {}

  void harmony::setup(MATTYPE &__Z, MATTYPE &__Phi, MATTYPE &__Phi_moe, VECTYPE __Pr_b,
                      VECTYPE __sigma, VECTYPE __theta, int __max_iter_kmeans,
                      float __epsilon_kmeans, float __epsilon_harmony,
                      int __K, float tau, float __block_size,
                      MATTYPE __lambda, bool __verbose, int seed)
  {
    engine.seed(seed);

    Z_corr = MATTYPE(__Z);
    Z_orig = MATTYPE(__Z);
    Z_cos = MATTYPE(Z_orig);
    cosine_normalize(Z_cos, 0, true); // normalize columns

    Phi = __Phi;
    Phi_moe = __Phi_moe;
    N = Z_corr.n_cols;
    Pr_b = __Pr_b;
    B = Phi.n_rows;
    d = Z_corr.n_rows;
    window_size = 3;
    epsilon_kmeans = __epsilon_kmeans;
    epsilon_harmony = __epsilon_harmony;

    lambda = __lambda;
    sigma = __sigma;
    sigma_prior = __sigma;
    block_size = __block_size;
    K = __K;
    max_iter_kmeans = __max_iter_kmeans;
    verbose = __verbose;

    theta = __theta;
    allocate_buffers();
    ran_setup = true;
  }

  void harmony::allocate_buffers()
  {
    _scale_dist = zeros<MATTYPE>(K, N);
    dist_mat = zeros<MATTYPE>(K, N);
    O = zeros<MATTYPE>(K, B);
    E = zeros<MATTYPE>(K, B);
    W = zeros<MATTYPE>(B + 1, d);
    Phi_Rk = zeros<MATTYPE>(B + 1, N);
  }

  void harmony::init_cluster_cpp(unsigned C)
  {
    // kmeans is called outside, in the R function
    cosine_normalize(Y, 0, false); // normalize columns

    // (2) ASSIGN CLUSTER PROBABILITIES
    // using a nice property of cosine distance,
    // compute squared distance directly with cross product
    dist_mat = 2 * (1 - Y.t() * Z_cos);

    // if C > 0, only initialize the clusters not set by the user
    // with cluster_prior
    if (C > 0 && C < K)
    {
      MATTYPE Rtmp = -dist_mat.rows(C, K - 1);
      Rtmp.each_col() /= sigma.rows(C, K - 1);
      Rtmp.each_row() -= max(Rtmp, 0);
      Rtmp = exp(Rtmp);
      Rtmp.each_row() /= sum(Rtmp, 0);
      R.rows(C, K - 1) = Rtmp;
    }
    else
    {
      R = -dist_mat;
      R.each_col() /= sigma;
      R.each_row() -= max(R, 0);
      R = exp(R);
      R.each_row() /= sum(R, 0);
    }

    // (3) BATCH DIVERSITY STATISTICS
    E = sum(R, 1) * Pr_b.t();
    O = R * Phi.t();

    compute_objective();
    objective_harmony.push_back(objective_kmeans.back());
    ran_init = true;
  }

  void harmony::compute_objective()
  {
    float kmeans_error = as_scalar(accu(R % dist_mat));
    float _entropy = as_scalar(accu(safe_entropy(R).each_col() % sigma)); // NEW: vector sigma
    float _cross_entropy;
    _cross_entropy = as_scalar(accu((R.each_col() % sigma) % ((arma::repmat(theta.t(), K, 1) % log((O + 1) / (E + 1))) * Phi)));
    objective_kmeans.push_back(kmeans_error + _entropy + _cross_entropy);
    objective_kmeans_dist.push_back(kmeans_error);
    objective_kmeans_entropy.push_back(_entropy);
    objective_kmeans_cross.push_back(_cross_entropy);
  }

  bool harmony::check_convergence(int type)
  {
    float obj_new, obj_old;
    switch (type)
    {
    case 0:
      // Clustering
      // compute new window mean
      obj_old = 0;
      obj_new = 0;
      for (int i = 0; i < window_size; i++)
      {
        obj_old += objective_kmeans[objective_kmeans.size() - 2 - i];
        obj_new += objective_kmeans[objective_kmeans.size() - 1 - i];
      }
      if ((obj_old - obj_new) / abs(obj_old) < epsilon_kmeans)
      {
        return (true);
      }
      else
      {
        return (false);
      }
    case 1:
      // Harmony
      obj_old = objective_harmony[objective_harmony.size() - 2];
      obj_new = objective_harmony[objective_harmony.size() - 1];
      if ((obj_old - obj_new) / abs(obj_old) < epsilon_harmony)
      {
        return (true);
      }
      else
      {
        return (false);
      }
    }

    // gives warning if we don't give default return value
    return (true);
  }

  int harmony::cluster_SPA()
  {
    SPA_results SPA_res = run_SPA(Z_corr, Y.n_cols);
    Y = Z_corr.cols(SPA_res.selected_columns);
    R = run_simplex_regression(Y, Z_corr, false);

    dist_mat = 2 * (1 - Y.t() * Z_cos); // Y was changed

    compute_objective();
    kmeans_rounds.push_back(1);
    objective_harmony.push_back(objective_kmeans.back());

    return (0);
  }

  int harmony::cluster_AA()
  {
    field<mat> AA_res = run_AA(Z_corr, Y, 50, 1e-300);
    Y = Z_corr * AA_res(0);
    R = AA_res(1);

    dist_mat = 2 * (1 - Y.t() * Z_cos); // Y was changed

    compute_objective();
    kmeans_rounds.push_back(1);
    objective_harmony.push_back(objective_kmeans.back());

    return (0);
  }

  int harmony::cluster_harmony()
  {
    int err_status = 0;
    int iter;

    // Z_cos has changed
    // R has assumed to not change
    // so update Y to match new integrated data
    dist_mat = 2 * (1 - Y.t() * Z_cos); // Z_cos was changed
    for (iter = 0; iter < max_iter_kmeans; iter++)
    {
      // STEP 1: Update Y
      Y = compute_Y(Z_cos, R);
      dist_mat = 2 * (1 - Y.t() * Z_cos); // Y was changed

      // STEP 3: Update R
      err_status = update_R();
      if (err_status != 0)
      {
        // Rcout << "Compute R failed. Exiting from clustering." << endl;
        return err_status;
      }

      // STEP 4: Check for convergence
      compute_objective();
      if (iter > window_size)
      {
        bool convergence_status = check_convergence(0);
        if (convergence_status)
        {
          //        Rcout << "... Breaking Clustering ..., status = " << convergence_status << endl;
          iter++;
          // Rcout << "Clustered for " << iter << " iterations" << endl;
          break;
        }
      }
    }
    kmeans_rounds.push_back(iter);
    objective_harmony.push_back(objective_kmeans.back());
    return 0;
  }

  int harmony::update_R()
  {
    vec idx = linspace<vec>(0, N - 1, N);
    aarand::shuffle(idx.memptr(), idx.n_elem, engine);
    update_order = conv_to<uvec>::from(idx);

    // update_order = shuffle(linspace<uvec>(0, N - 1, N));

    _scale_dist = -dist_mat;
    _scale_dist.each_col() /= sigma; // NEW: vector sigma
    _scale_dist.each_row() -= max(_scale_dist, 0);
    _scale_dist = exp(_scale_dist);

    // GENERAL CASE: online updates, in blocks of size (N * block_size)
    int n_blocks = (int)(ceil(1.0 / block_size));
    int cells_per_block = (N / n_blocks) + 1;
    for (int i = 0; i < n_blocks; i++)
    {
      // gather cell updates indices
      int idx_min = i * cells_per_block;
      int idx_max = min(idx_min + cells_per_block, N);
      if (idx_min > idx_max)
        break;
      uvec idx_list = linspace<uvec>(idx_min, idx_max - 1, idx_max - idx_min);
      cells_update = update_order.rows(idx_list);

      // Step 1: remove cells
      E -= sum(R.cols(cells_update), 1) * Pr_b.t();
      O -= R.cols(cells_update) * Phi.cols(cells_update).t();

      // Step 2: recompute R for removed cells
      R.cols(cells_update) = _scale_dist.cols(cells_update);
      R.cols(cells_update) = R.cols(cells_update) % (harmony_pow((E + 1) / (O + 1), theta) * Phi.cols(cells_update));
      R.cols(cells_update) = normalise(R.cols(cells_update), 1, 0); // L1 norm columns

      // Step 3: put cells back
      E += sum(R.cols(cells_update), 1) * Pr_b.t();
      O += R.cols(cells_update) * Phi.cols(cells_update).t();
    }
    return 0;
  }

  void harmony::moe_correct_ridge_cpp()
  {
    Z_corr = Z_orig;
    for (int k = 0; k < K; k++)
    {
      // Phi_Rk = Phi_moe * arma::diagmat(R.row(k));
      rowvec r = R.row(k);
      mat Phi_Rk = Phi_moe.each_row() % r;

      W = arma::inv(Phi_Rk * Phi_moe.t() + lambda) * Phi_Rk * Z_orig.t();
      W.row(0).zeros(); // do not remove the intercept
      Z_corr -= W.t() * Phi_Rk;
    }
    Z_cos = arma::normalise(Z_corr, 2, 0);
  }

  CUBETYPE harmony::moe_ridge_get_betas_cpp()
  {
    CUBETYPE W_cube(W.n_rows, W.n_cols, K); // rows, cols, slices
    for (unsigned k = 0; k < K; k++)
    {
      Phi_Rk = Phi_moe * arma::diagmat(R.row(k));
      W_cube.slice(k) = arma::inv(Phi_Rk * Phi_moe.t() + lambda) * Phi_Rk * Z_orig.t();
    }
    return W_cube;
  }

  mat run_harmony(mat &X, mat &W0, vec batch, int clustering_algorithm, double sigma_val, double theta_val, int max_iter_cluster, int max_iter_harmony, double eps_cluster, double eps_harmony, double tau, double block_size, double lambda_val, bool verbose, int seed)
  {
    int nclust = W0.n_cols;

    harmony harmonyObj(0);

    mat Phi = oneHot_encoding(batch);
    mat Phi_moe = join_vert(ones(1, X.n_cols), Phi);

    double N = X.n_cols;
    vec N_b = sum(Phi, 1);
    vec Pr_b = N_b / N;
    vec sigma = sigma_val * ones(nclust);
    vec theta = theta_val * ones(Phi.n_rows);
    vec lambda = lambda_val * ones(Phi.n_rows);
    vec lambda_ext = zeros(lambda.n_elem + 1);
    lambda_ext(span(1, lambda_ext.n_elem - 1)) = lambda;
    mat lambda_mat = diagmat(lambda_ext);

    harmonyObj.setup(
        X, Phi, Phi_moe,
        Pr_b, sigma, theta, max_iter_cluster, eps_cluster,
        eps_harmony, nclust, tau, block_size, lambda_mat, verbose, seed);

    harmonyObj.Y = W0;
    harmonyObj.init_cluster_cpp(0);

    // Main loop
    for (int iter = 0; iter < max_iter_harmony; iter++)
    {
      if (verbose)
      {
        stdout_printf("Harmony %d/%d\n", iter, max_iter_harmony);
      }

      // STEP 1: do clustering
      switch (clustering_algorithm)
      {
      case 1:
        harmonyObj.cluster_harmony();
        break;
      case 2:
        harmonyObj.cluster_SPA();
        break;
      case 3:
        harmonyObj.cluster_AA();
        break;
      default:
        harmonyObj.cluster_harmony();
      }

      // STEP 2 : regress out covariates
      harmonyObj.moe_correct_ridge_cpp();

      // STEP 3 : check for convergence
      if (harmonyObj.check_convergence(1))
      {
        if (verbose)
        {
          stdout_printf("Harmony converged after %d iterations\n", iter);
          break;
        }
      }
    }

    return (harmonyObj.Z_corr);
  }
}