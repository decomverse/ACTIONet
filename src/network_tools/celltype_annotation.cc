#include <ACTIONet.h>

namespace ACTIONet
{
  sp_mat scale_expression(sp_mat &S)
  {
    sp_mat T = S;

    sp_mat::iterator it = T.begin();
    sp_mat::iterator it_end = T.end();

    vec mu = vec(sum(T, 1)) / vec(sum(spones(T), 1));
    for (; it != it_end; ++it)
    {
      (*it) -= mu(it.row());
    }
    vec sigma = vec(sum(square(T), 1));

    T = S;
    for (; it != it_end; ++it)
    {
      (*it) /= sigma(it.row());
    }

    return (T);
  }

  mat compute_marker_aggregate_stats_basic_sum(sp_mat &S, sp_mat &marker_mat)
  {
    marker_mat = normalise(marker_mat, 1, 0);
    sp_mat X = trans(marker_mat);

    S = scale_expression(S);
    mat stats = mat(trans(X * S));

    return (stats);
  }

  mat compute_marker_aggregate_stats_basic_sum_perm(sp_mat &S, sp_mat &marker_mat, int perm_no = 100, int thread_no = 0)
  {
    marker_mat = normalise(marker_mat, 1, 0);
    mat X = trans(mat(marker_mat));

    // S = scale_expression(S);
    mat stats = mat(trans(sp_mat(X * S)));

    int N = X.n_cols;

    mat E = zeros(size(stats));
    mat Esq = zeros(size(stats));
    parallelFor(
        0, perm_no, [&](size_t i)
        {
                  uvec perm = randperm(N);
                  mat rand_stats = mat(trans(sp_mat(X.cols(perm) * S)));
                  mat shifted_vals = (rand_stats - stats);
                  E += shifted_vals;
                  Esq += square(shifted_vals); },
        thread_no);
    mat mu = E / perm_no + stats;
    mat sigma = sqrt((Esq - square(E) / perm_no) / (perm_no - 1));
    mat Z = (stats - mu) / sigma;

    return (Z);
  }

  mat compute_marker_aggregate_stats_basic_sum_perm_smoothed(sp_mat &G, sp_mat &S, sp_mat &marker_mat, double alpha = 0.85, int max_it = 5, int perm_no = 100, int thread_no = 0)
  {
    marker_mat = normalise(marker_mat, 1, 0);
    mat X = trans(mat(marker_mat));

    S = scale_expression(S);
    sp_mat raw_stats = trans(sp_mat(X * S));
    mat stats = compute_network_diffusion_fast(G, raw_stats, thread_no, alpha, max_it); // * diagmat(vec(trans(sum(raw_stats))));

    int N = X.n_cols;

    mat E = zeros(size(stats));
    mat Esq = zeros(size(stats));
    parallelFor(
        0, perm_no, [&](size_t i)
        {
                  uvec perm = randperm(N);
                  sp_mat raw_rand_stats = trans(sp_mat(X.cols(perm) * S));
                  mat rand_stats = compute_network_diffusion_fast(G, raw_rand_stats, 1, alpha, max_it); // * diagmat(vec(trans(sum(raw_rand_stats))));

                  mat shifted_vals = (rand_stats - stats);
                  E += shifted_vals;
                  Esq += square(shifted_vals); },
        thread_no);
    mat mu = E / perm_no + stats;
    mat sigma = sqrt((Esq - square(E) / perm_no) / (perm_no - 1));
    mat Z = (stats - mu) / sigma;

    return (Z);
  }

  sp_mat LSI(sp_mat &S, double size_factor = 100000)
  {
    sp_mat X = S;

    vec col_sum_vec = zeros(X.n_cols);
    vec row_sum_vec = zeros(X.n_rows);

    sp_mat::iterator it = X.begin();
    sp_mat::iterator it_end = X.end();
    for (; it != it_end; ++it)
    {
      col_sum_vec(it.col()) += (*it);
      row_sum_vec(it.row()) += (*it);
    }

    vec kappa = size_factor / col_sum_vec;
    vec IDF = log(1 + (X.n_cols / row_sum_vec));

    for (it = X.begin(); it != X.end(); ++it)
    {
      double x = (*it) * kappa(it.col());
      x = log(1 + x) * IDF(it.row());
      *it = x;
    }

    return (X);
  }

  mat compute_marker_aggregate_stats_basic_sum_smoothed(sp_mat &G, sp_mat &S, sp_mat &marker_mat, double alpha = 0.85, int max_it = 5, int perm_no = 100, int thread_no = 0)
  {
    marker_mat = normalise(marker_mat, 1, 0);
    mat X = trans(mat(marker_mat));

    sp_mat raw_stats = trans(sp_mat(X * S));
    mat stats = compute_network_diffusion_fast(G, raw_stats, thread_no, alpha, max_it) * diagmat(vec(trans(sum(raw_stats))));

    return (stats);
  }

  mat compute_marker_aggregate_stats_basic_sum_smoothed_normalized(sp_mat &G, sp_mat &S, sp_mat &marker_mat, double alpha = 0.85, int max_it = 5, int perm_no = 100, int thread_no = 0)
  {
    marker_mat = normalise(marker_mat, 1, 0);
    mat X = trans(mat(marker_mat));

    sp_mat raw_stats = trans(sp_mat(X * S));
    mat stats = compute_network_diffusion_fast(G, raw_stats, thread_no, alpha, max_it) * diagmat(vec(trans(sum(raw_stats))));

    sp_mat p = trans(sum(S));
    vec pr =
        compute_network_diffusion_fast(G, p, thread_no, alpha, max_it).col(0);

    for (int j = 0; j < stats.n_cols; j++)
    {
      vec ppr = stats.col(j);
      vec scores_norm = log2(ppr / pr);
      uvec zero_idx = find(ppr == 0);
      scores_norm(zero_idx).zeros();
      scores_norm = scores_norm % ppr;

      stats.col(j) = scores_norm;
    }

    return (stats);
  }

  mat compute_marker_aggregate_stats_TFIDF_sum_smoothed(sp_mat &G, sp_mat &S, sp_mat &marker_mat, double alpha = 0.85, int max_it = 5, int perm_no = 100, int thread_no = 0, int normalization = 1)
  {
    marker_mat = normalise(marker_mat, 1, 0);
    mat X = trans(mat(marker_mat));

    sp_mat T;
    if (normalization == 0)
    {
      T = S;
    }
    else if (normalization == 1)
    {
      T = LSI(S);
    }

    vec base = vec(trans(T.row(0)));

    sp_mat::iterator it = T.begin();
    sp_mat::iterator it_end = T.end();
    vec E = zeros(T.n_cols);
    vec Esq = zeros(T.n_cols);
    for (; it != it_end; ++it)
    {
      double x = *it - base(it.col());
      E(it.col()) += x;
      Esq(it.col()) += (x * x);
    }
    mat mu = E / T.n_rows + base;
    mat sigma = sqrt((Esq - square(E) / T.n_rows) / (T.n_rows - 1));

    vec w1 = vec(trans(sum(marker_mat, 0)));
    vec w2 = sqrt(vec(trans(sum(square(marker_mat), 0))));

    sp_mat raw_stats = trans(sp_mat(X * T));
    mat stats;
    if (alpha == 0)
    {
      stats = raw_stats;
    }
    else
    {
      stats = compute_network_diffusion_fast(G, raw_stats, thread_no, alpha, max_it) * diagmat(vec(trans(sum(raw_stats))));
    }

    for (int i = 0; i < stats.n_rows; i++)
    {
      for (int j = 0; j < stats.n_cols; j++)
      {
        double stat = stats(i, j);
        double z = (stat - mu(i) * w1(j)) / (sigma(i) * w2(j));
        stats(i, j) = z;
      }
    }

    return (stats);
  }

  mat compute_marker_aggregate_stats_basic_sum_perm_smoothed_v2(sp_mat &G, sp_mat &S, sp_mat &marker_mat, double alpha = 0.85, int max_it = 5, int perm_no = 100, int thread_no = 0)
  {
    marker_mat = normalise(marker_mat, 1, 0);
    mat X = trans(mat(marker_mat));

    sp_mat raw_stats = trans(sp_mat(X * S));
    mat stats = compute_network_diffusion_fast(G, raw_stats, thread_no, alpha, max_it) * diagmat(vec(trans(sum(raw_stats))));

    mat raw_stats_mat = mat(raw_stats);

    int N = X.n_cols;
    mat E = zeros(size(stats));
    mat Esq = zeros(size(stats));
    parallelFor(
        0, perm_no, [&](size_t i)
        {
                  uvec perm = randperm(N);

                  sp_mat raw_rand_stats = sp_mat(raw_stats_mat.rows(perm));
                  mat rand_stats = compute_network_diffusion_fast(G, raw_rand_stats, 1, alpha, max_it) * diagmat(vec(trans(sum(raw_rand_stats))));

                  E += rand_stats;
                  Esq += square(rand_stats); },
        thread_no);
    mat mu = E / perm_no;
    mat sigma = sqrt(Esq / perm_no - square(mu));
    mat Z = (stats - mu) / sigma;

    return (Z);
  }

  mat compute_marker_aggregate_stats(sp_mat &G, sp_mat &S, sp_mat &marker_mat,
                                     double alpha = 0.85, int max_it = 5,
                                     int thread_no = 0, bool ignore_baseline_expression = false)
  {
    mat stats = zeros(S.n_cols, marker_mat.n_cols);

    int n = G.n_rows;
    sp_mat o = sp_mat(ones(n, 1));
    vec pr = compute_network_diffusion_fast(G, o, thread_no, alpha, max_it).col(0);

    for (int i = 0; i < marker_mat.n_cols; i++)
    {
      int marker_count = (int)sum(sum(spones(marker_mat.col(i))));

      int idx = 0;
      vec w = zeros(marker_count);
      vec baseline = zeros(marker_count);
      sp_mat raw_expression(S.n_cols, marker_count);
      for (sp_mat::col_iterator it = marker_mat.begin_col(i);
           it != marker_mat.end_col(i); it++)
      {
        raw_expression.col(idx) = trans(S.row(it.row()));
        w(idx) = (*it);
        baseline(idx) = accu(raw_expression.col(idx));
        idx++;
      }
      if (!ignore_baseline_expression)
      {
        baseline = baseline / sum(baseline);
        w = w % baseline;
      }
      w = w / sqrt(sum(square(w)));

      mat imputed_expression = compute_network_diffusion_fast(
          G, raw_expression, thread_no, alpha, max_it);

      for (int j = 0; j < imputed_expression.n_cols; j++)
      {
        vec ppr = imputed_expression.col(j);
        vec scores = log2(ppr / pr);
        uvec zero_idx = find(ppr == 0);
        scores(zero_idx).zeros();
        scores = scores % ppr;

        stats.col(i) += w(j) * scores;
      }
    }

    return (stats);
  }

  mat compute_marker_aggregate_stats_nonparametric(mat &S, sp_mat &marker_mat, int thread_no)
  {
    mat St = trans(S);
    mat Z = RIN_transform(St, thread_no); // cell x gene

    mat stats = zeros(Z.n_rows, marker_mat.n_cols);
    for (int i = 0; i < marker_mat.n_cols; i++)
    {
      vec v = vec(marker_mat.col(i));
      uvec idx = find(v != 0);
      vec w = v(idx);
      double sigma = sqrt(sum(square(w)));
      stats.col(i) = sum(Z.cols(idx), 1) / sigma;
    }

    return (stats);
  }

  mat compute_markers_eigengene(mat &S, sp_mat &marker_mat, int normalization, int thread_no)
  {
    mat St = trans(S); // cell x gene

    mat Z;
    if (normalization == 0)
    {
      Z = zscore(St, thread_no);
    }
    else if (normalization == 1)
    {
      Z = RIN_transform(St, thread_no);
    }
    else // default to z-score
    {
      Z = zscore(St, thread_no);
    }

    mat stats = zeros(Z.n_rows, marker_mat.n_cols);
    parallelFor(
        0, marker_mat.n_cols, [&](size_t i)
        {
                  vec v = vec(marker_mat.col(i));
                  uvec idx = find(v != 0);
                  vec w = v(idx);
                  mat subZ = Z.cols(idx);
                  subZ.each_row() %= trans(w);
                  double denom = sqrt(sum(sum(cov(subZ))));
                  vec z = sum(subZ, 1) / denom;

                  field<mat> SVD_results = HalkoSVD(subZ, 1, 5, 0, 0);
                  vec u = SVD_results(0);
                  if (dot(u, z) < 0) // orient
                  {
                    u = -u;
                  }

                  u = u * stddev(z) / stddev(u);

                  stats.col(i) = u; },
        thread_no);

    return (stats);
  }

  sp_mat normalize_expression_profile(sp_mat &S, int normalization = 1)
  {
    sp_mat T;
    if (normalization == 0)
    {
      // No normalization
      T = S;
    }
    else if (normalization == 1)
    {
      // LSI normalization
      T = LSI(S);
    }

    return (T);
  }

  mat aggregate_genesets(sp_mat &G, sp_mat &S, sp_mat &marker_mat, int network_normalization_method, int expression_normalization_method, int gene_scaling_method, double diffusion_alpha, int thread_no)
  {
    if (S.n_rows != marker_mat.n_rows)
    {
      stderr_printf("Number of genes in the expression matrix (S) and marker matrix (marker_mat) do not match\n");
      FLUSH;
      return (mat());
    }
    if (S.n_cols != G.n_rows)
    {
      stderr_printf("Number of cell in the expression matrix (S) and cell network (G) do not match\n");
      FLUSH;
      return (mat());
    }

    sp_mat markers_mat_bin = spones(marker_mat);
    vec marker_counts = vec(trans(sum(markers_mat_bin)));

    // 0: no normalization, 1: TF/IDF
    sp_mat T = normalize_expression_profile(S, expression_normalization_method);

    // 0: pagerank, 2: sym_pagerank
    sp_mat P = normalize_adj(G, network_normalization_method);

    mat marker_stats(T.n_cols, marker_mat.n_cols);
    for (int j = 0; j < marker_mat.n_cols; j++)
    {
      mat marker_expr(T.n_cols, marker_counts(j));

      int idx = 0;
      for (sp_mat::col_iterator it = marker_mat.begin_col(j); it != marker_mat.end_col(j); it++)
      {
        double w = (*it);
        marker_expr.col(idx) = w * vec(trans(T.row(it.row())));
        idx++;
      }

      // 0: no normalization, 1: z-score, 2: RINT, 3: robust z-score
      mat marker_expr_scaled = normalize_scores(marker_expr, gene_scaling_method, thread_no);
      mat marker_expr_imputed = compute_network_diffusion_Chebyshev(P, marker_expr_scaled, thread_no);

      mat Sigma = cov(marker_expr_imputed);
      double norm_factor = sqrt(sum(Sigma.diag()));

      vec aggr_stats = sum(marker_expr_imputed, 1); // each column is a marker gene
      aggr_stats = aggr_stats / norm_factor;
      marker_stats.col(j) = aggr_stats;
    }
    mat marker_stats_smoothed = compute_network_diffusion_Chebyshev(P, marker_stats, thread_no);

    return (marker_stats_smoothed);
  }

  double F2z(double F, double d1, double d2)
  {
    double mu = d2 / (d2 - 2);                                                               // Only valud if d2 > 2
    double sigma_sq = (2 * d2 * d2 * (d1 + d2 - 2)) / (d1 * (d2 - 2) * (d2 - 2) * (d2 - 4)); // Only valid when d2 > 4

    double z = (F - mu) / sqrt(sigma_sq);
    return (z);
  }

  mat aggregate_genesets_mahalanobis_2archs(sp_mat &G, sp_mat &S, sp_mat &marker_mat, int network_normalization_method, int expression_normalization_method, int gene_scaling_method, double pre_alpha, double post_alpha, int thread_no)
  {
    if (S.n_rows != marker_mat.n_rows)
    {
      stderr_printf("Number of genes in the expression matrix (S) and marker matrix (marker_mat) do not match\n");
      FLUSH;
      return (mat());
    }
    if (S.n_cols != G.n_rows)
    {
      stderr_printf("Number of cell in the expression matrix (S) and cell network (G) do not match\n");
      FLUSH;
      return (mat());
    }

    // 0: pagerank, 2: sym_pagerank
    sp_mat P;
    if (pre_alpha != 0 || post_alpha != 0)
    {
      P = normalize_adj(G, network_normalization_method);
    }

    // 0: no normalization, 1: TF/IDF
    mat T = mat(normalize_expression_profile(S, expression_normalization_method));

    if (pre_alpha != 0)
    {
      mat T_t = trans(T);
      T = compute_network_diffusion_Chebyshev(P, T_t, thread_no, pre_alpha);
      T = trans(T);
    }

    mat marker_stats(T.n_cols, marker_mat.n_cols);
    parallelFor(
        0, marker_mat.n_cols, [&](int j)
        {

      vec w = vec(marker_mat.col(j));
      uvec nnz_idx = find(w != 0);
      if(nnz_idx.n_elem != 0) {

        mat T_scaled = T.rows(nnz_idx);
        //0: no normalization, 1: z-score, 2: RINT, 3: robust z-score
        if(gene_scaling_method != 0) {
          T_scaled = normalize_scores(T_scaled, gene_scaling_method, thread_no);
        }
        T_scaled = T_scaled.each_col() % w(nnz_idx);

        uvec idx(2);
        rowvec ss = sum(T_scaled);
        idx(0) = index_min(ss);
        idx(1) = index_max(ss);

        mat W0 = T_scaled.cols(idx);

        field<mat> AA_res = run_AA(T_scaled, W0, 100);
        mat C = AA_res(0);
        mat H = AA_res(1);
        mat W = T_scaled * C;
        uword selected_arch0 = index_min(sum(W));
        uword selected_arch1 = index_max(sum(W));
        vec mu = W.col(selected_arch0);

        double p = T_scaled.n_rows;
        double n = T_scaled.n_cols;

        mat Delta = T_scaled.each_col() - mu;

        mat sigma = Delta * trans(Delta) / (n-1);
        mat sigma_inv = pinv(sigma);

        for(int k = 0; k < n; k++) {
          vec delta = Delta.col(k);
          double dist = dot(delta, sigma_inv * delta);
          double z = (dist - p) / sqrt(2*p);
          z = z < 0? 0:z;

          marker_stats(k, j) = sign(mean(delta)) * z;
        }
      } },
        thread_no);

    marker_stats.replace(datum::nan, 0);

    mat marker_stats_smoothed = marker_stats; // zscore(marker_stats, thread_no);
    if (post_alpha != 0)
    {
      stdout_printf("Post-smoothing expression values ... ");
      marker_stats_smoothed = compute_network_diffusion_Chebyshev(P, marker_stats_smoothed, thread_no, post_alpha);
      stdout_printf("done\n");
      FLUSH;
    }

    return (marker_stats_smoothed);
  }

  mat aggregate_genesets_mahalanobis_2gmm(sp_mat &G, sp_mat &S, sp_mat &marker_mat, int network_normalization_method, int expression_normalization_method, int gene_scaling_method, double pre_alpha, double post_alpha, int thread_no)
  {
    if (S.n_rows != marker_mat.n_rows)
    {
      stderr_printf("Number of genes in the expression matrix (S) and marker matrix (marker_mat) do not match\n");
      FLUSH;
      return (mat());
    }
    if (S.n_cols != G.n_rows)
    {
      stderr_printf("Number of cell in the expression matrix (S) and cell network (G) do not match\n");
      FLUSH;
      return (mat());
    }

    // 0: pagerank, 2: sym_pagerank
    sp_mat P;
    if (pre_alpha != 0 || post_alpha != 0)
    {
      P = normalize_adj(G, network_normalization_method);
    }

    // 0: no normalization, 1: TF/IDF
    mat T = mat(normalize_expression_profile(S, expression_normalization_method));

    if (pre_alpha != 0)
    {
      mat T_t = trans(T);
      T = compute_network_diffusion_Chebyshev(P, T_t, thread_no, pre_alpha);
      T = trans(T);
    }

    mat marker_stats(T.n_cols, marker_mat.n_cols);
    parallelFor(
        0, marker_mat.n_cols, [&](int j)
        {

      vec w = vec(marker_mat.col(j));
      uvec nnz_idx = find(w != 0);
      if(nnz_idx.n_elem != 0) {

        mat T_scaled = T.rows(nnz_idx);
        //0: no normalization, 1: z-score, 2: RINT, 3: robust z-score
        if(gene_scaling_method != 0) {
          T_scaled = normalize_scores(T_scaled, gene_scaling_method, thread_no);
        }
        T_scaled = T_scaled.each_col() % w(nnz_idx);

        gmm_full model;

        bool status = model.learn(T_scaled, 2, maha_dist, static_spread, 10, 10, 1e-10, false);
        mat W = model.means;

        uword selected_arch0 = index_min(sum(W));
        uword selected_arch1 = index_max(sum(W));
        vec mu = W.col(selected_arch0);

        double p = T_scaled.n_rows;
        double n = T_scaled.n_cols;

        mat Delta = T_scaled.each_col() - mu;

        //mat sigma = cov(trans(T_scaled));
        mat sigma = Delta * trans(Delta) / (n-1);
        mat sigma_inv = pinv(sigma);

        for(int k = 0; k < n; k++) {
          vec delta = Delta.col(k);
          double dist = dot(delta, sigma_inv * delta);
          double z = (dist - p) / sqrt(2*p);
          z = z < 0? 0:z;

          marker_stats(k, j) = sign(mean(delta)) * z;
        }
      } },
        thread_no);

    marker_stats.replace(datum::nan, 0);

    mat marker_stats_smoothed = marker_stats; // zscore(marker_stats, thread_no);
    if (post_alpha != 0)
    {
      stdout_printf("Post-smoothing expression values ... ");
      marker_stats_smoothed = compute_network_diffusion_Chebyshev(P, marker_stats_smoothed, thread_no, post_alpha);
      stdout_printf("done\n");
      FLUSH;
    }

    return (marker_stats_smoothed);
  }

  mat doubleNorm(mat &X)
  {
    vec rs = sum(X, 1);
    vec cs = trans(sum(X, 0));

    mat Dr = diagmat(1 / sqrt(rs));
    mat Dc = diagmat(1 / sqrt(cs));

    mat Y = Dr * X * Dc;

    return (Y);
  }

  mat aggregate_genesets_weighted_enrichment(sp_mat &G, sp_mat &S, sp_mat &marker_mat, int network_normalization_method, int expression_normalization_method, int gene_scaling_method, double pre_alpha, double post_alpha, int thread_no)
  {
    if (S.n_rows != marker_mat.n_rows)
    {
      stdout_printf("Number of genes in the expression matrix (S) and marker matrix (marker_mat) do not match\n");
      FLUSH;
      return (mat());
    }
    if (S.n_cols != G.n_rows)
    {
      stdout_printf("Number of cell in the expression matrix (S) and cell network (G) do not match\n");
      FLUSH;
      return (mat());
    }

    // 0: pagerank, 2: sym_pagerank
    sp_mat P;
    if (pre_alpha != 0 || post_alpha != 0)
    {
      P = normalize_adj(G, network_normalization_method);
    }

    // 0: no normalization, 1: TF/IDF
    mat T = mat(normalize_expression_profile(S, expression_normalization_method));

    if (pre_alpha != 0)
    {
      mat T_t = trans(T);
      T = compute_network_diffusion_Chebyshev(P, T_t, thread_no, pre_alpha);
      T = trans(T);
    }

    if (gene_scaling_method != 0)
    {
      if (gene_scaling_method > 0)
      {
        T = normalise(T, gene_scaling_method, 1);
      }
      else
      {
        mat T_t = trans(T);
        T = normalize_scores(T_t, -gene_scaling_method, thread_no);
        T = trans(T);
      }
    }

    mat marker_stats;
    if (gene_scaling_method >= 0)
    {
      field<mat> res = assess_enrichment(T, marker_mat, thread_no);
      marker_stats = trans(res(0));
    }
    else
    {
      vec w = vec(sqrt(trans(sum(square(marker_mat), 0))));
      w.replace(0.0, 1.0);
      mat w_mat = diagmat(1.0 / w);
      sp_mat marker_mat_t = trans(marker_mat);
      marker_stats = trans(spmat_mat_product_parallel(marker_mat_t, T, thread_no)) * w_mat;
    }
    marker_stats.replace(datum::nan, 0);

    mat marker_stats_smoothed = marker_stats; // zscore(marker_stats, thread_no);
    if (post_alpha != 0)
    {
      stdout_printf("Post-smoothing expression values ... ");
      marker_stats_smoothed = compute_network_diffusion_Chebyshev(P, marker_stats_smoothed, thread_no, post_alpha);
      stdout_printf("done\n");
      FLUSH;
    }
    /*
        mat marker_stats_smoothed = zscore(marker_stats); // zscore(marker_stats, thread_no);
        field<vec> auto_out = autocorrelation_Moran(G, marker_stats_smoothed, 0, 0, thread_no);
        if (post_alpha != 0)
        {
          stdout_printf("Post-smoothing expression values ... ");
          marker_stats_smoothed = compute_network_diffusion_Chebyshev(P, marker_stats_smoothed, thread_no, post_alpha);
          stdout_printf("done\n");
          FLUSH;
        }
        vec locality = auto_out(0);
        locality.replace(datum::nan, 0);
        locality /= mean(locality);
        marker_stats_smoothed = marker_stats_smoothed * diagmat(locality);
    */

    return (marker_stats_smoothed);
  }

  mat aggregate_genesets_weighted_enrichment_permutation(sp_mat &G, sp_mat &S, sp_mat &marker_mat, int network_normalization_method, int expression_normalization_method, int gene_scaling_method, double pre_alpha, double post_alpha, int thread_no, int perm_no)
  {
    if (S.n_rows != marker_mat.n_rows)
    {
      stderr_printf("Number of genes in the expression matrix (S) and marker matrix (marker_mat) do not match\n");
      FLUSH;
      return (mat());
    }
    if (S.n_cols != G.n_rows)
    {
      stderr_printf("Number of cell in the expression matrix (S) and cell network (G) do not match\n");
      FLUSH;
      return (mat());
    }

    // 0: pagerank, 2: sym_pagerank
    sp_mat P;
    if (pre_alpha != 0 || post_alpha != 0)
    {
      P = normalize_adj(G, network_normalization_method);
    }

    // 0: no normalization, 1: TF/IDF
    mat T = mat(normalize_expression_profile(S, expression_normalization_method));

    if (pre_alpha != 0)
    {
      mat T_t = trans(T);
      T = compute_network_diffusion_Chebyshev(P, T_t, thread_no, pre_alpha);
      T = trans(T);
    }

    if (gene_scaling_method != 0)
    {
      if (gene_scaling_method > 0)
      {
        T = normalise(T, gene_scaling_method, 1);
      }
      else
      {
        mat T_t = trans(T);
        T = normalize_scores(T_t, -gene_scaling_method, thread_no);
        T = trans(T);
      }
    }

    sp_mat X = trans(marker_mat);

    mat Y = T;

    mat stats = spmat_mat_product(X, Y);

    mat E = zeros(size(stats));
    mat Esq = zeros(size(stats));
    for (int k = 0; k < perm_no; k++)
    {
      uvec perm = randperm(Y.n_rows);
      mat Y_perm = Y.rows(perm);
      mat rand_stats = spmat_mat_product(X, Y_perm);

      mat delta = (rand_stats - stats);
      E += delta;
      Esq += square(delta);
    }
    mat mu = stats + E / perm_no;
    mat sigma = sqrt((Esq - square(E) / perm_no) / (perm_no - 1));
    mat marker_stats = trans((stats - mu) / sigma);

    marker_stats.replace(datum::nan, 0);

    mat marker_stats_smoothed = marker_stats; // zscore(marker_stats, thread_no);
    if (post_alpha != 0)
    {
      stdout_printf("Post-smoothing expression values ... ");
      marker_stats_smoothed = compute_network_diffusion_Chebyshev(P, marker_stats_smoothed, thread_no, post_alpha);
      stdout_printf("done\n");
      FLUSH;
    }

    /*
        mat marker_stats_smoothed = zscore(marker_stats); // zscore(marker_stats, thread_no);
        field<vec> auto_out = autocorrelation_Moran(G, marker_stats_smoothed, 0, 0, thread_no);
        if (post_alpha != 0)
        {
          stdout_printf("Post-smoothing expression values ... ");
          marker_stats_smoothed = compute_network_diffusion_Chebyshev(P, marker_stats_smoothed, thread_no, post_alpha);
          stdout_printf("done\n");
          FLUSH;
        }
        vec locality = auto_out(0);
        locality.replace(datum::nan, 0);
        locality /= mean(locality);
        marker_stats_smoothed = marker_stats_smoothed * diagmat(locality);
    */

    return (marker_stats_smoothed);
  }

  field<mat> aggregate_genesets_vision(sp_mat &G, sp_mat &S, sp_mat &marker_mat, int network_normalization_method, double alpha, int thread_no)
  {
    field<mat> out(3);

    if (S.n_rows != marker_mat.n_rows)
    {
      stderr_printf("Number of genes in the expression matrix (S) and marker matrix (marker_mat) do not match\n");
      FLUSH;
      return (out);
    }
    if (S.n_cols != G.n_rows)
    {
      stderr_printf("Number of cell in the expression matrix (S) and cell network (G) do not match\n");
      FLUSH;
      return (out);
    }

    // 0: pagerank, 2: sym_pagerank
    sp_mat P;
    if (alpha != 0)
    {
      P = normalize_adj(G, network_normalization_method);
    }

    mat X = mat(marker_mat);
    sp_mat St = trans(S);

    mat stats = spmat_mat_product_parallel(St, X, thread_no);

    // Compute cell-specific stats to adjust for depth, etc.
    vec mu = zeros(S.n_cols);
    vec nnz = zeros(S.n_cols);
    for (sp_mat::const_iterator it = S.begin(); it != S.end(); ++it)
    {
      mu[it.col()] += (*it);
      nnz[it.col()]++;
    }
    mu /= S.n_rows;
    vec p_nnz = nnz / S.n_rows; // 1 - p_zero

    vec sigma_sq = zeros(S.n_cols);
    for (sp_mat::const_iterator it = S.begin(); it != S.end(); ++it)
    {
      float delta = mu[it.col()] - (*it);
      sigma_sq[it.col()] += delta * delta;
    }
    sigma_sq += (S.n_rows * (1 - p_nnz)) % square(mu); // Adjust for nnzs
    sigma_sq /= (S.n_rows - 1);

    // Standardize using sampling mean (from Vision: https://www.nature.com/articles/s41467-019-12235-0)
    rowvec k1 = rowvec(sum(marker_mat));
    rowvec k2 = rowvec(sum(square(marker_mat)));

    mat sampling_mu = mu * k1;
    mat sampling_sigma_sq = sigma_sq * k2;
    mat marker_stats = (stats - sampling_mu) / sqrt(sampling_sigma_sq);
    marker_stats.replace(datum::nan, 0);

    mat marker_stats_smoothed = marker_stats;
    if (alpha != 0)
    {
      stdout_printf("Smoothing geneset scores ... ");
      marker_stats_smoothed = compute_network_diffusion_Chebyshev(P, marker_stats_smoothed, thread_no, alpha);
      stdout_printf("done\n");
      FLUSH;
    }

    out(0) = marker_stats_smoothed;
    out(1) = marker_stats;
    out(2) = stats;

    return (out);
  }

} // namespace ACTIONet
