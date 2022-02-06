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

    //S = scale_expression(S);
    mat stats = mat(trans(sp_mat(X * S)));

    int N = X.n_cols;

    mat E = zeros(size(stats));
    mat Esq = zeros(size(stats));
    parallelFor(0, perm_no, [&] (size_t i)
                {
                  uvec perm = randperm(N);
                  mat rand_stats = mat(trans(sp_mat(X.cols(perm) * S)));
                  mat shifted_vals = (rand_stats - stats);
                  E += shifted_vals;
                  Esq += square(shifted_vals);
                }, thread_no);
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
    parallelFor(0, perm_no, [&] (size_t i)
                {
                  uvec perm = randperm(N);
                  sp_mat raw_rand_stats = trans(sp_mat(X.cols(perm) * S));
                  mat rand_stats = compute_network_diffusion_fast(G, raw_rand_stats, 1, alpha, max_it); // * diagmat(vec(trans(sum(raw_rand_stats))));

                  mat shifted_vals = (rand_stats - stats);
                  E += shifted_vals;
                  Esq += square(shifted_vals);
                }, thread_no);
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
    //mat stats = compute_network_diffusion_fast(G, raw_stats, thread_no, alpha, max_it) * diagmat(vec(trans(sum(raw_stats))));

    /*
	sp_mat p = trans(sum(S));
	vec pr =
	compute_network_diffusion_fast(G, p, thread_no, alpha, max_it).col(0);

	*/

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
    parallelFor(0, perm_no, [&] (size_t i)
                {
                  uvec perm = randperm(N);

                  sp_mat raw_rand_stats = sp_mat(raw_stats_mat.rows(perm));
                  mat rand_stats = compute_network_diffusion_fast(G, raw_rand_stats, 1, alpha, max_it) * diagmat(vec(trans(sum(raw_rand_stats))));

                  E += rand_stats;
                  Esq += square(rand_stats);
                }, thread_no);
    mat mu = E / perm_no;
    mat sigma = sqrt(Esq / perm_no - square(mu));
    mat Z = (stats - mu) / sigma;

    return (Z);
  }

  /*
mat compute_marker_aggregate_stats_basic_sum_smoothed(sp_mat &G, sp_mat &S, sp_mat &marker_mat, double alpha = 0.85, int max_it = 5, int perm_no = 100, int thread_no = 0) {
  mat X = trans(mat(marker_mat));
  X = normalize(X, 1, 0);
  
  sp_mat raw_stats = trans(sp_mat(X * S));
  mat stats = compute_network_diffusion_fast(G, raw_stats, thread_no, alpha, max_it) * diagmat(vec(trans(sum(raw_stats))));
  stats = normalize(stats, 1, 0);
 
  vec x = vec(trans(sum(S)));
  vec q = x / sum(x);

    for(int j = 0; j < stats.n_cols; j++) {
		vec p = stats.col(j);		
		
		vec scores_norm = log2(p / q);
		uvec zero_idx = find(p == 0);
		scores_norm(zero_idx).zeros();
		scores_norm = scores_norm % p;
		
		stats.col(j) = scores_norm;
	}
	 
 
  return (stats);
}
*/

  mat compute_marker_aggregate_stats(sp_mat &G, sp_mat &S, sp_mat &marker_mat,
                                     double alpha = 0.85, int max_it = 5,
                                     int thread_no = 0, bool ignore_baseline_expression = false)
  {
    mat stats = zeros(S.n_cols, marker_mat.n_cols);

    int n = G.n_rows;
    sp_mat o = sp_mat(ones(n, 1));
    // vec pr = compute_network_diffusion(G, o, thread_no, alpha, max_it).col(0);
    vec pr =
        compute_network_diffusion_fast(G, o, thread_no, alpha, max_it).col(0);

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
        baseline = baseline / sum(baseline); //sqrt(sum(square(baseline)));
        w = w % baseline;
      }
      w = w / sqrt(sum(square(w)));

      // mat imputed_expression = compute_network_diffusion(G, raw_expression,
      // thread_no, alpha, max_it);
      mat imputed_expression = compute_network_diffusion_fast(
          G, raw_expression, thread_no, alpha, max_it);

      for (int j = 0; j < imputed_expression.n_cols; j++)
      {
        vec ppr = imputed_expression.col(j);
        vec scores = log2(ppr / pr);
        uvec zero_idx = find(ppr == 0);
        scores(zero_idx).zeros();
        scores = scores % ppr;
        /*
		// Rank-based inverse normal transformation
		uvec row_perm_forward = stable_sort_index(scores);
		uvec row_perm = stable_sort_index(row_perm_forward);
		vec p = (row_perm + ones(size(row_perm))) / (row_perm.n_elem + 1);
		vec z = zeros(size(p));
		for (int j = 0; j < p.n_elem; j++) {
		  double norm_inv = r8_normal_01_cdf_inverse(p(j));
		  z(j) = norm_inv;
		}		
		*/
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
    parallelFor(0, marker_mat.n_cols, [&] (size_t i)
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

                  stats.col(i) = u;
                }, thread_no);
    /*
    mat Zt = normalise(trans(Z), 1, 0);
    mat basis = normalise(Zt * stats, 1, 0);
    mat weights = trans(run_simplex_regression(basis, Zt, false));
*/
    return (stats);
  }

/*
  mat annotate_cells_decoupled(mat &W, mat &H, sp_mat &marker_mat, int thread_no) {
    sp_mat It = trans(marker_mat);

    mat factor_stats = spmat_mat_product(It, W);

    factor_stats_z = ...
   
    mat cell_stats = factor_stats_z * H;
    
    rowvec w = sqrt(sum(square(H)));
    cell_stats.each_row() %= denom;

    return(cell_stats);
  }
*/
  sp_mat normalize_expression_profile(sp_mat &S, int normalization = 1) {
    sp_mat T;
    if (normalization == 0)
    {
      printf("No norm");
      T = S;
    }
    else if (normalization == 1)
    {
      printf("LSI");
      T = LSI(S);
    }

    return(T); 
  }

  mat aggregate_genesets(sp_mat &G, sp_mat &S, sp_mat &marker_mat, int network_normalization_method, int expression_normalization_method, int gene_scaling_method, double diffusion_alpha, int thread_no)
  {
    if(S.n_rows != marker_mat.n_rows) {
      fprintf(stderr, "Number of genes in the expression matrix (S) and marker matrix (marker_mat) do not match\n");
      return(mat());
    }
    if(S.n_cols != G.n_rows) {
      fprintf(stderr, "Number of cell in the expression matrix (S) and cell network (G) do not match\n");
      return(mat());
    }

    sp_mat markers_mat_bin = spones(marker_mat);
    vec marker_counts = vec(trans(sum(markers_mat_bin)));

    // 0: no normalization, 1: TF/IDF
    printf("Normalize expreesion profile\n");
    sp_mat T = normalize_expression_profile(S, expression_normalization_method);

    // 0: pagerank, 2: sym_pagerank
    printf("Normalize adjacency matrix\n");
    sp_mat P = normalize_adj(G, network_normalization_method);

    mat marker_stats(S.n_cols, marker_mat.n_cols);
    for(int j = 0; j < marker_mat.n_cols; j++) {
      printf("j = %d\n", j);
      mat marker_expr(S.n_cols, marker_counts(j));

      printf("subsetting expression\n");
      int idx = 0;
      for (sp_mat::col_iterator it = marker_mat.begin_col(j); it != marker_mat.end_col(j); it++) {
        double w = (*it);
        marker_expr.col(idx) = w*vec(trans(S.row(it.row())));        
        idx++;
      }

      //0: no normalization, 1: z-score, 2: RINT, 3: robust z-score
      printf("normalize_scores\n");      
      mat marker_expr_scaled = normalize_scores(marker_expr, gene_scaling_method, thread_no);

      printf("compute_network_diffusion_Chebyshev\n");
      mat marker_expr_imputed = compute_network_diffusion_Chebyshev(P, marker_expr_scaled, thread_no);

      mat Sigma = cov(marker_expr_imputed);
      //double norm_factor = sqrt(sum(sum(Sigma)));
      double norm_factor = sqrt(sum(Sigma.diag()));

      vec aggr_stats = sum(marker_expr_imputed, 1); // each column is a marker gene
      aggr_stats = aggr_stats / norm_factor;
      marker_stats.col(j) = aggr_stats;
    }
    mat marker_stats_smoothed = compute_network_diffusion_Chebyshev(P, marker_stats, thread_no);

    return(marker_stats_smoothed);
  }



} // namespace ACTIONet