#include <ACTIONet.h>

namespace ACTIONet
{

  mat normalize_scores(mat scores, int method, int thread_no)
  {
    mat normalized_scores(size(scores));
    switch (method)
    {
    case 0: //"none"
    {
      normalized_scores = scores;
      break;
    }
    case 1: //"zscore"
    {
      normalized_scores = zscore(scores, 0, thread_no);
      break;
    }
    case 2: //"RINT" (nonparametric)
    {
      normalized_scores = RIN_transform(scores, thread_no);
      break;
    }
    case 3: //"robust_zscore" (kinda hack!)
    {
      normalized_scores = robust_zscore(scores, 0, thread_no);
      break;
    }
    default:
      stderr_printf("Unknown normalization method\n");
      FLUSH;
      normalized_scores = scores;
    }
    return (normalized_scores);
  }

  // G is the cell-cell network, scores is a cell x geneset matrix
  field<vec> autocorrelation_Moran(mat G, mat scores, int normalization_method,
                                   int perm_no, int thread_no)
  {
    int nV = G.n_rows;
    int scores_no = scores.n_cols;
    stdout_printf("Normalizing scores (method=%d) ... ", normalization_method);
    mat normalized_scores =
        normalize_scores(scores, normalization_method, thread_no);
    stdout_printf("done\n");
    FLUSH;

    stdout_printf("Computing auto-correlation over network ... ", nV, scores_no);

    double W = sum(sum(G));
    vec norm_sq = vec(trans(sum(square(normalized_scores))));
    vec norm_factors = nV / (W * norm_sq);
    norm_factors.replace(datum::nan, 0); // replace each NaN with 0

    vec stat = zeros(scores_no);
    parallelFor(
        0, scores_no,
        [&](unsigned int i)
        {
          vec x = normalized_scores.col(i);
          double y = dot(x, G * x);
          stat(i) = y;
        },
        thread_no);

    vec mu = zeros(scores_no);
    vec sigma = zeros(scores_no);
    vec z = zeros(scores_no);
    if (0 < perm_no)
    {
      mat rand_stats = zeros(scores_no, perm_no);
      parallelFor(
          0, perm_no,
          [&](unsigned int j)
          {
            uvec perm = randperm(nV);
            mat score_permuted = normalized_scores.rows(perm);

            vec v = zeros(scores_no);
            for (int i = 0; i < scores_no; i++)
            {
              vec rand_x = score_permuted.col(i);
              v(i) = dot(rand_x, G * rand_x);
            }
            rand_stats.col(j) = v;
          },
          thread_no);

      mu = mean(rand_stats, 1);
      sigma = stddev(rand_stats, 0, 1);
      z = (stat - mu) / sigma;
      z.replace(datum::nan, 0);
    }
    stdout_printf("done\n");
    FLUSH;

    // Summary stats
    field<vec> results(4);
    results(0) = stat % norm_factors;
    results(1) = z;
    results(2) = mu;
    results(3) = sigma;
    return (results);
  }

  // G is the cell-cell network, scores is a cell x geneset matrix
  field<vec> autocorrelation_Moran(sp_mat G, mat scores, int normalization_method,
                                   int perm_no, int thread_no)
  {
    int nV = G.n_rows;
    int scores_no = scores.n_cols;

    mat normalized_scores =
        normalize_scores(scores, normalization_method, thread_no);

    stdout_printf("Computing auto-correlation over network ... ", nV, scores_no);
    double W = sum(sum(G));
    vec norm_sq = vec(trans(sum(square(normalized_scores))));
    vec norm_factors = nV / (W * norm_sq);
    norm_factors.replace(datum::nan, 0); // replace each NaN with 0

    vec stat = zeros(scores_no);
    parallelFor(
        0, scores_no,
        [&](unsigned int i)
        {
          vec x = normalized_scores.col(i);
          double y = dot(x, spmat_vec_product(G, x));
          stat(i) = y;
        },
        thread_no);

    vec mu = zeros(scores_no);
    vec sigma = zeros(scores_no);
    vec z = zeros(scores_no);
    if (0 < perm_no)
    {
      mat rand_stats = zeros(scores_no, perm_no);
      parallelFor(
          0, perm_no,
          [&](unsigned int j)
          {
            uvec perm = randperm(nV);
            mat score_permuted = normalized_scores.rows(perm);

            vec v = zeros(scores_no);
            for (int i = 0; i < scores_no; i++)
            {
              vec rand_x = score_permuted.col(i);
              v(i) = dot(rand_x, spmat_vec_product(G, rand_x));
            }
            rand_stats.col(j) = v;
          },
          thread_no);

      mu = mean(rand_stats, 1);
      sigma = stddev(rand_stats, 0, 1);
      z = (stat - mu) / sigma;
      z.replace(datum::nan, 0);
    }
    stdout_printf("done\n");
    FLUSH;

    // Summary stats
    field<vec> results(4);
    results(0) = stat % norm_factors;
    results(1) = z;
    results(2) = mu;
    results(3) = sigma;

    return (results);
  }

  // G is the cell-cell network, scores is a cell x geneset matrix
  field<vec> autocorrelation_Geary(mat G, mat scores, int normalization_method,
                                   int perm_no, int thread_no)
  {
    int nV = G.n_rows;
    int scores_no = scores.n_cols;
    mat normalized_scores =
        normalize_scores(scores, normalization_method, thread_no);

    stdout_printf("Computing auto-correlation over network ... ", nV, scores_no);
    double W = sum(sum(G));
    vec norm_sq = vec(trans(sum(square(normalized_scores))));
    vec norm_factors = (nV - 1) / ((2 * W) * norm_sq);
    norm_factors.replace(datum::nan, 0); // replace each NaN with 0

    // Compute graph Laplacian
    vec d = vec(trans(sum(G)));
    mat L(-G);
    L.diag() = d;

    vec stat = zeros(scores_no);
    parallelFor(
        0, scores_no,
        [&](unsigned int i)
        {
          vec x = normalized_scores.col(i);
          double y = dot(x, L * x);
          stat(i) = y;
        },
        thread_no);

    vec mu = zeros(scores_no);
    vec sigma = zeros(scores_no);
    vec z = zeros(scores_no);
    if (0 < perm_no)
    {
      mat rand_stats = zeros(scores_no, perm_no);
      parallelFor(
          0, perm_no,
          [&](unsigned int j)
          {
            uvec perm = randperm(nV);
            mat score_permuted = normalized_scores.rows(perm);

            vec v = zeros(scores_no);
            for (int i = 0; i < scores_no; i++)
            {
              vec rand_x = score_permuted.col(i);
              v(i) = dot(rand_x, L * rand_x);
            }
            rand_stats.col(j) = v;
          },
          thread_no);

      mu = mean(rand_stats, 1);
      sigma = stddev(rand_stats, 0, 1);
      z = (stat - mu) / sigma;
      z.replace(datum::nan, 0);
    }
    stdout_printf("done\n");
    FLUSH;

    // Summary stats
    field<vec> results(4);
    results(0) = stat % norm_factors;
    results(1) = -z;
    results(2) = mu;
    results(3) = sigma;
    return (results);
  }

  // G is the cell-cell network, scores is a cell x geneset matrix
  field<vec> autocorrelation_Geary(sp_mat G, mat scores, int normalization_method,
                                   int perm_no, int thread_no)
  {
    int nV = G.n_rows;
    int scores_no = scores.n_cols;
    mat normalized_scores =
        normalize_scores(scores, normalization_method, thread_no);

    stdout_printf("Computing auto-correlation over network ... ", nV, scores_no);
    double W = sum(sum(G));
    vec norm_sq = vec(trans(sum(square(normalized_scores))));
    vec norm_factors = (nV - 1) / ((2 * W) * norm_sq);
    norm_factors.replace(datum::nan, 0); // replace each NaN with 0

    // Compute graph Laplacian
    vec d = vec(trans(sum(G)));
    sp_mat L(-G);
    L.diag() = d;

    vec stat = zeros(scores_no);
    parallelFor(
        0, scores_no,
        [&](unsigned int i)
        {
          vec x = normalized_scores.col(i);
          double y = dot(x, spmat_vec_product(L, x));
          stat(i) = y;
        },
        thread_no);

    vec mu = zeros(scores_no);
    vec sigma = zeros(scores_no);
    vec z = zeros(scores_no);
    if (0 < perm_no)
    {
      mat rand_stats = zeros(scores_no, perm_no);
      parallelFor(
          0, perm_no,
          [&](unsigned int j)
          {
            uvec perm = randperm(nV);
            mat score_permuted = normalized_scores.rows(perm);

            vec v = zeros(scores_no);
            for (int i = 0; i < scores_no; i++)
            {
              vec rand_x = score_permuted.col(i);
              v(i) = dot(rand_x, spmat_vec_product(L, rand_x));
            }
            rand_stats.col(j) = v;
          },
          thread_no);

      mu = mean(rand_stats, 1);
      sigma = stddev(rand_stats, 0, 1);
      z = (stat - mu) / sigma;
      z.replace(datum::nan, 0);
    }
    stdout_printf("done\n");
    FLUSH;

    // Summary stats
    field<vec> results(4);
    results(0) = stat % norm_factors;
    results(1) = -z;
    results(2) = mu;
    results(3) = sigma;

    return (results);
  }

} // namespace ACTIONet
