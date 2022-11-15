#include <ACTIONet.h>

namespace ACTIONet
{

  mat assess_label_enrichment(sp_mat &H, mat &M, int thread_no)
  {
    mat Obs = spmat_mat_product_parallel(H, M, thread_no);

    rowvec p = mean(M, 0);
    mat Exp = sum(H, 1) * p;

    mat Lambda = Obs - Exp;

    mat Nu = (sum(square(H), 1) * p);
    vec a = vec(max(H, 1));

    mat Lambda_scaled = Lambda;
    for (int j = 0; j < Lambda_scaled.n_rows; j++)
    {
      Lambda_scaled.row(j) *= (a(j) / 3);
    }

    mat logPvals_upper = square(Lambda) / (2 * (Nu + Lambda_scaled));
    uvec lidx = find(Lambda <= 0);
    logPvals_upper(lidx) = zeros(lidx.n_elem);
    logPvals_upper.replace(datum::nan, 0); // replace each NaN with 0

    return logPvals_upper;
  }

  mat one_hot_encoding(vec labels)
  {
    int n = labels.n_elem;

    vec vals = unique(labels);

    uvec idx = find(0 <= vals);
    vals = vals(idx);

    int k = vals.n_elem;
    mat M = zeros(n, k);
    for (int i = 0; i < k; i++)
    {
      uvec idx = find(labels == vals(i));
      for (int j = 0; j < idx.n_elem; j++)
      {
        M(idx(j), i) = 1;
      }
    }

    return (M);
  }

  vec LPA(sp_mat &G, vec labels, double lambda, int iters,
          double sig_threshold, uvec fixed_labels, int thread_no)
  {
    int n = G.n_rows;

    vec updated_labels = labels;

    sp_mat H = G;
    H.diag().ones();
    H.diag() *= lambda;          // add "inertia"
    H = n * normalize_adj(H, 1); // row-normalize to n

    int it = 0;
    for (it; it < iters; it++)
    {
      vec vals = unique(updated_labels);
      uvec idx = find(0 <= vals);
      vals = vals(idx);

      mat M = one_hot_encoding(updated_labels);

      mat logPvals = assess_label_enrichment(H, M, thread_no);

      vec max_sig = max(logPvals, 1);
      vec new_labels = vals(index_max(logPvals, 1));

      // Only update vertices with significant enrichment in their neighborhood
      uvec sig_idx = find(sig_threshold < max_sig);
      updated_labels(sig_idx) = new_labels(sig_idx);

      // revert-back the vertices that need to be fixed
      updated_labels(fixed_labels) = labels(fixed_labels);
    }
    stdout_printf("\r\tLPA iteration %d/%d", it, iters);
    FLUSH;

    return (updated_labels);
  }

} // namespace ACTIONet
