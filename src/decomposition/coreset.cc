#include "ACTIONet.h"

namespace ACTIONet {

Coreset compute_AA_coreset(sp_mat &S, int m = 5000) {
  vec mu = vec(mean(S, 1));

  vec q = zeros(S.n_cols);
  sp_mat::const_iterator it = S.begin();
  sp_mat::const_iterator it_end = S.end();
  for (; it != it_end; ++it) {
    double d = (*it) - mu[it.row()];
    q[it.col()] += d * d;
  }
  vec p = q / sum(q);  // sampling probability

  if (m == 0) {
    double p_sum = sum(p);
    double p_sq_sum = sum(square(p));
    m = (int)((p_sum * p_sum) / (p_sq_sum));
  }
  m = std::min(m, (int)S.n_cols);

  vec p_sorted = sort(p, "descend");
  double p_threshold = p_sorted(m - 1);

  uvec sample_idx = find(p_threshold <= p);
  m = sample_idx.n_elem;  // Can do actual sampling, maybe!

  mat S_coreset(S.n_rows, m);
  for (int j = 0; j < m; j++) {
    S_coreset.col(j) = vec(S.col(sample_idx(j)));
  }
  vec w_coreset = 1.0 / (m * p(sample_idx));

  Coreset coreset;
  coreset.S_coreset = S_coreset;
  coreset.w_coreset = w_coreset;
  coreset.index = sample_idx;

  return (coreset);
}
}  // namespace ACTIONet
