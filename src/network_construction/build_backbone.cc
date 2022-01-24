#include <ACTIONet.h>
#include <hnswlib.h>
#include <atomic>
#include <thread>
namespace ACTIONet {

double threshold_vector(vec &dist, double LC) {
  vec dist_sorted = sort(dist, "ascend");
  vec beta = join_vert(LC * dist_sorted, datum::inf * ones(1));
  vec beta_sq = square(beta);

  double lambda = beta(0) + 1, k = 0, Sum_beta = 0, Sum_beta_square = 0;

  for (; k < beta.n_elem - 1; k++) {
    Sum_beta += beta[k];
    Sum_beta_square += (beta_sq[k]);
    lambda = (1.0 / k) *
             (Sum_beta + sqrt(k + Sum_beta * Sum_beta - k * Sum_beta_square));

    if (lambda <= beta(k + 1)) break;
  }

  double dist_threshold = dist_sorted[(int)k];

  return (dist_threshold);
}

// k^{*}-Nearest Neighbors: From Global to Local (NIPS 2016)
mat Prune_PageRank(mat &U, double density = 1.0) {
  double LC = 1.0 / density;

  mat Dist = -log(U);

  mat U_thresholded(size(U));
  for (int i = 0; i < U.n_cols; i++) {
    vec d = Dist.col(i);
    double dist_threshold = threshold_vector(d, LC);

    uvec filter_idx = find(dist_threshold <= Dist.col(i));
    vec u = U.col(i);
    u(filter_idx).zeros();

    U_thresholded.col(i) = u;
  }

  U_thresholded = normalise(U_thresholded, 1, 0);

  return (U_thresholded);
}
}  // namespace ACTIONet
