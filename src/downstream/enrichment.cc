#include <ACTIONet.h>

namespace ACTIONet {

field<mat> assess_enrichment(mat &scores, sp_mat &associations,
                             int thread_no = 0) {
  field<mat> res(3);

  if (scores.n_rows != associations.n_rows) {
    fprintf(stderr,
            "Number of rows in scores and association matrices should both "
            "match the number of features\n");
  }

  associations = spones(associations);

  mat sorted_scores = sort(scores, "descend");
  vec a_max = trans(sorted_scores.row(0));
  umat perms(size(scores));
  for (int j = 0; j < scores.n_cols; j++) {
    perms.col(j) =
        stable_sort_index(stable_sort_index(scores.col(j), "descend"));
  }

  vec n_success = vec(trans(sum(associations, 0)));
  vec p_success = n_success / (double)associations.n_rows;

  mat Acumsum = cumsum(sorted_scores);
  mat A2cumsum = cumsum(square(sorted_scores));

  mat logPvals = zeros(associations.n_cols, scores.n_cols);
  mat thresholds = zeros(associations.n_cols, scores.n_cols);

  // for(int k = 0; k < associations.n_cols; k++) {
  parallelFor(0, associations.n_cols, [&] (size_t k) {
        int n_k = n_success(k);
        if (n_k > 1) {
          double p_k = p_success(k);

          mat O = zeros(n_k, scores.n_cols);
          mat E = zeros(n_k, scores.n_cols);
          mat Nu = zeros(n_k, scores.n_cols);
          mat rows = zeros(n_k, scores.n_cols);

          for (int j = 0; j < scores.n_cols; j++) {
            uvec perm = perms.col(j);

            uvec sorted_rows(n_k);
            sp_mat::const_col_iterator it = associations.begin_col(k);
            sp_mat::const_col_iterator it_end = associations.end_col(k);
            for (int idx = 0; it != it_end; ++it, idx++) {
              sorted_rows[idx] = perm[it.row()];
            }
            sorted_rows = sort(sorted_rows);

            for (int idx = 0; idx < n_k; idx++) {
              int ii = sorted_rows(idx);

              O(idx, j) = sorted_scores(ii, j);
              E(idx, j) = Acumsum(ii, j) * p_k;
              Nu(idx, j) = A2cumsum(ii, j) * p_k;
              rows(idx, j) = ii;
            }
          }
          O = cumsum(O);

          mat Lambda = O - E;
          mat aLambda = Lambda;
          for (int j = 0; j < aLambda.n_cols; j++) {
            aLambda.col(j) *= a_max(j);
          }

          mat logPvals_k = square(Lambda) / (2.0 * (Nu + (aLambda / 3.0)));
          uvec idx = find(Lambda <= 0);
          logPvals_k(idx) = zeros(idx.n_elem);
          logPvals_k.replace(datum::nan, 0);
          for (int j = 0; j < logPvals_k.n_cols; j++) {
            vec v = logPvals_k.col(j);
            logPvals(k, j) = max(v);
            thresholds(k, j) = rows[v.index_max(), j];
          }
        }
      }, thread_no);

  field<mat> output(2);
  output(0) = logPvals;
  output(1) = thresholds;

  return (output);
}

}  // namespace ACTIONet
