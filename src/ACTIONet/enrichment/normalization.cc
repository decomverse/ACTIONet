#include <ACTIONet.h>

namespace ACTIONet {
sp_mat
renormalize_input_matrix(sp_mat &S,
                         arma::Col<unsigned long long> sample_assignments) {
  printf("Computing and normalizing pseudobulk profiles ... ");
  mat pb = compute_pseudo_bulk_per_cluster(S, sample_assignments);
  rowvec mean_pb_t = trans(mean(pb, 1));

  // Aligns pseudo-bulk profiles with each other
  for (int j = 0; j < pb.n_cols; j++) {
    vec x = pb.col(j);

    double num = sum(mean_pb_t * x);
    double denom = sum(square(x));

    pb.col(j) *= (num / denom);
  }
  printf("done\n");

  // Align individual columns now
  printf("Normalizing individual columns ... ");

  sp_mat S_norm = S;
  mat pb_t = trans(pb);
  for (int j = 0; j < S.n_cols; j++) {
    vec x = vec(S.col(j));
    rowvec y_t = pb_t.row(sample_assignments(j) - 1);

    double num = sum(y_t * x);
    double denom = sum(square(x));

    S_norm.col(j) *= (num / denom);
  }
  printf("done\n");

  return (S_norm);
}

mat renormalize_input_matrix(mat &S,
                             arma::Col<unsigned long long> sample_assignments) {
  printf("Computing and normalizing pseudobulk profiles ... ");
  mat pb = compute_pseudo_bulk_per_cluster(S, sample_assignments);
  rowvec mean_pb_t = trans(mean(pb, 1));

  // Aligns pseudo-bulk profiles with each other
  for (int j = 0; j < pb.n_cols; j++) {
    vec x = pb.col(j);

    double num = sum(mean_pb_t * x);
    double denom = sum(square(x));

    pb.col(j) *= (num / denom);
  }
  printf("done\n");

  // Align individual columns now
  printf("Normalizing individual columns ... ");

  mat S_norm = S;
  mat pb_t = trans(pb);
  for (int j = 0; j < S.n_cols; j++) {
    vec x = vec(S.col(j));
    rowvec y_t = pb_t.row(sample_assignments(j) - 1);

    double num = sum(y_t * x);
    double denom = sum(square(x));

    S_norm.col(j) *= (num / denom);
  }
  printf("done\n");

  return (S_norm);
}

} // namespace ACTIONet
