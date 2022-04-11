#include <ACTIONet.h>

namespace ACTIONet {

mat compute_grouped_rowsums(sp_mat &S, arma::Col<unsigned long long> sample_assignments) {

  uvec lv_vec = conv_to<uvec>::from(unique(sample_assignments));
  mat pb = zeros(S.n_rows, lv_vec.n_elem);

  sp_mat::const_iterator it = S.begin();
  sp_mat::const_iterator it_end = S.end();
  for (; it != it_end; ++it) {
    int i = it.row();
    int j = sample_assignments[it.col()] - 1;
    pb(i, j) += (*it);
  }

  return (pb);
}

mat compute_grouped_rowsums(mat &S, arma::Col<unsigned long long> sample_assignments) {

  uvec lv_vec = conv_to<uvec>::from(unique(sample_assignments));
  mat pb = zeros(S.n_rows, lv_vec.n_elem);
  
  for (int j = 0; j < pb.n_cols; j++) {
    uvec idx = find(sample_assignments == (j + 1));
    if (idx.n_elem == 0) continue;

    if (idx.n_elem > 1) {
      mat subS = S.cols(idx);
      pb.col(j) = sum(subS, 1);
    } else {
      pb.col(j) = S.col(idx(0));
    }
  }

  return (pb);
}

mat compute_grouped_rowmeans(sp_mat &S, arma::Col<unsigned long long> sample_assignments) {

  mat pb = compute_grouped_rowsums(S, sample_assignments);

  for (int j = 0; j < pb.n_cols; j++) {
    uvec idx = find(sample_assignments == (j + 1));
    pb.col(j) /= max(1, (int)idx.n_elem);
  }

  return (pb);
}

mat compute_grouped_rowmeans(mat &S, arma::Col<unsigned long long> sample_assignments) {
  
  mat pb = compute_grouped_rowsums(S, sample_assignments);

  for (int j = 0; j < pb.n_cols; j++) {
    uvec idx = find(sample_assignments == (j + 1));
    pb.col(j) /= max(1, (int)idx.n_elem);
  }

  return (pb);
}

mat compute_grouped_rowvars(sp_mat &S, arma::Col<unsigned long long> sample_assignments) {

  uvec lv_vec = conv_to<uvec>::from(unique(sample_assignments));
  mat pb = zeros(S.n_rows, lv_vec.n_elem);

  // mat pb_mean = compute_grouped_rowmeans(S, sample_assignments);

  for (int i = 0; i < S.n_rows; i++) {
    // rowvec r = arma::conv_to<rowvec>::from(S.row(i));
    sp_mat::const_row_iterator it = S.begin_row(i);
    sp_mat::const_row_iterator it_end = S.end_row(i);
    vec r = zeros(S.n_cols);

    for (; it != it_end; ++it) {
      // int k = sample_assignments[it.col()] - 1;
      r(it.col()) = (*it);
    }

    for (int j = 0; j < pb.n_cols; j++){
      uvec idx = find(sample_assignments == (j + 1));
      pb(i, j) = var(r(idx));
    }
  }

  // sp_mat::const_iterator it = S.begin();
  // sp_mat::const_iterator it_end = S.end();
  // for (; it != it_end; ++it) {
  //   int i = it.row();
  //   int j = sample_assignments[it.col()] - 1;
  //   pb(i, j) += (*it);
  // }

  return (pb);
}

mat compute_pseudo_bulk_per_archetype(sp_mat &S, mat &H) {
  mat H_norm = trans(H);

  vec col_means = trans(mean(H, 0));
  for (int i = 0; i < H_norm.n_cols; i++) {
		double denom = col_means(i);
		H_norm.row(i) /= denom == 0?1:denom;
  }

  mat pb = mat(S * H_norm);

  return (pb);
}

mat compute_pseudo_bulk_per_archetype(mat &S, mat &H) {
  mat H_norm = trans(H);

  vec col_means = trans(mean(H, 0));
  for (int i = 0; i < H_norm.n_cols; i++) {
	  double denom = col_means(i);
    H_norm.row(i) /= denom == 0?1:denom;
  }

  mat pb = mat(S * H_norm);

  return (pb);
}


// H: archs x cells
// individuals: 0 ... (ind-1)
field<mat> compute_pseudo_bulk_per_archetype_and_ind(
    sp_mat &S, arma::mat &H,
    arma::Col<unsigned long long> individuals) {
  mat H_norm = trans(H); // cell x archs
  vec col_means = trans(mean(H, 0));
  for (int i = 0; i < H_norm.n_cols; i++) {
	  double denom = col_means(i);
    H_norm.row(i) /= denom == 0?1:denom;
  }

  int arch_no = H_norm.n_cols;				
  field<mat> pbs(arch_no);
  for (int k = 0; k < arch_no; k++) {
    pbs(k) = zeros(S.n_rows, max(individuals));
  }

  sp_mat::const_iterator it = S.begin();
  sp_mat::const_iterator it_end = S.end();
  for (; it != it_end; ++it) {
    int i = it.row();
    int j = individuals[it.col()] - 1;

	for (int k = 0; k < arch_no; k++) {
		pbs(k)(i, j) += (H_norm(j, k) * (*it));
	}
  }

  return (pbs);
}

field<mat> compute_pseudo_bulk_per_archetype_and_ind(
    mat &S, mat &H,
    arma::Col<unsigned long long> individuals) {
  mat H_norm = trans(H); // cell x archs
  vec col_means = trans(mean(H, 0));
  for (int i = 0; i < H_norm.n_cols; i++) {
	  double denom = col_means(i);
    H_norm.row(i) /= denom == 0?1:denom;
  }

  int arch_no = H_norm.n_cols;		
  int ind_no = max(individuals);		
  field<mat> pbs(arch_no);
  for (int k = 0; k < arch_no; k++) {
    pbs(k) = zeros(S.n_rows, ind_no);
  }    

  for (int j = 0; j < ind_no; j++) {
    uvec idx = find((individuals == (j + 1)));
    mat subS = S.cols(idx);
    mat subH = H_norm.rows(idx);
    
    mat X = (subS * subH);
	for (int k = 0; k < arch_no; k++) {
		pbs(k).col(j) += X.col(k);
	}
  }

  return (pbs);
}


}  // namespace ACTIONet
