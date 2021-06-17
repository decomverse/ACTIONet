#include <ACTIONet.h>

namespace ACTIONet {
mat compute_pseudo_bulk_per_cluster(
    sp_mat &S, arma::Col<unsigned long long> sample_assignments) {
  mat pb = zeros(S.n_rows, max(sample_assignments));

  sp_mat::const_iterator it = S.begin();
  sp_mat::const_iterator it_end = S.end();
  for (; it != it_end; ++it) {
    int i = it.row();
    int j = sample_assignments[it.col()] - 1;
    pb(i, j) += (*it);
  }

  for (int j = 0; j < pb.n_cols; j++) {
    uvec idx = find(sample_assignments == (j + 1));
    pb.col(j) /= max(1, idx.n_elem);
  }

  return (pb);
}

mat compute_pseudo_bulk_per_cluster(
    mat &S, arma::Col<unsigned long long> sample_assignments) {
  mat pb = zeros(S.n_rows, max(sample_assignments));

  for (int j = 0; j < pb.n_cols; j++) {
    uvec idx = find(sample_assignments == (j + 1));
    if (idx.n_elem == 0) continue;

    if (idx.n_elem > 1) {
      mat subS = S.cols(idx);
      pb.col(j) = mean(subS, 1);
    } else {
      pb.col(j) = S.col(idx(0));
    }
  }

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

field<mat> compute_pseudo_bulk_per_cluster_and_ind(
    sp_mat &S, arma::Col<unsigned long long> sample_assignments,
    arma::Col<unsigned long long> individuals) {
  field<mat> pbs(max(sample_assignments));
  for (int k = 0; k < max(sample_assignments); k++) {
    pbs(k) = zeros(S.n_rows, max(individuals));
  }

  sp_mat::const_iterator it = S.begin();
  sp_mat::const_iterator it_end = S.end();
  for (; it != it_end; ++it) {
    int i = it.row();
    int j = individuals[it.col()] - 1;
    int k = sample_assignments[it.col()] - 1;

    pbs(k)(i, j) += (*it);
  }

  for (int j = 0; j < max(individuals); j++) {
    for (int k = 0; k < max(sample_assignments); k++) {
      uvec idx = intersect(find((sample_assignments == (k + 1))),
                           find((individuals == (j + 1))));

      pbs(k).col(j) /= max(1, idx.n_elem);
    }
  }

  return (pbs);
}

field<mat> compute_pseudo_bulk_per_cluster_and_ind(
    mat &S, arma::Col<unsigned long long> sample_assignments,
    arma::Col<unsigned long long> individuals) {
  field<mat> pbs(max(sample_assignments));
  for (int k = 0; k < max(sample_assignments); k++) {
    pbs(k) = zeros(S.n_rows, max(individuals));
  }

  for (int j = 0; j < max(individuals); j++) {
    for (int k = 0; k < max(sample_assignments); k++) {
      uvec idx = intersect(find((sample_assignments == (k + 1))),
                           find((individuals == (j + 1))));
      mat subS = S.cols(idx);
      pbs(k).col(j) = mean(subS, 1);
    }
  }

  return (pbs);
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
