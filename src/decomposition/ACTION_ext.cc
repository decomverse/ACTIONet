#include "ACTIONet.h"

namespace ACTIONet {
field<mat> run_AA(mat &A, mat &W0, int max_it, double min_delta);
// Solves the weighted Archetypal Analysis (AA) problem
field<mat> run_weighted_AA(mat &A, mat &W0, vec w, int max_it = 50,
                           double min_delta = 0.01) {
  int N = A.n_cols;
  field<mat> decomposition(2);

  if (N != w.n_elem) {
    stdout_printf("Number of elements in the weight vector should match the total number of samples (columns in A)\n");
    FLUSH;
    return (decomposition);
  }

  w = clamp(w, 0, 1);
  mat A_scaled = A;
  for (int i = 0; i < N; i++) {
    A_scaled.col(i) *= w[i];
  }
  decomposition = run_AA(A_scaled, W0, max_it, min_delta);

  mat C = decomposition(0);
  mat weighted_archs = A_scaled * C;
  mat H = run_simplex_regression(weighted_archs, A, false);
  decomposition(1) = H;

  return (decomposition);
}

ACTION_results run_weighted_ACTION(mat &S_r, vec w, int k_min, int k_max,
                                   int thread_no, int max_it = 50,
                                   double min_delta = 1e-16) {
  int feature_no = S_r.n_rows;

  stdout_printf("Running weighted ACTION (%d threads):", thread_no); FLUSH;

  if (k_max == -1) k_max = (int)S_r.n_cols;

  k_min = std::max(k_min, 2);
  k_max = std::min(k_max, (int)S_r.n_cols);

  ACTION_results trace;

  trace.H = field<mat>(k_max + 1);
  trace.C = field<mat>(k_max + 1);
  trace.selected_cols = field<uvec>(k_max + 1);

  int N = S_r.n_cols;
  if (N != w.n_elem) {
    stderr_printf("Number of elements in the weight vector should match the total number of samples (columns in S_r)\n");
    FLUSH;
    return (trace);
  }

  w = clamp(w, 0, 1);
  mat X_r = normalise(S_r, 1);  // ATTENTION!

  mat X_r_scaled = X_r;
  for (int i = 0; i < N; i++) {
    X_r_scaled.col(i) *= w[i];
  }

  int current_k = 0;
  char status_msg[50];

  sprintf(status_msg, "Iterating from k = %d ... %d:", k_min, k_max);
  stderr_printf("\n\t%s %d/%d finished", status_msg, current_k, (k_max - k_min + 1));
  FLUSH;

  parallelFor(k_min, k_max + 1, [&] (size_t kk) {
    SPA_results SPA_res = run_SPA(X_r_scaled, kk);
    trace.selected_cols[kk] = SPA_res.selected_columns;

    mat W = X_r_scaled.cols(trace.selected_cols[kk]);

    field<mat> AA_res = run_AA(X_r_scaled, W, max_it, min_delta);

    mat C = AA_res(0);
    mat weighted_archs = X_r_scaled * C;
    mat H = run_simplex_regression(weighted_archs, X_r, false);
    AA_res(1) = H;

    trace.C[kk] = AA_res(0);
    trace.H[kk] = AA_res(1);
    current_k++;

    stderr_printf("\r\t%s %d/%d finished", status_msg, current_k,
                  (k_max - k_min + 1));
    FLUSH;

  }, thread_no);

  stdout_printf("\r\t%s %d/%d finished\n", status_msg, current_k,
                (k_max - k_min + 1)); FLUSH;

  return trace;
}


field<mat> Online_update_AA(mat &Xt, mat &D, mat &A, mat &B) {
  // Compute archetype coefficients using the last learned dictionary
  mat Ct = run_simplex_regression(D, Xt, false);

  // Just in case!
  Ct = clamp(Ct, 0, 1);
  Ct = normalise(Ct, 1);
  mat Ct_T = trans(Ct);

  // Update sufficient statistics
  mat At = A + Ct * Ct_T;
  mat Bt = B + Xt * Ct_T;

  // Update the dictionary using block-coordinate-descent (BCD)
  mat Dt(size(D));
  for (int j = 0; j < D.n_cols; j++) {
    vec u = D.col(j) + (1.0 / At(j, j)) * (Bt.col(j) - D * At.col(j));
    Dt.col(j) = u / std::max(norm(u, 2), 1.0);
  }

  field<mat> decomposition(4, 1);

  decomposition(0) = At;
  decomposition(1) = Bt;
  decomposition(2) = Ct;
  decomposition(3) = Dt;

  return (decomposition);
}

field<mat> run_online_AA(mat &X, mat &D0, field<uvec> samples) {
  int m = X.n_rows;
  int n = X.n_cols;
  int k = D0.n_cols;

  mat At = zeros(k, k);
  mat Bt = zeros(m, k);
  mat Ct(k, n), Ct_T(n, k);

  // Just in case
  X = normalise(X, 2, 0);
  mat Dt = normalise(D0, 2, 0);

  for (int t = 0; t < samples.n_elem; t++) {
    // Extract the next batch
    uvec idx = samples(t);
    mat Xt = X.cols(idx);

    // Compute archetype coefficients using the last learned dictionary
    Ct = run_simplex_regression(Dt, Xt, false);
    Ct_T = trans(Ct);

    // Update sufficient statistics
    At += Ct * Ct_T;
    Bt += Xt * Ct_T;

    // Update the dictionary using block-coordinate-descent (BCD)
    for (int j = 0; j < k; j++) {
      vec u = Dt.col(j) + (1.0 / At(j, j)) * (Bt.col(j) - Dt * At.col(j));
      Dt.col(j) = u / std::max(norm(u, 2), 1.0);
    }
  }

  // Just in case!
  Ct = clamp(Ct, 0, 1);
  Ct = normalise(Ct, 1);

  field<mat> decomposition(4, 1);

  decomposition(0) = At;
  decomposition(1) = Bt;
  decomposition(2) = Ct;
  decomposition(3) = Dt;

  return (decomposition);
}

Online_ACTION_results run_online_ACTION(mat &S_r, field<uvec> samples,
                                        int k_min, int k_max, int thread_no) {
  int feature_no = S_r.n_rows;

  stdout_printf("Running online ACTION (%d threads):", thread_no); FLUSH;

  if (k_max == -1) k_max = (int)S_r.n_cols;

  k_min = std::max(k_min, 2);
  k_max = std::min(k_max, (int)S_r.n_cols);

  Online_ACTION_results trace;

  trace.A = field<mat>(k_max + 1);
  trace.B = field<mat>(k_max + 1);
  trace.C = field<mat>(k_max + 1);
  trace.D = field<mat>(k_max + 1);
  trace.selected_cols = field<uvec>(k_max + 1);

  mat X_r_L1 = normalise(S_r, 1, 0);
  mat X_r_L2 = normalise(S_r, 2, 0);

  int current_k = 0;
  char status_msg[50];

  sprintf(status_msg, "Iterating from k = %d ... %d:", k_min, k_max);
  stderr_printf("\n\t%s %d/%d finished", status_msg, current_k,
                (k_max - k_min + 1));
  FLUSH;

  parallelFor(k_min, k_max + 1, [&] (size_t kk) {
    SPA_results SPA_res = run_SPA(X_r_L1, kk);
    trace.selected_cols[kk] = SPA_res.selected_columns;

    mat W = X_r_L2.cols(trace.selected_cols[kk]);

    field<mat> AA_res;
    AA_res = run_online_AA(X_r_L2, W, samples);

    trace.A[kk] = AA_res(0);
    trace.B[kk] = AA_res(1);
    trace.C[kk] = AA_res(2);
    trace.D[kk] = AA_res(3);

    current_k++;

    stderr_printf("\r\t%s %d/%d finished", status_msg, current_k,
                  (k_max - k_min + 1));
    FLUSH;
  }, thread_no);

  stdout_printf("\r\t%s %d/%d finished\n", status_msg, current_k,
                (k_max - k_min + 1)); FLUSH;

  return trace;
}


	mat oneHot_encoding(vec batches) {
		vec uniue_batches = sort(unique(batches));

		mat encoding = zeros(uniue_batches.n_elem, batches.n_elem);
		for(int i = 0; i < uniue_batches.n_elem; i++) {
			uvec idx = find(batches == uniue_batches[i]);
			vec batch_encoding = zeros(batches.n_elem);
			batch_encoding.elem(idx).ones();

			encoding.row(i) = trans(batch_encoding);
		}

		return(encoding);
	}

field<mat> run_AA_with_batch_correction(mat &Z, mat &W0, vec batch, int max_it = 100, int max_correction_rounds = 10, double lambda = 1, double min_delta = 1e-6) {
	mat Phi = oneHot_encoding(batch);
	mat Phi_moe = join_vert(ones(1, Z.n_cols), Phi);

int sample_no = Z.n_cols, k = W0.n_cols;

  mat C = zeros(sample_no, k);
  mat H = zeros(k, sample_no);

	// First round is just AA using raw input
	mat W = W0;
	mat Z_corr  = Z;

	for(int correction_round = 0; correction_round < max_correction_rounds; correction_round++) {
		field<mat> AA_res = run_AA(Z_corr, W, max_it, min_delta);
		C = AA_res(0);
		H = AA_res(1);

		// Correction using mixture of experts -- Adopted from the Harmony method
		Z_corr = Z;
		for(int k = 0; k < H.n_rows; k++) {
			rowvec h = H.row(k);
			//mat Phi_Rk = Phi_moe * arma::diagmat(h);
			mat Phi_Rk = Phi_moe.each_row() % h;

			mat beta = arma::inv(Phi_Rk * Phi_moe.t() + lambda) * Phi_Rk * Z.t();
			beta.row(0).zeros(); // do not remove the intercept
			Z_corr -= beta.t() * Phi_Rk;
		}
		W = Z_corr * C; // Start the next round of AA from current state
	}


  C = clamp(C, 0, 1);
  C = normalise(C, 1);
  H = clamp(H, 0, 1);
  H = normalise(H, 1);

  field<mat> decomposition(3);
  decomposition(0) = C;
  decomposition(1) = H;
  decomposition(2) = Z_corr;


  return decomposition;
}


ACTION_results run_ACTION_with_batch_correction(mat &S_r, vec batch, int k_min, int k_max, int thread_no,
                          int max_it = 100, int max_correction_rounds = 10, double lambda = 1, double min_delta = 1e-6) {
  if (thread_no <= 0) {
    thread_no = SYS_THREADS_DEF;
  }

  int feature_no = S_r.n_rows;

  stdout_printf("Running ACTION (%d threads, with batch correction):", thread_no);
  FLUSH;

  if (k_max == -1) k_max = (int)S_r.n_cols;

  k_min = std::max(k_min, 2);
  k_max = std::min(k_max, (int)S_r.n_cols);

  ACTION_results trace;

  trace.H = field<mat>(k_max + 1);
  trace.C = field<mat>(k_max + 1);
  trace.selected_cols = field<uvec>(k_max + 1);

  mat X_r = normalise(S_r, 1);  // ATTENTION!

  int current_k = 0;
  char status_msg[50];

  sprintf(status_msg, "Iterating from k = %d ... %d:", k_min, k_max);
  stderr_printf("\n\t%s %d/%d finished", status_msg, current_k,
                (k_max - k_min + 1));
  FLUSH;

  parallelFor(k_min, k_max + 1, [&] (size_t kk) {
    SPA_results SPA_res = run_SPA(X_r, kk);
    trace.selected_cols[kk] = SPA_res.selected_columns;

    mat W = X_r.cols(trace.selected_cols[kk]);
    field<mat> AA_res = run_AA_with_batch_correction(X_r, W, batch, max_it, max_correction_rounds, lambda, min_delta);

    trace.C[kk] = AA_res(0);
    trace.H[kk] = AA_res(1);
    current_k++;

    stderr_printf("\r\t%s %d/%d finished", status_msg, current_k,
                  (k_max - k_min + 1));
    FLUSH;

  }, thread_no);
  stdout_printf("\r\t%s %d/%d finished\n", status_msg, current_k,
                (k_max - k_min + 1)); FLUSH;

  return trace;

}


}  // namespace ACTIONet
