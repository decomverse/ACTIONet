#include "ACTIONet.h"
#include <cholmod.h>
sp_mat &as_arma_sparse(cholmod_sparse *chol_A, sp_mat &A,
                       cholmod_common *chol_c);

void dsdmult(char transpose, int m, int n, void *a, double *b, double *c,
             cholmod_common *chol_cp);

cholmod_sparse *as_cholmod_sparse(sp_mat &A, cholmod_sparse *chol_A,
                                  cholmod_common *chol_c);

namespace ACTIONet
{

  // Solves the standard Archetypal Analysis (AA) problem
  field<mat> run_AA(mat &A, mat &W0, int max_it, double min_delta)
  {
    int sample_no = A.n_cols;
    int d = A.n_rows;  // input dimension
    int k = W0.n_cols; // AA components

    mat C = zeros(sample_no, k);
    mat H = zeros(k, sample_no);

    mat W = W0;
    vec c(sample_no);

    double old_RSS = 0;

    for (int it = 0; it < max_it; it++)
    {
      mat C_old = C;
      mat H_old = H;
      double A_norm = norm(A, "fro");
      H = run_simplex_regression(W, A, true);
      // H = run_simplex_regression_FW(W, A, 30);

      mat R = A - W * H;
      mat Ht = trans(H);
      for (int i = 0; i < k; i++)
      {
        vec w = W.col(i);
        vec h = Ht.col(i);

        double norm_sq = arma::dot(h, h);
        if (norm_sq < double(10e-8))
        {
          // singular
          int max_res_idx = index_max(rowvec(sum(square(R), 0)));
          W.col(i) = A.col(max_res_idx);
          c.zeros();
          c(max_res_idx) = 1;
          C.col(i) = c;
        }
        else
        {
          vec b = w;
          cblas_dgemv(CblasColMajor, CblasNoTrans, R.n_rows, R.n_cols,
                      (1.0 / norm_sq), R.memptr(), R.n_rows, h.memptr(), 1, 1,
                      b.memptr(), 1);

          C.col(i) = run_simplex_regression(A, b, false);
          // C.col(i) = run_simplex_regression_FW(A, b, 30);

          vec w_new = A * C.col(i);
          vec delta = (w - w_new);

          // Rank-1 update: R += delta*h
          cblas_dger(CblasColMajor, R.n_rows, R.n_cols, 1.0, delta.memptr(), 1,
                     h.memptr(), 1, R.memptr(), R.n_rows);

          W.col(i) = w_new;
        }
      }
      double RSS = norm(R, "fro");
      double delta_RSS = abs(RSS - old_RSS) / A_norm;
      old_RSS = RSS;

      if (delta_RSS < min_delta)
        break;
    }

    C = clamp(C, 0, 1);
    C = normalise(C, 1);
    H = clamp(H, 0, 1);
    H = normalise(H, 1);

    field<mat> decomposition(2, 1);
    decomposition(0) = C;
    decomposition(1) = H;

    return decomposition;
  }

  // Solves separable NMF problem
  SPA_results run_SPA_rows_sparse(sp_mat &A, int k)
  {
    int m = A.n_rows;
    int n = A.n_cols;
    sp_mat A_sq = square(A);

    cholmod_common chol_c;
    cholmod_start(&chol_c);
    chol_c.final_ll = 1; /* LL' form of simplicial factorization */

    cholmod_sparse *AS = as_cholmod_sparse(A, AS, &chol_c);
    cholmod_sparse *AS_sq = as_cholmod_sparse(A_sq, AS_sq, &chol_c);

    SPA_results res;

    uvec K(k); // selected columns from A

    vec o = ones(n);
    vec normM(m);
    dsdmult('n', m, n, AS_sq, o.memptr(), normM.memptr(), &chol_c);
    vec normM1 = normM;
    mat U(n, k);

    vec norm_trace = zeros(k);
    double eps = 1e-6;
    for (int i = 0; i < k; i++)
    {
      // Find the column with maximum norm. In case of having more than one column
      // with almost very small diff in norm, pick the one that originally had the
      // largest norm
      double a = max(normM);
      norm_trace(i) = a;

      uvec b = find((a * ones(m, 1) - normM) / a <= eps);

      if (b.n_elem > 1)
      {
        uword idx = index_max(normM1(b));
        K(i) = b(idx);
      }
      else
      {
        K(i) = b(0);
      }

      // Pick row
      U.col(i) = vec(trans(A.row(K(i))));

      // Orthogonalize with respect to current basis
      for (int j = 0; j < i - 1; j++)
      {
        U.col(i) = U.col(i) - dot(U.col(j), U.col(i)) * U.col(j);
      }
      U.col(i) = U.col(i) / norm(U.col(i), 2);

      // Update column norms
      vec u = U.col(i);
      for (int j = i - 1; 0 <= j; j--)
      {
        u = u - dot(U.col(j), u) * U.col(j);
      }
      vec r(m);
      dsdmult('n', m, n, AS, u.memptr(), r.memptr(), &chol_c);

      uvec idx = find(U > 0);
      double perc = 100 * idx.n_elem / U.n_elem;
      stdout_printf("\t%d- res_norm = %f, U_density = %.2f%% (%d nnz)\n", i, a,
                    perc, idx.n_elem);
      FLUSH;

      normM = normM - (r % r);
    }

    res.selected_columns = K;
    res.column_norms = norm_trace;

    cholmod_free_sparse(&AS, &chol_c);
    cholmod_free_sparse(&AS_sq, &chol_c);
    cholmod_finish(&chol_c);

    return res;
  }

  // Solves separable NMF problem
  SPA_results run_SPA(mat &A, int k)
  {
    SPA_results res;

    int n = A.n_cols;
    uvec K(k); // selected columns from A
    K.zeros();

    rowvec normM = sum(square(A), 0);
    rowvec normM1 = normM;

    mat U(A.n_rows, k);

    vec norm_trace = zeros(k);
    double eps = 1e-16;

    for (int i = 1; i <= k; i++)
    {
      // Find the column with maximum norm. In case of having more than one column
      // with almost very small diff in norm, pick the one that originally had the
      // largest norm
      double a = max(normM);
      norm_trace(i - 1) = a;

      uvec b = find((a * ones(1, n) - normM) / a <= eps);
      if (b.n_elem == 0)
      {
        break;
      }
      else if (b.n_elem > 1)
      {
        uword idx = index_max(normM1(b));
        K(i - 1) = b(idx);
      }
      else
      {
        K(i - 1) = b(0);
      }

      // Pick column
      U.col(i - 1) = A.col(K(i - 1));

      // Orthogonalize with respect to current basis
      if (i > 1)
      {
        for (int j = 1; j <= i - 1; j++)
        {
          U.col(i - 1) =
              U.col(i - 1) - sum(U.col(j - 1) % U.col(i - 1)) * U.col(j - 1);
        }
      }
      double nm = norm(U.col(i - 1), 2);
      if (nm > 0)
        U.col(i - 1) /= nm;

      // Update column norms
      vec u = U.col(i - 1);
      if (i > 1)
      {
        for (int j = i - 1; 1 <= j; j--)
        {
          u = u - sum(U.col(j - 1) % u) * U.col(j - 1);
        }
      }
      normM = normM - square(u.t() * A);
      normM.transform([](double val)
                      { return (val < 0 ? 0 : val); });
    }

    res.selected_columns = K;
    res.column_norms = norm_trace;

    return res;
  }

  ACTION_results run_ACTION(mat &S_r, int k_min, int k_max, int thread_no,
                            int max_it = 100, double min_delta = 1e-6,
                            int normalization = 0)
  {
    if (thread_no <= 0)
    {
      thread_no = SYS_THREADS_DEF;
    }

    int feature_no = S_r.n_rows;

    stdout_printf("Running ACTION (%d threads):", thread_no);
    FLUSH;

    if (k_max == -1)
      k_max = (int)S_r.n_cols;

    k_min = std::max(k_min, 2);
    k_max = std::min(k_max, (int)S_r.n_cols);

    ACTION_results trace;

    trace.H = field<mat>(k_max + 1);
    trace.C = field<mat>(k_max + 1);
    trace.selected_cols = field<uvec>(k_max + 1);

    // ATTENTION!
    mat X_r = normalize_mat(S_r, normalization, 0);

    int current_k = 0;
    char status_msg[50];

    sprintf(status_msg, "Iterating from k = %d ... %d:", k_min, k_max);
    stderr_printf("\n\t%s %d/%d finished", status_msg, current_k,
                  (k_max - k_min + 1));
    FLUSH;

    parallelFor(
        k_min, k_max + 1,
        [&](size_t kk)
        {
          SPA_results SPA_res = run_SPA(X_r, kk);
          trace.selected_cols[kk] = SPA_res.selected_columns;

          mat W = X_r.cols(trace.selected_cols[kk]);

          field<mat> AA_res;

          AA_res = run_AA(X_r, W, max_it, min_delta);
          trace.C[kk] = AA_res(0);
          trace.H[kk] = AA_res(1);
          current_k++;

          stderr_printf("\r\t%s %d/%d finished", status_msg, current_k,
                        (k_max - k_min + 1));
          FLUSH;
        },
        thread_no);

    stdout_printf("\r\t%s %d/%d finished\n", status_msg, current_k,
                  (k_max - k_min + 1));
    FLUSH;

    return trace;
  }

  ACTION_results run_ACTION_plus(mat &S_r, int k_min, int k_max, int max_it = 100,
                                 double min_delta = 1e-6, int max_trial = 3)
  {
    stdout_printf("Running ACTION++ (%d threads):");
    FLUSH;

    int D = std::min((int)S_r.n_rows, (int)S_r.n_cols);
    if (k_max == -1)
      k_max = D;

    k_min = std::max(k_min, 2);
    k_max = std::min(k_max, D);

    ACTION_results trace;

    trace.H = field<mat>(k_max + 1, 1);
    trace.C = field<mat>(k_max + 1, 1);
    trace.selected_cols = field<uvec>(k_max + 1, 1);

    mat X_r = normalise(S_r, 1); // ATTENTION!
    SPA_results SPA_res = run_SPA(X_r, D);
    uvec selected_cols = SPA_res.selected_columns;

    mat W = mat(X_r.col(selected_cols(0)));

    field<mat> AA_res;
    int cur_idx = 0, jj, kk;
    stdout_printf("Iterating from k=%d ... %d (max trial = %d)\n", k_min, k_max,
                  max_trial);
    FLUSH;
    for (kk = k_min; kk <= k_max; kk++)
    {
      stdout_printf("\tk = %d\n", kk);
      FLUSH;

      for (jj = 0; jj < max_trial; jj++)
      {
        cur_idx++;
        stdout_printf("\t\tTrial %d: candidate %d = %d ... ", jj + 1, cur_idx + 1,
                      selected_cols(cur_idx));
        FLUSH;
        mat W_tmp = join_rows(W, X_r.col(selected_cols(cur_idx)));

        AA_res = run_AA(X_r, W_tmp, max_it, min_delta);

        vec influential_cells = vec(trans(sum(spones(sp_mat(AA_res(0))), 0)));
        int trivial_counts = (int)sum(influential_cells <= 1);

        if ((trivial_counts == 0))
        {
          stdout_printf("success\n");
          FLUSH;
          selected_cols(kk - 1) = selected_cols(cur_idx);
          break;
        }

        stdout_printf("failed\n");
        FLUSH;
        if ((cur_idx == (D - 1)))
        {
          stdout_printf("Reached end of the line!\n");
          FLUSH;
          break;
        }
      }

      if ((jj == max_trial) || (cur_idx == (D - 1)))
      {
        break;
      }

      trace.C[kk] = AA_res(0);
      trace.H[kk] = AA_res(1);
      trace.selected_cols(kk) = selected_cols(span(0, kk - 1));

      W = X_r * AA_res(0);
    }

    trace.C = trace.C.rows(0, kk - 1);
    trace.H = trace.H.rows(0, kk - 1);
    trace.selected_cols = trace.selected_cols.rows(0, kk - 1);

    return trace;
  }

  field<mat> run_AA_with_prior(mat &A, mat &W0, mat &W_prior, int max_it = 50,
                               double min_delta = 1e-16)
  {
    int sample_no = A.n_cols;
    int d = A.n_rows;  // input dimension
    int k = W0.n_cols; // AA components

    mat C = zeros(sample_no, k);
    mat H = zeros(k, sample_no);

    mat W = W0;
    vec c(sample_no);

    for (int it = 0; it < max_it; it++)
    {
      mat combined_W = join_rows(W, W_prior);
      mat combined_H = run_simplex_regression(combined_W, A, true);

      H = combined_H.rows(span(0, k - 1));

      mat R = A - W * H;
      mat Ht = trans(H);
      for (int i = 0; i < k; i++)
      {
        vec w = W.col(i);
        vec h = Ht.col(i);

        double norm_sq = arma::dot(h, h);
        if (norm_sq < double(10e-8))
        {
          // singular
          int max_res_idx = index_max(rowvec(sum(square(R), 0)));
          W.col(i) = A.col(max_res_idx);
          c.zeros();
          c(max_res_idx) = 1;
          C.col(i) = c;
        }
        else
        {
          vec b = w;
          cblas_dgemv(CblasColMajor, CblasNoTrans, R.n_rows, R.n_cols,
                      (1.0 / norm_sq), R.memptr(), R.n_rows, h.memptr(), 1, 1,
                      b.memptr(), 1);

          C.col(i) = run_simplex_regression(A, b, false);

          vec w_new = A * C.col(i);
          vec delta = (w - w_new);

          // Rank-1 update: R += delta*h
          cblas_dger(CblasColMajor, R.n_rows, R.n_cols, 1.0, delta.memptr(), 1,
                     h.memptr(), 1, R.memptr(), R.n_rows);

          W.col(i) = w_new;
        }
      }
    }

    C = clamp(C, 0, 1);
    C = normalise(C, 1);
    H = clamp(H, 0, 1);
    H = normalise(H, 1);

    field<mat> decomposition(2, 1);
    decomposition(0) = C;
    decomposition(1) = H;

    return decomposition;
  }

  ACTION_results run_subACTION(mat &S_r, mat &W_parent, mat &H_parent, int kk,
                               int k_min, int k_max, int thread_no,
                               int max_it = 50, double min_delta = 1e-16)
  {
    int feature_no = S_r.n_rows;

    stdout_printf("Running subACTION (%d threads) for parent archetype %d\n",
                  thread_no, kk + 1);
    FLUSH;

    if (k_max == -1)
      k_max = (int)S_r.n_cols;

    k_min = std::max(k_min, 2);
    k_max = std::min(k_max, (int)S_r.n_cols);

    mat X_r = normalise(S_r, 1); // ATTENTION!

    vec h = vec(trans(H_parent.row(kk)));
    mat W_prior = W_parent;
    W_prior.shed_col(kk);

    mat X_r_scaled = X_r; // To deflate or not deflate!

    for (int i = 0; i < X_r_scaled.n_cols; i++)
    {
      X_r_scaled.col(i) *= h[i];
    }

    ACTION_results trace;
    trace.H = field<mat>(k_max + 1);
    trace.C = field<mat>(k_max + 1);
    trace.selected_cols = field<uvec>(k_max + 1);

    int current_k = 0;
    char status_msg[50];

    sprintf(status_msg, "Iterating from k = %d ... %d:", k_min, k_max);
    stderr_printf("\n\t%s %d/%d finished", status_msg, current_k,
                  (k_max - k_min + 1));
    FLUSH;

    stdout_printf("Iterating from k=%d ... %d\n", k_min, k_max);
    FLUSH;
    parallelFor(
        k_min, k_max + 1,
        [&](size_t kkk)
        {
          SPA_results SPA_res = run_SPA(X_r_scaled, kkk);
          trace.selected_cols[kkk] = SPA_res.selected_columns;

          mat W = X_r.cols(trace.selected_cols[kkk]);
          field<mat> AA_res;
          AA_res = run_AA_with_prior(X_r_scaled, W, W_prior, max_it, min_delta);

          trace.C[kkk] = AA_res(0);
          trace.H[kkk] = AA_res(1);
          current_k++;

          stderr_printf("\r\t%s %d/%d finished", status_msg, current_k,
                        (k_max - k_min + 1));
          FLUSH;
        },
        thread_no);

    stdout_printf("\r\t%s %d/%d finished\n", status_msg, current_k,
                  (k_max - k_min + 1));
    FLUSH;

    return trace;
  }

  /* alpha: sums to one and indicates relative importance of each view
   *
   */
  void findConsensus(vector<mat> S, full_trace &run_trace, int arch_no, vec alpha,
                     double lambda, int max_it, int)
  {
    int i;
    int ds_no = run_trace.indiv_trace[arch_no]
                    .H_primary.size(); // number of datasets ( ~ 2)
    int cell_no = S[0].n_cols;
    vec c(cell_no);

    // make sure it's a convex vector
    alpha.transform([](double val)
                    { return (max(0.0, val)); });
    alpha = normalise(alpha, 1);

    run_trace.indiv_trace[arch_no].C_secondary[0] =
        run_trace.indiv_trace[arch_no].C_primary[0];
    run_trace.indiv_trace[arch_no].H_secondary[0] =
        run_trace.indiv_trace[arch_no].H_primary[0];
    for (int ds = 1; ds < ds_no; ds++)
    {
      mat G = 1 + cor(trans(run_trace.indiv_trace[arch_no].H_primary[0]),
                      trans(run_trace.indiv_trace[arch_no].H_primary[ds]));
      mat G_matched = MWM_hungarian(G);
      uvec perm = index_max(G_matched, 1);

      run_trace.indiv_trace[arch_no].C_secondary[ds] =
          run_trace.indiv_trace[arch_no].C_primary[ds].cols(perm);
      run_trace.indiv_trace[arch_no].H_secondary[ds] =
          run_trace.indiv_trace[arch_no].H_primary[ds].rows(perm);
    }

    // Compute initial H_consensus
    run_trace.H_consensus[arch_no] =
        zeros(size(run_trace.indiv_trace[arch_no].H_primary[0]));
    for (int ds = 0; ds < ds_no; ds++)
    {
      run_trace.H_consensus[arch_no] +=
          alpha(ds) * run_trace.indiv_trace[arch_no].H_secondary[ds];
    }

    // Estimate relative ratio of error terms
    mat H_hat = run_trace.H_consensus[arch_no];
    double a = 0.0, b = 0.0, x, y;
    for (int ds = 0; ds < ds_no; ds++)
    {
      mat W = run_trace.indiv_trace[arch_no].C_secondary[ds];
      mat H = run_trace.indiv_trace[arch_no].H_secondary[ds];

      x = norm(S[ds] - S[ds] * W * H, "fro");
      y = norm(H - H_hat, "fro");
      a += (x * x);
      b += (alpha(ds) * y * y);
    }
    double ratio = a / b;
    lambda *= ratio;

    // Main loop
    for (int it = 0; it < max_it; it++)
    {
      // Permute rows
      for (int ds = 1; ds < ds_no; ds++)
      {
        mat G = 1 + cor(trans(run_trace.indiv_trace[arch_no].H_secondary[0]),
                        trans(run_trace.indiv_trace[arch_no].H_secondary[ds]));
        mat G_matched = MWM_hungarian(G);
        uvec perm = index_max(G_matched, 1);

        run_trace.indiv_trace[arch_no].C_secondary[ds] =
            run_trace.indiv_trace[arch_no].C_secondary[ds].cols(perm);
        run_trace.indiv_trace[arch_no].H_secondary[ds] =
            run_trace.indiv_trace[arch_no].H_secondary[ds].rows(perm);
      }

      // Compute shared subspace
      run_trace.H_consensus[arch_no] =
          zeros(size(run_trace.indiv_trace[arch_no].H_primary[0]));
      for (int ds = 0; ds < ds_no; ds++)
      {
        run_trace.H_consensus[arch_no] +=
            alpha(ds) * run_trace.indiv_trace[arch_no].H_secondary[ds];
      }

      // Recompute H_i
      for (int ds = 0; ds < ds_no; ds++)
      {
        mat I = eye(arch_no, arch_no);
        mat Z =
            S[ds] *
            run_trace.indiv_trace[arch_no].C_secondary[ds]; // Archetype matrix
        double weight = lambda * alpha[ds];

        mat A = join_vert(trans(Z) * Z, weight * I);
        mat B =
            join_vert(trans(Z) * S[ds], weight * run_trace.H_consensus[arch_no]);

        run_trace.indiv_trace[arch_no].H_secondary[ds] =
            run_simplex_regression(A, B, false);
      }

      // Recompute C_i
      for (int ds = 0; ds < ds_no; ds++)
      {
        mat W = S[ds] * run_trace.indiv_trace[arch_no].C_secondary[ds];
        mat H = run_trace.indiv_trace[arch_no].H_secondary[ds];
        mat R = S[ds] - W * H;
        for (int j = 0; j < arch_no; j++)
        {
          double norm_sq = sum(square(H.row(j)));
          vec h = trans(H.row(j)) / norm_sq;
          vec b = R * h + W.col(j);

          c = run_simplex_regression(S[ds], b, false);

          R += (W.col(j) - S[ds] * c) * H.row(j);
          W.col(j) = S[ds] * c;
          run_trace.indiv_trace[arch_no].C_secondary[ds].col(j) = c;
        }
      }
    }
  }

  full_trace runACTION_muV(vector<mat> S_r, int k_min, int k_max, vec alpha,
                           double lambda, int AA_iters, int Opt_iters,
                           int thread_no)
  {
    stdout_printf("Running ACTION muV (%d threads):", thread_no);
    FLUSH;

    double lambda2 = 1e-5, epsilon = 1e-5;

    k_min = std::max(k_min, 2);
    k_max = std::min(k_max, (int)(S_r[0].n_cols));

    int cell_no = S_r[0].n_cols;
    vec c(cell_no);

    full_trace run_trace;
    run_trace.H_consensus.resize(k_max + 1);
    run_trace.indiv_trace.resize(k_max + 1);
    for (int kk = 0; kk <= k_max; kk++)
    {
      run_trace.indiv_trace[kk].selected_cols.resize(S_r.size());
      run_trace.indiv_trace[kk].H_primary.resize(S_r.size());
      run_trace.indiv_trace[kk].C_primary.resize(S_r.size());
      run_trace.indiv_trace[kk].H_secondary.resize(S_r.size());
      run_trace.indiv_trace[kk].C_secondary.resize(S_r.size());
      run_trace.indiv_trace[kk].C_consensus.resize(S_r.size());
    }

    // Normalize signature profiles
    for (int i = 0; i < S_r.size(); i++)
    {
      S_r[i] = normalise(S_r[i], 1, 0); // norm-1 normalize across columns --
                                        // particularly important for SPA
    }

    field<mat> AA_res(2, 1);
    int current_k = 0;
    char status_msg[50];

    sprintf(status_msg, "Iterating from k = %d ... %d:", k_min, k_max);
    stderr_printf("\n\t%s %d/%d finished", status_msg, current_k,
                  (k_max - k_min + 1));
    FLUSH;

    parallelFor(
        k_min, k_max + 1,
        [&](size_t kk)
        {
          // Solve ACTION for a fixed-k to "jump-start" the joint optimization
          // problem.
          for (int i = 0; i < S_r.size(); i++)
          {
            SPA_results SPA_res = run_SPA(S_r[i], kk);
            run_trace.indiv_trace[kk].selected_cols[i] = SPA_res.selected_columns;

            mat W = S_r[i].cols(run_trace.indiv_trace[kk].selected_cols[i]);

            AA_res = run_AA(S_r[i], W, AA_iters);

            mat C0 = AA_res(0);
            C0.transform([](double val)
                         { return (min(1.0, max(0.0, val))); });
            C0 = normalise(C0, 1);
            run_trace.indiv_trace[kk].C_primary[i] = C0;

            mat H0 = AA_res(1);
            H0.transform([](double val)
                         { return (min(1.0, max(0.0, val))); });
            H0 = normalise(H0, 1);
            run_trace.indiv_trace[kk].H_primary[i] = H0;
          }

          // Compute consensus latent subspace, H^*
          findConsensus(S_r, run_trace, kk, alpha, lambda, Opt_iters,
                        thread_no); // sets secondary and consensus objects

          // decouple to find individual consensus C matrices
          for (int i = 0; i < S_r.size(); i++)
          {
            mat S = S_r[i];
            mat C = run_trace.indiv_trace[kk].C_secondary[i];
            mat H = run_trace.indiv_trace[kk].H_secondary[i];
            mat W = S * C;

            mat R = S - W * H;

            run_trace.indiv_trace[kk].C_consensus[i] = zeros(cell_no, kk);
            for (int j = 0; j < kk; j++)
            {
              vec w = W.col(j);
              vec h = trans(H.row(j));

              double norm_sq = arma::dot(h, h);
              if (norm_sq < double(10e-8))
              {
                // singular
                int max_res_idx = index_max(rowvec(sum(square(R), 0)));
                W.col(j) = S.col(max_res_idx);
                c.zeros();
                c(max_res_idx) = 1;
                C.col(j) = c;
              }
              else
              {
                vec b = w;
                cblas_dgemv(CblasColMajor, CblasNoTrans, R.n_rows, R.n_cols,
                            (1.0 / norm_sq), R.memptr(), R.n_rows, h.memptr(), 1,
                            1, b.memptr(), 1);

                C.col(j) = run_simplex_regression(S, b, false);

                vec w_new = S * C.col(j);
                vec delta = (w - w_new);

                // Rank-1 update: R += delta*h
                cblas_dger(CblasColMajor, R.n_rows, R.n_cols, 1.0, delta.memptr(),
                           1, h.memptr(), 1, R.memptr(), R.n_rows);

                W.col(j) = w_new;
              }

              run_trace.indiv_trace[kk].C_consensus[i].col(j) = C.col(j);
            }
          }
          current_k++;
          stderr_printf("\r\t%s %d/%d finished", status_msg, current_k,
                        (k_max - k_min + 1));
          FLUSH;
        },
        thread_no);

    stdout_printf("\r\t%s %d/%d finished\n", status_msg, current_k,
                  (k_max - k_min + 1));
    FLUSH;

    return run_trace;
  }

} // namespace ACTIONet
