#include <ACTIONet.h>

namespace ACTIONet
{

  field<mat> assess_enrichment(mat &scores, sp_mat &associations,
                               int thread_no = 0)
  {
    field<mat> res(3);

    if (scores.n_rows != associations.n_rows)
    {
      stderr_printf("Number of rows in scores and association matrices should both match the number of features\n");
      FLUSH;
      return (res);
    }

    associations = spones(associations);

    mat sorted_scores = sort(scores, "descend");
    vec a_max = trans(sorted_scores.row(0));
    umat perms(size(scores));
    for (int j = 0; j < scores.n_cols; j++)
    {
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
    parallelFor(
        0, associations.n_cols, [&](size_t k)
        {
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
        } },
        thread_no);

    field<mat> output(2);
    output(0) = logPvals;
    output(1) = thresholds;

    return (output);
  }

  vec xicor(vec xvec, vec yvec, bool compute_pval, int seed)
  {
    vec out(2);

    std::mt19937_64 engine(seed);

    vec idx = regspace(0, xvec.n_elem - 1);
    aarand::shuffle(idx.memptr(), idx.n_elem, engine);
    uvec perm = conv_to<uvec>::from(idx);

    xvec = xvec(perm);
    yvec = yvec(perm);

    int n = xvec.n_elem;

    vec er = rank_vec(xvec);
    vec fr = rank_vec(yvec, 1) / n;
    vec gr = rank_vec(-yvec, 1) / n;

    uvec ord = sort_index(er);
    fr = fr(ord);

    // Calculate xi
    double A1 = sum(abs(fr(span(0, n - 2)) - fr(span(1, n - 1)))) / (2.0 * n);
    double CU = mean(gr % (1.0 - gr));
    double xi = 1 - (A1 / CU);
    out(0) = xi;

    if (compute_pval == true)
    {
      // Calculate p-values
      vec qfr = sort(fr);
      vec ind = regspace(0, n - 1);
      vec ind2 = 2 * n - 2 * ind + 1;
      double ai = mean(ind2 % square(qfr)) / n;
      double ci = mean(ind2 % qfr) / n;
      vec cq = cumsum(qfr);
      vec m = (cq + (n - ind) % qfr) / n;
      double b = mean(square(m));
      double v = (ai - 2 * b + (ci * ci)) / (CU * CU);

      double xi_sd = sqrt(v / n);
      double z = sqrt(n) * xi / sqrt(v);

      out(1) = z;
    }
    else
    {
      out(1) = 0;
    }

    return (out);
  }

  field<mat> XICOR(mat &X, mat &Y, bool compute_pval, int seed, int thread_no)
  {
    field<mat> out(2);

    // arma_rng::set_seed(seed);
    // uvec perm = randperm(X.n_rows);

    // X = X.rows(perm);
    // Y = Y.rows(perm);

    bool swapped = false;
    int n1 = X.n_cols, n2 = Y.n_cols;
    if (n1 < n2)
    {
      swapped = true;
      mat Z = X;
      X = Y;
      Y = Z;
    }

    mat XI = zeros(X.n_cols, Y.n_cols);
    mat XI_Z = zeros(X.n_cols, Y.n_cols);

    // for(int k = 0; k < associations.n_cols; k++) {
    parallelFor(
        0, X.n_cols, [&](size_t i)
        {
            vec x = X.col(i);
            for(int j = 0; j < Y.n_cols; j++) {
              vec y = Y.col(j);
              vec out = xicor(x, y, compute_pval, seed);
              XI(i, j) = out(0);
              XI_Z(i, j) = out(1);
            } },
        thread_no);

    if (swapped)
    {
      XI = trans(XI);
      XI_Z = trans(XI_Z);
    }

    out(0) = XI;
    out(1) = XI_Z;

    return (out);
  }

} // namespace ACTIONet
