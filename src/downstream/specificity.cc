#include <ACTIONet.h>
#include "cholmod.h"

#include <atomic>
#include <thread>

template <class Function>
inline void ParallelFor(size_t start, size_t end, size_t numThreads,
                        Function fn)
{
  if (numThreads <= 0)
  {
    numThreads = SYS_THREADS_DEF;
  }

  if (numThreads == 1)
  {
    for (size_t id = start; id < end; id++)
    {
      fn(id, 0);
    }
  }
  else
  {
    std::vector<std::thread> threads;
    std::atomic<size_t> current(start);

    // keep track of exceptions in threads
    // https://stackoverflow.com/a/32428427/1713196
    std::exception_ptr lastException = nullptr;
    std::mutex lastExceptMutex;

    for (size_t threadId = 0; threadId < numThreads; ++threadId)
    {
      threads.push_back(std::thread([&, threadId]
                                    {
                                      while (true)
                                      {
                                        size_t id = current.fetch_add(1);

                                        if ((id >= end))
                                        {
                                          break;
                                        }

                                        try
                                        {
                                          fn(id, threadId);
                                        }
                                        catch (...)
                                        {
                                          std::unique_lock<std::mutex> lastExcepLock(lastExceptMutex);
                                          lastException = std::current_exception();
                                          /*
             * This will work even when current is the largest value that
             * size_t can fit, because fetch_add returns the previous value
             * before the increment (what will result in overflow
             * and produce 0 instead of current + 1).
             */
                                          current = end;
                                          break;
                                        }
                                      }
                                    }));
    }
    for (auto &thread : threads)
    {
      thread.join();
    }
    if (lastException)
    {
      std::rethrow_exception(lastException);
    }
  }
}

namespace ACTIONet
{
  field<mat> compute_feature_specificity_bin(sp_mat &Sb, mat &H, int thread_no = 0)
  {
    stdout_printf("Computing feature specificity ... ");
    field<mat> res(3);

    mat Ht = trans(H);
    Ht.each_col([](vec &h)
                {
                  double mu = mean(h);
                  h /= (mu == 0) ? 1 : mu;
                }); // For numerical stability

    // printf("Compute stats ... ");
    // vec row_p = vec(mean(Sb, 1));
    // rowvec col_p = rowvec(mean(Sb, 0));

    // Heuristic optimization! Shall add parallel for later on
    vec row_p = zeros(Sb.n_rows);
    vec col_p = zeros(Sb.n_cols);
    sp_mat::const_iterator it = Sb.begin();
    sp_mat::const_iterator it_end = Sb.end();
    for (; it != it_end; ++it)
    {
      row_p[it.row()]++;
      col_p[it.col()]++;
    }
    row_p /= Sb.n_cols;
    col_p /= Sb.n_rows;
    // printf("done\n");

    // printf("Computing observation statistics ... ");
    cholmod_common chol_c;
    cholmod_start(&chol_c);
    chol_c.final_ll = 1; /* LL' form of simplicial factorization */
    cholmod_sparse_struct *Schol = new cholmod_sparse_struct;
    as_cholmod_sparse(Schol, Sb);
    mat Obs = zeros(Sb.n_rows, Ht.n_cols);
    ParallelFor(0, Ht.n_cols, thread_no, [&](size_t i, size_t threadId)
                { dsdmult('n', Sb.n_rows, Sb.n_cols, Schol, Ht.colptr(i), Obs.colptr(i), &chol_c); });
    // Free up matrices
    cholmod_free_sparse(&Schol, &chol_c);
    cholmod_finish(&chol_c);
    // printf("done\n");

    // printf("Computing expectation statistics ... ");
    double rho = mean(col_p);
    vec beta = col_p / rho; // Relative density compared to the overall density
    mat Gamma = Ht;
    vec a(H.n_rows);
    for (int i = 0; i < H.n_rows; i++)
    {
      Gamma.col(i) %= beta;
      a(i) = max(Gamma.col(i));
    }

    mat Exp = row_p * sum(Gamma, 0);
    mat Nu = row_p * sum(square(Gamma), 0);
    // printf("done\n");

    // printf("Computing significance ... ");
    mat Lambda = Obs - Exp;

    mat logPvals_lower = square(Lambda) / (2 * Nu);
    uvec uidx = find(Lambda >= 0);
    logPvals_lower(uidx) = zeros(uidx.n_elem);
    logPvals_lower.replace(datum::nan, 0); // replace each NaN with 0

    mat Lambda_scaled = Lambda;
    for (int j = 0; j < Lambda_scaled.n_cols; j++)
    {
      Lambda_scaled.col(j) *= (a(j) / 3);
    }
    mat logPvals_upper = square(Lambda) / (2 * (Nu + Lambda_scaled));
    uvec lidx = find(Lambda <= 0);
    logPvals_upper(lidx) = zeros(lidx.n_elem);
    logPvals_upper.replace(datum::nan, 0); // replace each NaN with 0

    logPvals_lower /= log(10);
    logPvals_upper /= log(10);
    stdout_printf("done\n");

    res(0) = Obs / Ht.n_rows;
    res(1) = logPvals_upper;
    res(2) = logPvals_lower;

    return (res);
  }

  field<mat> compute_feature_specificity(sp_mat &S, mat &H, int thread_no = 0)
  {
    stdout_printf("Computing feature specificity ... ");
    field<mat> res(3);

    mat Ht = trans(H);
    Ht.each_col([](vec &h)
                {
                  double mu = mean(h);
                  h /= (mu == 0) ? 1 : mu;
                }); // For numerical stability

    // make sure all values are positive
    double min_val = S.min();
    S.for_each([min_val](mat::elem_type &val)
               { val -= min_val; });

    //stdout_printf("Compute stats ... ");
    // sp_mat Sb = spones(S);

    // Heuristic optimization! Shall add parallel for later on
    vec row_p = zeros(S.n_rows);
    vec col_p = zeros(S.n_cols);
    vec row_factor = zeros(S.n_rows);

    sp_mat::const_iterator it = S.begin();
    sp_mat::const_iterator it_end = S.end();
    for (; it != it_end; ++it)
    {
      col_p[it.col()]++;
      row_p[it.row()]++;
      row_factor[it.row()] += (*it);
    }
    row_factor /= row_p;
    row_p /= S.n_cols;
    col_p /= S.n_rows;

    //stdout_printf("done\n");

    //stdout_printf("Computing observation statistics ... ");
    cholmod_common chol_c;
    cholmod_start(&chol_c);
    chol_c.final_ll = 1; /* LL' form of simplicial factorization */
    cholmod_sparse_struct *Schol = new cholmod_sparse_struct;
    as_cholmod_sparse(Schol, S);
    mat Obs = zeros(S.n_rows, Ht.n_cols);
    ParallelFor(0, Ht.n_cols, thread_no, [&](size_t i, size_t threadId)
                { dsdmult('n', S.n_rows, S.n_cols, Schol, Ht.colptr(i), Obs.colptr(i), &chol_c); });
    // Free up matrices
    cholmod_free_sparse(&Schol, &chol_c);
    cholmod_finish(&chol_c);
    //stdout_printf("done\n");

    //stdout_printf("Computing expectation statistics ... ");
    double rho = mean(col_p);
    vec beta = col_p / rho; // Relative density compared to the overall density
    mat Gamma = Ht;
    vec a(H.n_rows);
    for (int i = 0; i < H.n_rows; i++)
    {
      Gamma.col(i) %= beta;
      a(i) = max(Gamma.col(i));
    }

    mat Exp = (row_p % row_factor) * sum(Gamma, 0);
    mat Nu = (row_p % square(row_factor)) * sum(square(Gamma), 0);
    mat A = (row_factor * trans(a));
    //stdout_printf("done\n");

    //stdout_printf("Computing significance ... ");
    mat Lambda = Obs - Exp;

    mat logPvals_lower = square(Lambda) / (2 * Nu);
    uvec uidx = find(Lambda >= 0);
    logPvals_lower(uidx) = zeros(uidx.n_elem);
    logPvals_lower.replace(datum::nan, 0); // replace each NaN with 0

    mat logPvals_upper = square(Lambda) / (2 * (Nu + (Lambda % A / 3)));
    uvec lidx = find(Lambda <= 0);
    logPvals_upper(lidx) = zeros(lidx.n_elem);
    logPvals_upper.replace(datum::nan, 0); // replace each NaN with 0

    logPvals_lower /= log(10);
    logPvals_upper /= log(10);
    stdout_printf("done\n");

    res(0) = Obs / Ht.n_rows;
    res(1) = logPvals_upper;
    res(2) = logPvals_lower;

    return (res);
  }

  field<mat> compute_feature_specificity(mat &S, mat &H, int thread_no = 0)
  {
    stdout_printf("Computing feature specificity ... ");
    field<mat> res(3);

    // make sure all values are positive
    double min_val = S.min();
    S.for_each([min_val](mat::elem_type &val)
               { val -= min_val; });

    mat Sb = S;
    uvec nnz_idx = find(Sb > 0);
    (Sb(nnz_idx)).ones();

    mat Ht = trans(H);
    Ht.each_col([](vec &h)
                {
                  double mu = mean(h);
                  h /= (mu == 0) ? 1 : mu;
                }); // For numerical stability

    // printf("Compute stats ... ");
    vec row_p = vec(sum(Sb, 1));
    vec row_factor = vec(sum(S, 1)) / row_p; // mean of nonzero elements
    row_p /= Sb.n_cols;
    vec col_p = vec(trans(mean(Sb, 0)));
    // printf("done\n");

    // printf("Computing observation statistics ... ");
    mat Obs = zeros(S.n_rows, Ht.n_cols);
    ParallelFor(0, Ht.n_cols, thread_no, [&](size_t i, size_t threadId)
                { Obs.col(i) = S * Ht.col(i); });
    // printf("done\n");

    // printf("Computing expectation statistics ... ");
    double rho = mean(col_p);
    vec beta = col_p / rho; // Relative density compared to the overall density
    mat Gamma = Ht;
    vec a(H.n_rows);
    for (int i = 0; i < H.n_rows; i++)
    {
      Gamma.col(i) %= beta;
      a(i) = max(Gamma.col(i));
    }

    mat Exp = (row_p % row_factor) * sum(Gamma, 0);
    mat Nu = (row_p % square(row_factor)) * sum(square(Gamma), 0);
    mat A = (row_factor * trans(a));
    // printf("done\n");

    // printf("Computing significance ... ");
    mat Lambda = Obs - Exp;

    mat logPvals_lower = square(Lambda) / (2 * Nu);
    uvec uidx = find(Lambda >= 0);
    logPvals_lower(uidx) = zeros(uidx.n_elem);
    logPvals_lower.replace(datum::nan, 0); // replace each NaN with 0

    mat logPvals_upper = square(Lambda) / (2 * (Nu + (Lambda % A / 3)));
    uvec lidx = find(Lambda <= 0);
    logPvals_upper(lidx) = zeros(lidx.n_elem);
    logPvals_upper.replace(datum::nan, 0); // replace each NaN with 0

    logPvals_lower /= log(10);
    logPvals_upper /= log(10);
    stdout_printf("done\n");

    res(0) = Obs / Ht.n_rows;
    res(1) = logPvals_upper;
    res(2) = logPvals_lower;

    return (res);
  }

  field<mat> compute_feature_specificity(sp_mat &S, uvec sample_assignments, int thread_no = 0)
  {
    mat H(max(sample_assignments), S.n_cols);

    for (int i = 1; i <= max(sample_assignments); i++)
    {
      vec v = zeros(S.n_cols);
      uvec idx = find(sample_assignments == i);
      v(idx) = ones(idx.n_elem);
      H.row(i - 1) = trans(v);
    }

    field<mat> res = compute_feature_specificity(S, H, thread_no);

    return (res);
  }

  field<mat> compute_feature_specificity(mat &S, uvec sample_assignments, int thread_no = 0)
  {
    mat H(max(sample_assignments), S.n_cols);

    for (int i = 1; i <= max(sample_assignments); i++)
    {
      vec v = zeros(S.n_cols);
      uvec idx = find(sample_assignments == i);
      v(idx) = ones(idx.n_elem);
      H.row(i - 1) = trans(v);
    }

    field<mat> res = compute_feature_specificity(S, H, thread_no);

    return (res);
  }

} // namespace ACTIONet
