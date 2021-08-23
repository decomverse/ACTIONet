#include <ACTIONet.h>
#include "cholmod.h"

#include <atomic>
#include <thread>

namespace ACTIONet
{
    template <class Function>
    inline void ParallelFor(size_t start, size_t end, size_t thread_no,
                            Function fn)
    {
        if (thread_no <= 0)
        {
            thread_no = std::thread::hardware_concurrency();
        }

        if (thread_no == 1)
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

            for (size_t threadId = 0; threadId < thread_no; ++threadId)
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

    mat normalize_scores(mat scores, int method = 1, int thread_no = 0)
    {
        mat normalized_scores(size(scores));
        switch (method)
        {
        case 0: //"none"
        {
            normalized_scores = scores;
            break;
        }
        case 1: //"zscore"
        {
            normalized_scores = zscore(scores, thread_no);
            break;
        }
        case 2: //"RINT" (nonparametric)
        {
            normalized_scores = RIN_transform(scores, thread_no);
            break;
        }
        case 3: //"robust_zscore" (kinda hack!)
        {
            normalized_scores = robust_zscore(scores, thread_no);
            break;
        }
        default:
            fprintf(stderr, "Unknown normalization method\n");
            normalized_scores = scores;
        }
        return (normalized_scores);
    }

    // G is the cell-cell network, scores is a cell x geneset matrix
    field<vec> autocorrelation_Moran(mat G, mat scores, int normalization_method, int perm_no, int thread_no)
    {
        int nV = G.n_rows;
        int scores_no = scores.n_cols;
        stdout_printf("Normalizizing scores (method=%d) ... ", normalization_method);
        mat normalized_scores = normalize_scores(scores, normalization_method, thread_no);
        stdout_printf("done\n");

        stdout_printf("Computing auto-correlation over network of %d samples for %d scores\n", nV, scores_no);
        double W = sum(sum(G));
        vec norm_sq = vec(trans(sum(square(normalized_scores))));
        vec norm_factors = nV / (W * norm_sq);

        vec stat = zeros(scores_no);
        ParallelFor(0, scores_no, thread_no, [&](size_t i, size_t threadId)
                    {
                        vec x = normalized_scores.col(i);
                        stat(i) = dot(x, G * x);
                    });

        vec mu = zeros(scores_no);
        vec sigma = zeros(scores_no);
        vec z = zeros(scores_no);
        if (0 < perm_no)
        {
            stdout_printf("Computing permutations ... ");

            mat rand_stats = zeros(scores_no, perm_no);
            ParallelFor(0, perm_no, thread_no, [&](size_t j, size_t threadId)
                        {
                            uvec perm = randperm(nV);
                            mat score_permuted = normalized_scores.rows(perm);

                            ParallelFor(0, scores_no, 1, [&](size_t i, size_t threadId)
                                        {
                                            vec rand_x = score_permuted.col(i);
                                            rand_stats(i, j) = dot(rand_x, G * rand_x);
                                        });
                        });
            stdout_printf("Done\n");

            mu = mean(rand_stats, 1);
            sigma = stddev(rand_stats, 0, 1);
            z = (stat - mu) / sigma;
            z.replace(datum::nan, 0);
        }
        // Summary stats
        stdout_printf("done\n");

        field<vec> results(4);
        results(0) = stat % norm_factors;
        results(1) = z;
        results(2) = mu;
        results(3) = sigma;
        return (results);
    }

    // G is the cell-cell network, scores is a cell x geneset matrix
    field<vec> autocorrelation_Moran(sp_mat G, mat scores, int normalization_method, int perm_no, int thread_no)
    {
        int nV = G.n_rows;
        int scores_no = scores.n_cols;

        stdout_printf("Normalizizing scores (method=%d) ... ", normalization_method);
        mat normalized_scores = normalize_scores(scores, normalization_method, thread_no);
        stdout_printf("done\n");

        stdout_printf("Computing auto-correlation over network of %d samples for %d scores\n", nV, scores_no);
        double W = sum(sum(G));
        vec norm_sq = vec(trans(sum(square(normalized_scores))));
        vec norm_factors = nV / (W * norm_sq);

        vec stat = zeros(scores_no);
        ParallelFor(0, scores_no, thread_no, [&](size_t i, size_t threadId)
                    {
                        vec x = normalized_scores.col(i);
                        stat(i) = dot(x, spmat_vec_product(G, x));
                    });

        vec mu = zeros(scores_no);
        vec sigma = zeros(scores_no);
        vec z = zeros(scores_no);
        if (0 < perm_no)
        {
            stdout_printf("Computing permutations ... ");

            mat rand_stats = zeros(scores_no, perm_no);
            ParallelFor(0, perm_no, thread_no, [&](size_t j, size_t threadId)
                        {
                            uvec perm = randperm(nV);
                            mat score_permuted = normalized_scores.rows(perm);

                            ParallelFor(0, scores_no, 1, [&](size_t i, size_t threadId)
                                        {
                                            vec rand_x = score_permuted.col(i);
                                            rand_stats(i, j) = dot(rand_x, spmat_vec_product(G, rand_x));
                                        });
                        });
            stdout_printf("Done\n");

            mu = mean(rand_stats, 1);
            sigma = stddev(rand_stats, 0, 1);
            z = (stat - mu) / sigma;
            z.replace(datum::nan, 0);
        }
        // Summary stats
        stdout_printf("done\n");

        field<vec> results(4);
        results(0) = stat % norm_factors;
        results(1) = z;
        results(2) = mu;
        results(3) = sigma;

        return (results);
    }

    // G is the cell-cell network, scores is a cell x geneset matrix
    field<vec> autocorrelation_Geary(mat G, mat scores, int normalization_method, int perm_no, int thread_no)
    {
        int nV = G.n_rows;
        int scores_no = scores.n_cols;
        stdout_printf("Normalizizing scores (method=%d) ... ", normalization_method);
        mat normalized_scores = normalize_scores(scores, normalization_method, thread_no);
        stdout_printf("done\n");

        stdout_printf("Computing auto-correlation over network of %d samples for %d scores\n", nV, scores_no);
        double W = sum(sum(G));
        vec norm_sq = vec(trans(sum(square(normalized_scores))));
        vec norm_factors = (nV - 1) / ((2 * W) * norm_sq);

        // Compute graph Laplacian
        vec d = vec(trans(sum(G)));
        mat L(-G);
        L.diag() = d;

        vec stat = zeros(scores_no);
        ParallelFor(0, scores_no, thread_no, [&](size_t i, size_t threadId)
                    {
                        vec x = normalized_scores.col(i);
                        stat(i) = dot(x, L * x);
                    });

        vec mu = zeros(scores_no);
        vec sigma = zeros(scores_no);
        vec z = zeros(scores_no);
        if (0 < perm_no)
        {
            stdout_printf("Computing permutations ... ");

            mat rand_stats = zeros(scores_no, perm_no);
            ParallelFor(0, perm_no, thread_no, [&](size_t j, size_t threadId)
                        {
                            uvec perm = randperm(nV);
                            mat score_permuted = normalized_scores.rows(perm);

                            ParallelFor(0, scores_no, 1, [&](size_t i, size_t threadId)
                                        {
                                            vec rand_x = score_permuted.col(i);
                                            rand_stats(i, j) = dot(rand_x, L * rand_x);
                                        });
                        });
            stdout_printf("Done\n");

            mu = mean(rand_stats, 1);
            sigma = stddev(rand_stats, 0, 1);
            z = (stat - mu) / sigma;
            z.replace(datum::nan, 0);
        }
        // Summary stats
        stdout_printf("done\n");

        field<vec> results(4);
        results(0) = stat % norm_factors;
        results(1) = -z;
        results(2) = mu;
        results(3) = sigma;
        return (results);
    }

    // G is the cell-cell network, scores is a cell x geneset matrix
    field<vec> autocorrelation_Geary(sp_mat G, mat scores, int normalization_method, int perm_no, int thread_no)
    {
        int nV = G.n_rows;
        int scores_no = scores.n_cols;

        stdout_printf("Normalizizing scores (method=%d) ... ", normalization_method);
        mat normalized_scores = normalize_scores(scores, normalization_method, thread_no);
        stdout_printf("done\n");

        stdout_printf("Computing auto-correlation over network of %d samples for %d scores\n", nV, scores_no);
        double W = sum(sum(G));
        vec norm_sq = vec(trans(sum(square(normalized_scores))));
        vec norm_factors = (nV - 1) / ((2 * W) * norm_sq);

        // Compute graph Laplacian
        vec d = vec(trans(sum(G)));
        sp_mat L(-G);
        L.diag() = d;

        vec stat = zeros(scores_no);
        ParallelFor(0, scores_no, thread_no, [&](size_t i, size_t threadId)
                    {
                        vec x = normalized_scores.col(i);
                        stat(i) = dot(x, spmat_vec_product(L, x));
                    });

        vec mu = zeros(scores_no);
        vec sigma = zeros(scores_no);
        vec z = zeros(scores_no);
        if (0 < perm_no)
        {
            stdout_printf("Computing permutations ... ");

            mat rand_stats = zeros(scores_no, perm_no);
            ParallelFor(0, perm_no, thread_no, [&](size_t j, size_t threadId)
                        {
                            uvec perm = randperm(nV);
                            mat score_permuted = normalized_scores.rows(perm);

                            ParallelFor(0, scores_no, 1, [&](size_t i, size_t threadId)
                                        {
                                            vec rand_x = score_permuted.col(i);
                                            rand_stats(i, j) = dot(rand_x, spmat_vec_product(L, rand_x));
                                        });
                        });
            stdout_printf("Done\n");

            mu = mean(rand_stats, 1);
            sigma = stddev(rand_stats, 0, 1);
            z = (stat - mu) / sigma;
            z.replace(datum::nan, 0);
        }
        // Summary stats
        stdout_printf("done\n");

        field<vec> results(4);
        results(0) = stat % norm_factors;
        results(1) = -z;
        results(2) = mu;
        results(3) = sigma;

        return (results);
    }

} // namespace ACTIONet
