#include <ACTIONet.h>
#include "cholmod.h"

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

    // G is the cell-cell network, scores is a cell x geneset matrix
    field<vec> computeAutocorrelation_Geary_sparse(sp_mat &G, mat &scores, int perm_no = 30, int thread_no = 0)
    {
        int nV = G.n_rows;
        int feature_set_no = scores.n_cols;

        printf("Computing auto-correlation over network of %d samples for %d scores\n", nV, feature_set_no);
        int nnz = G.n_nonzero;
        int idx = 0;

        vec vals(nnz);
        umat subs(2, nnz);
        for (sp_mat::iterator it = G.begin(); it != G.end(); ++it)
        {
            vals(idx) = *it;
            subs(0, idx) = it.row();
            subs(1, idx) = it.col();
            idx++;
        }

        //double scale_factor = sum(sum(G)) / (sample_no-1); // 2W / (N-1)
        double total_weight = sum(vals);

        // Compute graph Laplacian
        vec d = vec(trans(sum(G)));

        sp_mat L(nV, nV);
        L.diag() = d;
        L -= G;

        cholmod_common chol_c;
        cholmod_start(&chol_c);
        chol_c.final_ll = 1; /* LL' form of simplicial factorization */
        cholmod_sparse_struct *Ls = new cholmod_sparse_struct;
        as_cholmod_sparse(Ls, L);

        printf("Computing autocorrelations ...");
        fflush(stdout);
        vec Cstat = zeros(feature_set_no);
        ParallelFor(0, feature_set_no, thread_no, [&](size_t i, size_t threadId)
                    {
                        vec x = scores.col(i);
                        vec Lx(L.n_rows);

                        dsdmult('n', L.n_rows, L.n_cols, Ls, x.colptr(), Lx.colptr(i), &chol_c);
                        double stat = dot(x, Lx);
                        double norm_fact = var(x) * total_weight;
                        Cstat(i) = 1 - (stat / norm_fact);
                    });
        // Free up matrices
        if (0 <= perm_no)
        {
            mat Cstat_rand = zeros(feature_set_no);
            ParallelFor(0, perm_no, thread_no, [&](size_t j, size_t threadId)
                        {
                            uvec perm = randperm(scores.n_row);
                            mat score_permuted = scores.rows(perm);

                            ParallelFor(0, feature_set_no, 1, [&](size_t i, size_t threadId)
                                        {
                                            vec x = score_permuted.col(i);
                                            vec Lx(L.n_rows);

                                            dsdmult('n', L.n_rows, L.n_cols, Ls, x.colptr(), Lx.colptr(i), &chol_c);
                                            double stat = dot(x, Lx);
                                            double norm_fact = var(x) * total_weight;
                                            Cstat_rand(i, j) = 1 - (stat / norm_fact);
                                        });
                        });
            mu = trans(mean(Cstat_rand));
            sigma = trans(stddev(Cstat_rand));
            Cstat_Z = (Cstat - mu) / sigma;
        }
        else
        {
            mu = zeros(feature_set_no);
            sigma = zeros(feature_set_no);
            Cstat_Z = zeros(feature_set_no);
        }
        cholmod_free_sparse(&Ls, &chol_c);
        cholmod_finish(&chol_c);

        // Summary stats
        printf("done\n");
        fflush(stdout);

        field<vec> results(4);
        results(0) = Cstat;
        results(1) = mu;
        results(2) = sigma;
        results(3) = Cstat_Z;

        return (results);
    }
} // namespace ACTIONet
