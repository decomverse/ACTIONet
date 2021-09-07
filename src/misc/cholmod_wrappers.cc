#include "ACTIONet.h"
#include <cholmod.h>

#include <atomic>
#include <thread>

std::mutex mtx; // mutex for critical section

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

    // Mat-vec product (Ax)
    //dsdmult('n', A.n_rows, A.n_cols, A_as_cholmod, x.memptr(), out.memptr(), &chol_c);
    void dsdmult(char transpose, int n_rows, int n_cols, void *A, double *x, double *out,
                 cholmod_common *chol_cp)
    {
        int t = transpose == 't' ? 1 : 0; // 'n': computes Ax, 't': computes A'x
        cholmod_sparse *cha = (cholmod_sparse *)A;

        // x
        cholmod_dense chb;
        chb.nrow = t ? n_rows : n_cols;
        chb.d = chb.nrow;
        chb.ncol = 1;
        chb.nzmax = chb.nrow;
        chb.xtype = cha->xtype;
        chb.dtype = 0;
        chb.x = (void *)x;
        chb.z = (void *)NULL;

        // out
        cholmod_dense chc;
        chc.nrow = t ? n_cols : n_rows;
        chc.d = chc.nrow;
        chc.ncol = 1;
        chc.nzmax = chc.nrow;
        chc.xtype = cha->xtype;
        chc.dtype = 0;
        chc.x = (void *)out;
        chc.z = (void *)NULL;

        double one[] = {1, 0}, zero[] = {0, 0};

        cholmod_sdmult(cha, t, one, zero, &chb, &chc, chol_cp);
    }

    cholmod_sparse *as_cholmod_sparse(arma::sp_mat &A, cholmod_sparse *chol_A,
                                      cholmod_common *chol_c)
    {
        int nrow = A.n_rows, ncol = A.n_cols, nz = A.n_nonzero;
        cholmod_allocate_work(0, max(nrow, ncol), 0, chol_c);
        chol_A = cholmod_allocate_sparse(nrow, ncol, nz, 1 /*sorted*/, 1 /*packed*/, 0 /*NOT symmetric*/, CHOLMOD_REAL, chol_c);

        int *ptr = (int *)chol_A->p;
        double *x_ptr = (double *)chol_A->x;
        int *i_ptr = (int *)chol_A->i;

        mtx.lock();
        A.sync();
        {
            for (int k = 0; k < A.n_nonzero; k++)
            {
                x_ptr[k] = (A.values)[k];
                i_ptr[k] = (A.row_indices)[k];
            }
            for (int k = 0; k < A.n_cols + 1; k++)
            {
                ptr[k] = A.col_ptrs[k];
            }
        }
        mtx.unlock();

        return chol_A;
    }

    sp_mat &as_arma_sparse(cholmod_sparse *chol_A, arma::sp_mat &A, cholmod_common *chol_c)
    {
        // Allocate space
        A = arma::sp_mat(chol_A->nrow, chol_A->ncol);
        A.mem_resize(static_cast<unsigned>(chol_A->nzmax));

        mtx.lock();
        A.sync();
        double *in_x_ptr = (double *)chol_A->x;
        int *in_i_ptr = (int *)chol_A->i;
        double *out_x_ptr = (double *)arma::access::rwp(A.values);
        arma::uword *out_i_ptr = arma::access::rwp(A.row_indices);
        for (int k = 0; k < chol_A->nzmax; k++)
        {
            out_x_ptr[k] = in_x_ptr[k];
            out_i_ptr[k] = in_i_ptr[k];
        }

        int *in_p_ptr = (int *)chol_A->p;
        arma::uword *out_p_ptr = arma::access::rwp(A.col_ptrs);
        for (int k = 0; k < chol_A->ncol; k++)
        {
            out_p_ptr[k] = in_p_ptr[k];
        }

        // important: set the sentinel as well
        arma::access::rwp(A.col_ptrs)[chol_A->ncol] = chol_A->nzmax;

        // set the number of non-zero elements
        arma::access::rw(A.n_nonzero) = chol_A->nzmax;
        mtx.unlock();

        return A;
    }

    vec spmat_vec_product(sp_mat &A, vec &x)
    {
        mat X = mat(x);
        mat Ax = spmat_mat_product(A, X);

        return (Ax.col(0));
    }

    mat spmat_mat_product(sp_mat &A, mat &B)
    {
        cholmod_common chol_c;
        cholmod_start(&chol_c);

        cholmod_sparse *chol_A = as_cholmod_sparse(A, chol_A, &chol_c);

        cholmod_dense *chol_B = cholmod_allocate_dense(B.n_rows, B.n_cols, B.n_rows, CHOLMOD_REAL, &chol_c);
        chol_B->x = (void *)B.memptr();
        chol_B->z = (void *)NULL;

        mat res = zeros(A.n_rows, B.n_cols);
        cholmod_dense *out = cholmod_allocate_dense(A.n_rows, B.n_cols, A.n_rows, CHOLMOD_REAL, &chol_c);
        out->x = (void *)res.memptr();
        out->z = (void *)NULL;

        double one[] = {1, 0}, zero[] = {0, 0};
        cholmod_sdmult(chol_A, 0, one, zero, chol_B, out, &chol_c);

        cholmod_free_sparse(&chol_A, &chol_c);
        cholmod_finish(&chol_c);
        return (res);
    }

    sp_mat spmat_spmat_product(sp_mat &A, sp_mat &B)
    {
        arma::sp_mat res;

        if (A.n_cols != B.n_rows)
        {
            fprintf(stderr, "spmat_spmat_product:: Inner dimension of matrices should match\n.");
            return (res);
        }
        cholmod_common chol_c;
        cholmod_start(&chol_c);

        cholmod_sparse *chol_A = as_cholmod_sparse(A, chol_A, &chol_c);
        cholmod_sparse *chol_B = as_cholmod_sparse(B, chol_B, &chol_c);

        cholmod_sparse *chol_res = cholmod_ssmult(chol_A, chol_B, 0, true,
                                                  true, &chol_c);

        res = as_arma_sparse(chol_res, res, &chol_c);

        cholmod_free_sparse(&chol_A, &chol_c);
        cholmod_free_sparse(&chol_B, &chol_c);
        cholmod_finish(&chol_c);

        return (res);
    }

    mat spmat_mat_product_parallel(sp_mat &A, mat &B, int thread_no)
    {
        if (thread_no <= 0)
        {
            thread_no = SYS_THREADS_DEF;
        }
        int M = A.n_rows;
        int N = B.n_cols;
        mat res = zeros(M, N);

        if (thread_no > N)
        {
            thread_no = N;
            ParallelFor(0, B.n_cols, thread_no, [&](unsigned int k, size_t threadId)
                        {
                            vec u = B.col(k);
                            vec v = spmat_vec_product(A, u);
                            mtx.lock();
                            res.col(k) = v;
                            mtx.unlock();
                        });
        }
        else
        {
            int slice_size = ceil((double)N / thread_no);

            ParallelFor(0, thread_no, thread_no, [&](unsigned int k, size_t threadId)
                        {
                            int i = k * slice_size;
                            if (i <= (N - 1))
                            {
                                int j = (k + 1) * slice_size - 1;
                                if (j > (N - 1))
                                    j = N - 1;

                                mat subB = B.cols(i, j);
                                mat subC = spmat_mat_product(A, subB);
                                mtx.lock();
                                res.cols(i, j) = subC;
                                mtx.unlock();
                            }
                        });
        }

        return (res);
    }

    mat mat_mat_product_parallel(mat &A, mat &B, int thread_no)
    {
        if (thread_no <= 0)
        {
            thread_no = SYS_THREADS_DEF;
        }

        int M = A.n_rows;
        int N = B.n_cols;
        mat res = zeros(M, N);

        if (thread_no > N)
        {
            thread_no = N;
            ParallelFor(0, B.n_cols, thread_no, [&](size_t k, size_t threadId)
                        {
                            vec u = B.col(k);
                            vec v = A * u;
                            mtx.lock();
                            res.col(k) = v;
                            mtx.unlock();
                        });
        }
        else
        {
            int slice_size = ceil((double)N / thread_no);

            ParallelFor(0, thread_no, thread_no, [&](size_t k, size_t threadId)
                        {
                            int i = k * slice_size;
                            if (i <= (N - 1))
                            {
                                int j = (k + 1) * slice_size - 1;
                                if (j > (N - 1))
                                    j = N - 1;

                                mat subB = B.cols(i, j);
                                mat subC = A * subB;
                                mtx.lock();
                                res.cols(i, j) = subC;
                                mtx.unlock();
                            }
                        });
        }

        return (res);
    }

}