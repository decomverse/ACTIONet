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
    cholmod_sparse *as_cholmod_sparse(cholmod_sparse *ans,
                                      sp_mat &A)
    {
        ans->itype = CHOLMOD_INT; /* characteristics of the system */
        ans->dtype = CHOLMOD_DOUBLE;
        ans->packed = true;

        ans->i = new int[A.n_nonzero];
        ans->x = new double[A.n_nonzero];
        ans->p = new int[(A.n_cols + 2)];
        int *ptr = (int *)ans->p;
        double *x_ptr = (double *)ans->x;
        int *i_ptr = (int *)ans->i;

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

        ans->nrow = A.n_rows;
        ans->ncol = A.n_cols;
        ans->nzmax = A.n_nonzero;

        ans->xtype = CHOLMOD_REAL;
        ans->stype = 0;
        ans->dtype = 0;

        ans->sorted = 1;

        return ans;
    }

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

    vec spmat_vec_product(sp_mat &A, vec &x)
    {
        cholmod_common chol_c;
        cholmod_start(&chol_c);
        chol_c.final_ll = 1;

        cholmod_sparse As;
        as_cholmod_sparse(&As, A);
        vec Ax = zeros(A.n_rows);
        dsdmult('n', A.n_rows, A.n_cols, &As, x.memptr(), Ax.memptr(), &chol_c);

        cholmod_finish(&chol_c);

        return (Ax);
    }

    mat spmat_mat_product(sp_mat &A, mat &B)
    {
        cholmod_common chol_c;
        cholmod_start(&chol_c);
        chol_c.final_ll = 1;

        cholmod_sparse As;
        as_cholmod_sparse(&As, A);

        cholmod_dense Bd;
        Bd.nrow = B.n_rows;
        Bd.ncol = B.n_cols;
        Bd.d = B.n_rows;
        Bd.nzmax = B.n_rows * B.n_cols;
        Bd.xtype = CHOLMOD_REAL;
        Bd.dtype = 0;
        Bd.x = (void *)B.memptr();
        Bd.z = (void *)NULL;

        mat res = zeros(A.n_rows, B.n_cols);
        cholmod_dense out;
        out.nrow = A.n_rows;
        out.ncol = B.n_cols;
        out.d = A.n_rows;
        out.nzmax = A.n_rows * B.n_cols;
        out.xtype = CHOLMOD_REAL;
        out.dtype = 0;
        out.x = (void *)res.memptr();
        out.z = (void *)NULL;

        double one[] = {1, 0}, zero[] = {0, 0};
        cholmod_sdmult(&As, 0, one, zero, &Bd, &out, &chol_c);

        cholmod_finish(&chol_c);
        return (res);
    }

    sp_mat spmat_spmat_product(sp_mat &A, sp_mat &B)
    {

        cholmod_common chol_c;
        cholmod_start(&chol_c);
        chol_c.final_ll = 1;

        cholmod_sparse As, Bs;
        as_cholmod_sparse(&As, A);
        as_cholmod_sparse(&Bs, B);

        int nrow = A.n_rows, ncol = B.n_cols;

        cholmod_sparse *ans = cholmod_ssmult(&As, &Bs, 0, true,
                                             true, &chol_c);

        arma::sp_mat res(A.n_rows, B.n_cols);
        res.mem_resize(static_cast<unsigned>(ans->nzmax));

        mtx.lock();
        res.sync();
        double *res_x_ptr = (double *)ans->x;
        int *res_i_ptr = (int *)ans->i;
        double *out_x_ptr = (double *)arma::access::rwp(res.values);
        arma::uword *out_i_ptr = arma::access::rwp(res.row_indices);
        for (int k = 0; k < ans->nzmax; k++)
        {
            out_x_ptr[k] = res_x_ptr[k];
            out_i_ptr[k] = res_i_ptr[k];
        }

        int *res_p_ptr = (int *)ans->p;
        arma::uword *out_p_ptr = arma::access::rwp(res.col_ptrs);
        for (int k = 0; k < B.n_cols; k++)
        {
            out_p_ptr[k] = res_p_ptr[k];
        }

        // important: set the sentinel as well
        arma::access::rwp(res.col_ptrs)[B.n_cols] = ans->nzmax;

        // set the number of non-zero elements
        arma::access::rw(res.n_nonzero) = ans->nzmax;
        mtx.unlock();

        cholmod_finish(&chol_c);

        return (res);
    }

    mat spmat_mat_product_parallel(sp_mat &A, mat &B, int thread_no)
    {

        int M = A.n_rows;
        int N = B.n_cols;
        mat res = zeros(M, N);

        if (thread_no > N)
        {
            thread_no = N;
            ParallelFor(0, B.n_cols, thread_no, [&](size_t k, size_t threadId)
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
            double *out = (double *)malloc(M * N * sizeof(double));

            ParallelFor(0, thread_no, thread_no, [&](size_t k, size_t threadId)
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
            double *out = (double *)malloc(M * N * sizeof(double));

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