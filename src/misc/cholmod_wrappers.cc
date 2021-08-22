#include "ACTIONet.h"
#include <cholmod.h>

namespace ACTIONet
{
    cholmod_sparse *as_cholmod_sparse(cholmod_sparse *ans,
                                      sp_mat &A)
    {
        ans->itype = CHOLMOD_INT; /* characteristics of the system */
        ans->dtype = CHOLMOD_DOUBLE;
        ans->packed = true;

        A.sync();
        /* DOES NOT WORK!!! arma uses (arma::uword *) whereas CHOLMOD uses (int *)
        void *x_ptr0 = (void *)const_cast<double *>(reinterpret_cast<const double *>(A.values));
        void *i_ptr0 = (void *)const_cast<arma::uword *>(reinterpret_cast<const arma::uword *>(A.row_indices));
        void *p_ptr0 = (void *)const_cast<arma::uword *>(reinterpret_cast<const arma::uword *>(A.col_ptrs));

        ans->x = x_ptr0;
        ans->i = i_ptr0;
        ans->p = p_ptr0;
    
        // DOES NOT WORK! Same as above
        std::vector<double> x(A.values, A.values + A.n_nonzero);
        std::vector<int> i(A.row_indices, A.row_indices + A.n_nonzero);
        std::vector<int> p(A.col_ptrs, A.col_ptrs + A.n_cols + 1);

        ans->x = x.data();
        ans->i = i.data();
        ans->p = p.data();
        */

        ans->i = new int[A.n_nonzero];
        ans->x = new double[A.n_nonzero];
        double *x_ptr = (double *)ans->x;
        int *i_ptr = (int *)ans->i;
        for (int k = 0; k < A.n_nonzero; k++)
        {
            x_ptr[k] = (A.values)[k];
            i_ptr[k] = (A.row_indices)[k];
        }

        ans->p = new int[(A.n_cols + 1)];
        int *ptr = (int *)ans->p;
        for (int k = 0; k < A.n_cols + 1; k++)
        {
            ptr[k] = A.col_ptrs[k];
        }

        ans->nrow = A.n_rows;
        ans->ncol = A.n_cols;
        ans->nzmax = A.n_nonzero;

        ans->xtype = CHOLMOD_REAL;
        ans->stype = 0;
        ans->dtype = 0;

        ans->sorted = 1;

        return ans;
    }

    // Mat-vec product
    //dsdmult('n', A.n_rows, A.n_cols, A_as_cholmod, x.memptr(), out.memptr(), &chol_c);
    void dsdmult(char transpose, int n_rows, int n_cols, void *A, double *x, double *out,
                 cholmod_common *chol_cp)
    {
        int t = transpose == 'n' ? 0 : 1; // 'n': computes Ax, 't': computes A'x
        cholmod_sparse *cha = (cholmod_sparse *)A;

        cholmod_dense chb;
        chb.nrow = transpose == 'n' ? n_rows : n_cols;
        chb.d = chb.nrow;
        chb.ncol = 1;
        chb.nzmax = chb.nrow;
        chb.xtype = cha->xtype;
        chb.dtype = 0;
        chb.x = (void *)x;
        chb.z = (void *)NULL;

        cholmod_dense chc;
        chc.nrow = transpose == 'n' ? n_cols : n_rows;
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

        cholmod_sparse *As = new cholmod_sparse[1];
        as_cholmod_sparse(As, A);

        vec Ax = zeros(A.n_rows);
        dsdmult('n', A.n_rows, A.n_cols, As, x.memptr(), Ax.memptr(), &chol_c);

        cholmod_free_sparse(&As, &chol_c);
        cholmod_finish(&chol_c);

        return (Ax);
    }

    mat spmat_mat_product(sp_mat &A, mat &B)
    {
        cholmod_common chol_c;
        cholmod_start(&chol_c);
        chol_c.final_ll = 1;

        cholmod_sparse *As = new cholmod_sparse[1];
        as_cholmod_sparse(As, A);

        cholmod_dense *Bd = new cholmod_dense[1];
        Bd->nrow = B.n_rows;
        Bd->ncol = B.n_cols;
        Bd->d = B.n_rows;
        Bd->nzmax = B.n_rows * B.n_cols;
        Bd->xtype = CHOLMOD_REAL;
        Bd->dtype = 0;
        Bd->x = (void *)B.memptr();
        Bd->z = (void *)NULL;

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
        cholmod_sdmult(As, 0, one, zero, Bd, &out, &chol_c);

        cholmod_free_sparse(&As, &chol_c);
        //cholmod_free_dense(&Bd, &chol_c);
        cholmod_finish(&chol_c);

        return (res);
    }

    sp_mat spmat_spmat_product(sp_mat &A, sp_mat &B)
    {
        cholmod_common chol_c;
        cholmod_start(&chol_c);
        chol_c.final_ll = 1;

        cholmod_sparse *As = new cholmod_sparse[1];
        as_cholmod_sparse(As, A);

        cholmod_sparse *Bs = new cholmod_sparse[1];
        as_cholmod_sparse(Bs, B);

        int nrow = A.n_rows, ncol = B.n_cols;
        arma::sp_mat res(A.n_rows, B.n_cols);

        cholmod_sparse *ans = cholmod_ssmult(As, Bs, 0, true,
                                             true, &chol_c);

        res.mem_resize(static_cast<unsigned>(ans->nzmax));
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

        cholmod_free_sparse(&As, &chol_c);
        cholmod_free_sparse(&Bs, &chol_c);

        cholmod_finish(&chol_c);

        return (res);
    }
}