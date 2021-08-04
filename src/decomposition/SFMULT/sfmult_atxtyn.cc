//==============================================================================
//=== sfmult_AT_XT_YN ==========================================================
//==============================================================================

// y = A'*x'	    A is m-by-n, x is k-by-m, y is n-by-k

// compare with sfmult_AT_XN_YN for kernel usage
// compare with sfmult_AT_XT_YT for outer loop structure but different kernels

#include "sfmult.h"

arma::mat &sfmult_AT_XT_YN // y = A'*x'
    (
        const arma::sp_mat &A,
        const arma::mat &X,
        int ac, // if true: conj(A)   if false: A. ignored if A real
        int xc, // if true: conj(x)   if false: x. ignored if x real
        int yc  // if true: conj(y)   if false: y. ignored if y real
    )
{
    double *Ax, *Az, *Xx, *Xz, *Yx, *Yz, *Wx, *Wz;
    Int *Ap, *Ai;
    Int m, n, k, k1, i;

    //--------------------------------------------------------------------------
    // get inputs
    //--------------------------------------------------------------------------

    m = A.n_rows;
    n = A.n_cols;
    k = X.n_cols;
    static mat Y(m, k);

    if (n != X.n_rows)
    {
        stderr_printf("Number of rows in X and columns in A do not match.");
        return (Y);
    }

    A.sync();

    Ap = arma::access::rwp(A.col_ptrs);
    Ai = arma::access::rwp(A.row_indices);
    Ax = arma::access::rwp(A.values);
    Az = NULL;
    Xx = (double *)X.memptr();
    Xz = NULL;

    //--------------------------------------------------------------------------
    // allocate result
    //--------------------------------------------------------------------------
    Yx = Y.memptr();
    Yz = NULL;

    //--------------------------------------------------------------------------
    // special cases
    //--------------------------------------------------------------------------

    if (k == 0 || m == 0 || n == 0 || Ap[n] == 0)
    {
        // Y = 0
        return (sfmult_yzero(Y));
    }
    if (k == 1)
    {
        // Y = A' * X
        sfmult_AT_x_1(Yx, Yz, Ap, Ai, Ax, Az, m, n, Xx, Xz, ac, xc, yc);
        return (Y);
    }
    if (k == 2)
    {
        // Y = A' * X'
        sfmult_AT_XT_YN_2(Yx, Yz, Ap, Ai, Ax, Az, m, n, Xx, Xz, ac, xc, yc);
        return (Y);
    }
    if (k == 4)
    {
        // Y = A' * X'
        sfmult_AT_XT_YN_4(Yx, Yz, Ap, Ai, Ax, Az, m, n, Xx, Xz, ac, xc, yc);
        return (Y);
    }

    //--------------------------------------------------------------------------
    // allocate workspace
    //--------------------------------------------------------------------------

    sfmult_walloc(4, m, &Wx, &Wz);

    //--------------------------------------------------------------------------
    // Y = A'*X', in blocks of up to 4 columns of X using sfmult_atxtyn
    //--------------------------------------------------------------------------

    k1 = k % 4;
    if (k1 == 1)
    {
        // W = X (1,:)'
        for (i = 0; i < m; i++)
        {
            Wx[i] = Xx[k * i];
        }
        // Y (:,1) = A' * W
        sfmult_AT_x_1(Yx, Yz, Ap, Ai, Ax, Az, m, n, Wx, Wz, ac, xc, yc);
        Yx += n;
        Yz += n;
        Xx += 1;
        Xz += 1;
    }
    else if (k1 == 2)
    {
        // W = X (1:2,:)
        for (i = 0; i < m; i++)
        {
            Wx[2 * i] = Xx[k * i];
            Wx[2 * i + 1] = Xx[k * i + 1];
        }
        // Y = A' * W'
        sfmult_AT_XT_YN_2(Yx, Yz, Ap, Ai, Ax, Az, m, n, Wx, Wz, ac, xc, yc);
        Yx += 2 * n;
        Yz += 2 * n;
        Xx += 2;
        Xz += 2;
    }
    else if (k1 == 3)
    {
        // W = X (1:3,:)
        for (i = 0; i < m; i++)
        {
            Wx[4 * i] = Xx[k * i];
            Wx[4 * i + 1] = Xx[k * i + 1];
            Wx[4 * i + 2] = Xx[k * i + 2];
        }
        // Y = A' * W'
        sfmult_AT_XT_YN_3(Yx, Yz, Ap, Ai, Ax, Az, m, n, Wx, Wz, ac, xc, yc);
        Yx += 3 * n;
        Yz += 3 * n;
        Xx += 3;
        Xz += 3;
    }
    for (; k1 < k; k1 += 4)
    {
        // W = X (1:4,:)
        for (i = 0; i < m; i++)
        {
            Wx[4 * i] = Xx[k * i];
            Wx[4 * i + 1] = Xx[k * i + 1];
            Wx[4 * i + 2] = Xx[k * i + 2];
            Wx[4 * i + 3] = Xx[k * i + 3];
        }
        // Y (:,k1+(1:4)) = A' * W'
        sfmult_AT_XT_YN_4(Yx, Yz, Ap, Ai, Ax, Az, m, n, Wx, Wz, ac, xc, yc);
        Yx += 4 * n;
        Yz += 4 * n;
        Xx += 4;
        Xz += 4;
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    free(Wx);
    return (Y);
}
