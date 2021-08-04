//==============================================================================
//=== sfmult_AN_XT_YN ==========================================================
//==============================================================================

// y = A*x'	    A is m-by-n, x is k-by-n, y is m-by-k

// compare with sfmult_AN_XT_YT for kernel usage
// compare with sfmult_AN_XN_YN for outer loop structure but different kernels

#include "sfmult.h"

arma::mat &sfmult_AN_XT_YN // y = A*x'
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
        // Y = A*X
        sfmult_AN_x_1(Yx, Yz, Ap, Ai, Ax, Az, m, n, Xx, Xz, ac, xc, yc);
        return (Y);
    }

    //--------------------------------------------------------------------------
    // allocate workspace
    //--------------------------------------------------------------------------

    sfmult_walloc((k == 2) ? 2 : 4, m, &Wx, &Wz);

    //--------------------------------------------------------------------------
    // Y = A*X', in blocks of up to 4 columns of X, using sfmult_anxtyt
    //--------------------------------------------------------------------------

    k1 = k % 4;
    if (k1 == 1)
    {
        // Y (:,1) = A * X(1,:)'
        sfmult_AN_xk_1(Yx, Yz, Ap, Ai, Ax, Az, m, n, Xx, Xz, ac, xc, yc, k);
        Yx += m;
        Yz += m;
        Xx += 1;
        Xz += 1;
    }
    else if (k1 == 2)
    {
        // W = (A * X(1:2,:)')'
        sfmult_AN_XT_YT_2(Wx, Wz, Ap, Ai, Ax, Az, m, n, Xx, Xz, ac, xc, yc, k);
        // Y (:,1:2) = W'

        for (i = 0; i < m; i++)
            Yx[i] = Wx[2 * i];
        Yx += m;
        Yz += m;
        for (i = 0; i < m; i++)
            Yx[i] = Wx[2 * i + 1];
        Yx += m;
        Yz += m;

#if 0
	for (i = 0 ; i < m ; i++)
	{
	    Yx [i  ] = Wx [2*i  ] ;
	    Yx [i+m] = Wx [2*i+1] ;
	}
	Yx += 2*m ;
	Yz += 2*m ;
#endif

        Xx += 2;
        Xz += 2;
    }
    else if (k1 == 3)
    {
        // W = (A * X(1:3,:)')'
        sfmult_AN_XT_YT_3(Wx, Wz, Ap, Ai, Ax, Az, m, n, Xx, Xz, ac, xc, yc, k);
        // Y (:,1:3) = W'

        for (i = 0; i < m; i++)
            Yx[i] = Wx[4 * i];
        Yx += m;
        Yz += m;
        for (i = 0; i < m; i++)
            Yx[i] = Wx[4 * i + 1];
        Yx += m;
        Yz += m;
        for (i = 0; i < m; i++)
            Yx[i] = Wx[4 * i + 2];
        Yx += m;
        Yz += m;

#if 0
	for (i = 0 ; i < m ; i++)
	{
	    Yx [i    ] = Wx [4*i  ] ;
	    Yx [i+  m] = Wx [4*i+1] ;
	    Yx [i+2*m] = Wx [4*i+2] ;
	}
	Yx += 3*m ;
	Yz += 3*m ;
#endif

        Xx += 3;
        Xz += 3;
    }
    for (; k1 < k; k1 += 4)
    {
        // W = (A * X(1:4,:)')'
        sfmult_AN_XT_YT_4(Wx, Wz, Ap, Ai, Ax, Az, m, n, Xx, Xz, ac, xc, yc, k);
        // Y (:,k1+(1:4)) = W'

        for (i = 0; i < m; i++)
            Yx[i] = Wx[4 * i];
        Yx += m;
        Yz += m;
        for (i = 0; i < m; i++)
            Yx[i] = Wx[4 * i + 1];
        Yx += m;
        Yz += m;
        for (i = 0; i < m; i++)
            Yx[i] = Wx[4 * i + 2];
        Yx += m;
        Yz += m;
        for (i = 0; i < m; i++)
            Yx[i] = Wx[4 * i + 3];
        Yx += m;
        Yz += m;

#if 0
	for (i = 0 ; i < m ; i++)
	{
	    Yx [i    ] = Wx [4*i  ] ;
	    Yx [i+  m] = Wx [4*i+1] ;
	    Yx [i+2*m] = Wx [4*i+2] ;
	    Yx [i+3*m] = Wx [4*i+3] ;
	}
	Yx += 4*m ;
	Yz += 4*m ;
#endif

        Xx += 4;
        Xz += 4;
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    free(Wx);
    return (Y);
}
