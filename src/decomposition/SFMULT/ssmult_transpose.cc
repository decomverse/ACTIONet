//==============================================================================
// ssmult_transpose
//==============================================================================

// C = A' or A.' where the input matrix A may have unsorted columns.  The output
// C is always returned with sorted columns.

#include "sfmult.h"

arma::sp_mat &ssmult_transpose // returns C = A' or A.'
    (
        const arma::sp_mat &A,
        int conj // compute A' if true, compute A.' if false
    )
{
    Int *Cp, *Ci, *Ap, *Ai, *W;
    double *Cx, *Cz, *Ax, *Az; // (TO DO): do single too
    Int p, pend, q, i, j, n, m, anz, cnz;

    //--------------------------------------------------------------------------
    // get inputs
    //--------------------------------------------------------------------------

    m = A.n_rows;
    n = A.n_cols;
    A.sync();
    Ap = arma::access::rwp(A.col_ptrs);
    Ai = arma::access::rwp(A.row_indices);
    Ax = arma::access::rwp(A.values);
    Az = NULL;

    anz = Ap[n];

    //--------------------------------------------------------------------------
    // allocate C but do not initialize it
    //--------------------------------------------------------------------------

    cnz = MAX(anz, 1);

    Cp = (Int *)malloc((m + 1) * sizeof(Int));
    Ci = (Int *)malloc(MAX(cnz, 1) * sizeof(Int));
    Cx = (double *)malloc(MAX(cnz, 1) * sizeof(double));

    //--------------------------------------------------------------------------
    // allocate workspace
    //--------------------------------------------------------------------------

    W = (Int *)malloc(MAX(m, 1) * sizeof(Int));

    //--------------------------------------------------------------------------
    // compute row counts
    //--------------------------------------------------------------------------

    for (p = 0; p < anz; p++)
    {
        W[Ai[p]]++;
    }

    //--------------------------------------------------------------------------
    // compute column pointers of C and copy back into W
    //--------------------------------------------------------------------------

    for (p = 0, i = 0; i < m; i++)
    {
        Cp[i] = p;
        p += W[i];
        W[i] = Cp[i];
    }
    Cp[m] = p;

    //--------------------------------------------------------------------------
    // C = A'
    //--------------------------------------------------------------------------

    p = 0;
    // C = A' (real case)
    for (j = 0; j < n; j++)
    {
        pend = Ap[j + 1];
        for (; p < pend; p++)
        {
            q = W[Ai[p]]++; // find position for C(j,i)
            Ci[q] = j;      // place A(i,j) as entry C(j,i)
            Cx[q] = Ax[p];
        }
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    free(W);

    static sp_mat C(n, m);
    C.sync();

    // Making space for the elements
    C.mem_resize(static_cast<unsigned>(cnz));

    // Copying elements
    memcpy(arma::access::rwp(C.col_ptrs), Cp, (m + 1) * sizeof(Int));
    memcpy(arma::access::rwp(C.row_indices), Ci, MAX(cnz, 1) * sizeof(Int));
    memcpy(arma::access::rwp(C.values), Cx, MAX(cnz, 1) * sizeof(double));

    return (C);
}
