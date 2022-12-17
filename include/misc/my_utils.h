#ifndef MY_UTILS_H
#define MY_UTILS_H

#include <ACTIONet.h>
#include <pcg_random.hpp>

namespace ACTIONet
{
    mat sampleUnif(int l, int m, double a, double b, int seed);
    void gram_schmidt(mat &A);
    field<mat> eigSVD(mat A);
    mat randNorm(int l, int m, int seed);

    mat zscore(mat &A, int dim = 0, int thread_no = 1);
    mat tzscoret(mat &A);
    mat robust_zscore(mat &A, int dim = 0, int thread_no = 1);
    mat RIN_transform(mat &A, int thread_no = 1);

    // Used in IRLB
    void randNorm_inplace(int n, double *out, int seed);
    void convtests(int Bsz, int n, double tol, double svtol, double Smax,
                   double *svratio, double *residuals, int *k, int *converged,
                   double S);
    void orthog(double *X, double *Y, double *T, int xm, int xn, int yn);

    uint32_t lfsr113(uint64_t **state);
    void lfsr113_seed(uint32_t seed, uint64_t **state);

    void randN_Marsaglia(double *values, int n, pcg32 rng);
    void randN_BM(double *values, int n, uint64_t **state);
    void randN_normsinv(double *values, int n);

    vec spmat_vec_product(sp_mat &A, vec &x);
    mat spmat_mat_product(sp_mat &A, mat &B);
    sp_mat spmat_spmat_product(sp_mat &A, sp_mat &B);

    mat spmat_mat_product_parallel(sp_mat &A, mat &B, int thread_no);
    mat mat_mat_product_parallel(mat &A, mat &B, int thread_no);

} // namespace ACTIONet

#endif
