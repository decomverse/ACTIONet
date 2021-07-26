#ifndef MY_UTILS_H
#define MY_UTILS_H

#include <math/pcg/pcg_random.hpp>
#include <cholmod.h>

namespace ACTIONet {
mat sampleUnif(int l, int m, double a, double b, int seed);
void gram_schmidt(mat &A);
field<mat> eigSVD(mat A);
mat randNorm(int l, int m, int seed);

mat zscore(mat A);
mat robust_zscore(mat A);

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

void dsdmult(char transpose, int m, int n, void *a, double *b, double *c,
             cholmod_common *chol_cp);
cholmod_sparse_struct *as_cholmod_sparse(cholmod_sparse_struct *ans, sp_mat &A);

}  // namespace ACTIONet

#endif
