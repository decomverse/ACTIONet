#ifndef CBLAS_ALT_TEMPLATE_H
#define CBLAS_ALT_TEMPLATE_H

#ifndef CBLAS_ALT_TEMPLATE
#define CBLAS_ALCBLAS_ALTT_TEMPLATE


//#include <cblas.h>
#include <stddef.h>
//#include <blas.h>
#ifdef small
#undef small
#endif
//#include <lapack.h>
#include <cblas_defvar.h>


#define NEW_MATLAB_BLAS

#ifdef NEW_MATLAB_BLAS
#define INTT ptrdiff_t
#else
#ifdef OLD_MATLAB_BLAS
#define INTT int
#else
#ifdef MKL_INT
#define INTT MKL_INT
#else
#ifdef INT_64BITS
#define INTT long long int
#else
#define INTT int
#endif
#endif
#endif
#endif

#ifdef INT_64BITS
#define INTM INTT
#else
#define INTM int
#endif


/// a few static variables for lapack
static char lower='L';
static char upper='u';
static INTT info=0;
static char incr='I';
static char decr='D';
static char no='N';
static char reduced='S';
static char allV='V';

#ifdef ADD_ATL
#define dnrm2_ ATL_dnrm2
#define snrm2_ ATL_snrm2
#define dcopy_ ATL_dcopy
#define scopy_ ATL_scopy
#define daxpy_ ATL_daxpy
#define saxpy_ ATL_saxpy
#define daxpby_ ATL_daxpby
#define saxpby_ ATL_saxpby
#define dscal_ ATL_dscal
#define sscal_ ATL_sscal
#define dasum_ ATL_dasum
#define sasum_ ATL_sasum
#define ddot_ ATL_ddot
#define sdot_ ATL_sdot
#define dgemv_ ATL_dgemv
#define sgemv_ ATL_sgemv
#define dger_ ATL_dger
#define sger_ ATL_sger
#define dtrmv_ ATL_dtrmv
#define strmv_ ATL_strmv
#define dsyr_ ATL_dsyr
#define ssyr_ ATL_ssyr
#define dsymv_ ATL_dsymv
#define ssymv_ ATL_ssymv
#define dgemm_ ATL_dgemm
#define sgemm_ ATL_sgemm
#define dsyrk_ ATL_dsyrk
#define ssyrk_ ATL_ssyrk
#define dtrmm_ ATL_dtrmm
#define strmm_ ATL_strmm
#define dtrtri_ ATL_dtrtri
#define strtri_ ATL_strtri
#define idamax_ ATL_idamax
#define isamax_ ATL_isamax
#endif

#ifdef ACCELERATE
	#define REMOVE_
#endif

#ifdef REMOVE_
#define dnrm2_ dnrm2
#define snrm2_ snrm2
#define dcopy_ dcopy
#define scopy_ scopy
#define daxpy_ daxpy
#define saxpy_ saxpy
#define daxpby_ daxpby
#define saxpby_ saxpby
#define dscal_ dscal
#define sscal_ sscal
#define dasum_ dasum
#define sasum_ sasum
#define ddot_ ddot
#define sdot_ sdot
#define ddoti_ ddoti
#define sdoti_ sdoti
#define dgemv_ dgemv
#define sgemv_ sgemv
#define dger_ dger
#define sger_ sger
#define dtrmv_ dtrmv
#define strmv_ strmv
#define dsyr_ dsyr
#define ssyr_ ssyr
#define dsymv_ dsymv
#define ssymv_ ssymv
#define dgemm_ dgemm
#define sgemm_ sgemm
#define dsyrk_ dsyrk
#define ssyrk_ ssyrk
#define dtrmm_ dtrmm
#define strmm_ strmm
#define dtrtri_ dtrtri
#define strtri_ strtri
#define idamax_ idamax
#define isamax_ isamax
#define dsytrf_ dsytrf
#define ssytrf_ ssytrf
#define dsytri_ dsytri
#define ssytri_ ssytri
#define dlasrt_ dlasrt
#define slasrt_ slasrt
#define dgesvd_ dgesvd
#define sgesvd_ sgesvd
#define dsyev_ dsyev
#define ssyev_ ssyev
#endif

/// external functions
//#ifdef HAVE_MKL   // obsolete
//extern "C" {
//#endif
INTT cblas_idamin( INTT n,  double* X,  INTT incX);
INTT cblas_isamin( INTT n,  float* X,  INTT incX);
//#ifdef HAVE_MKL
//};
//#endif

//#define HAVE_MKL
/*#ifdef HAVE_MKL   // obsolete, do not use

#define idamin_ idamin
#define isamin_ isamin
extern "C" {
   void vdSqr( int n,  double* vecIn, double* vecOut);
   void vsSqr( int n,  float* vecIn, float* vecOut);
   void vdSqrt( int n,  double* vecIn, double* vecOut);
   void vsSqrt( int n,  float* vecIn, float* vecOut);
   void vdInvSqrt( int n,  double* vecIn, double* vecOut);
   void vsInvSqrt( int n,  float* vecIn, float* vecOut);
   void vdSub( int n,  double* vecIn,  double* vecIn2, double* vecOut);
   void vsSub( int n,  float* vecIn,  float* vecIn2, float* vecOut);
   void vdDiv( int n,  double* vecIn,  double* vecIn2, double* vecOut);
   void vsDiv( int n,  float* vecIn,  float* vecIn2, float* vecOut);
   void vdExp( int n,  double* vecIn, double* vecOut);
   void vsExp( int n,  float* vecIn, float* vecOut);
   void vdInv( int n,  double* vecIn, double* vecOut);
   void vsInv( int n,  float* vecIn, float* vecOut);
   void vdAdd( int n,  double* vecIn,  double* vecIn2, double* vecOut);
   void vsAdd( int n,  float* vecIn,  float* vecIn2, float* vecOut);
   void vdMul( int n,  double* vecIn,  double* vecIn2, double* vecOut);
   void vsMul( int n,  float* vecIn,  float* vecIn2, float* vecOut);
   void vdAbs( int n,  double* vecIn, double* vecOut);
   void vsAbs( int n,  float* vecIn, float* vecOut);
}
#endif*/

 // INTTerfaces to a few BLAS function, Level 1
 /// INTTerface to cblas_*nrm2
 template <typename T> T inline cblas_nrm2( INTT n,  T* X,  INTT incX);
 /// INTTerface to cblas_*copy
 template <typename T> void inline cblas_copy( INTT n,  T* X,  INTT incX,
                                               T* Y,  INTT incY);
 /// INTTerface to cblas_*axpy
 template <typename T> void inline cblas_axpy( INTT n,  T a,  T* X,
                                               INTT incX, T* Y,  INTT incY);
 template <typename T> void inline cblas_axpby( INTT n,  T a,  T* X,
                                                INTT incX, T b,  T* Y,  INTT incY);
 /// INTTerface to cblas_*scal
 template <typename T> void inline cblas_scal( INTT n,  T a, T* X,
                                               INTT incX);
 /// INTTerface to cblas_*asum
 template <typename T> T inline cblas_asum( INTT n,  T* X,  INTT incX);
 /// INTTerface to cblas_*adot
 template <typename T> T inline cblas_dot( INTT n,  T* X,  INTT incX,
                                           T* Y, INTT incY);
 /// interface to cblas_i*amin
 template <typename T> INTT inline cblas_iamin( INTT n,  T* X,  INTT incX);
 /// interface to cblas_i*amax
 template <typename T> INTT inline cblas_iamax( INTT n,  T* X,  INTT incX);
 
 // INTTerfaces to a few BLAS function, Level 2
 
 /// INTTerface to cblas_*gemv
 template <typename T> void cblas_gemv( CBLAS_ORDER order,
                                        CBLAS_TRANSPOSE TransA,  INTT M,
                                        INTT N,  T alpha,  T *A,  INTT lda,  T *X,
                                        INTT incX,  T beta,T *Y,   INTT incY);
 /// INTTerface to cblas_*trmv
 template <typename T> void inline cblas_trmv( CBLAS_ORDER order,  CBLAS_UPLO Uplo,
                                               CBLAS_TRANSPOSE TransA,  CBLAS_DIAG Diag,  INTT N,
                                               T *A,  INTT lda, T *X,  INTT incX);
 /// INTTerface to cblas_*syr
 template <typename T> void inline cblas_syr( CBLAS_ORDER order,
                                              CBLAS_UPLO Uplo,  INTT N,  T alpha,
                                              T *X,  INTT incX, T *A,  INTT lda);
 
 /// INTTerface to cblas_*symv
 template <typename T> inline void cblas_symv( CBLAS_ORDER order,
                                               CBLAS_UPLO Uplo,  INTT N,
                                               T alpha,  T *A,  INTT lda,  T *X,
                                               INTT incX,  T beta,T *Y,   INTT incY);
 
 
 // INTTerfaces to a few BLAS function, Level 3
 /// INTTerface to cblas_*gemm
 template <typename T> void cblas_gemm( CBLAS_ORDER order,
                                        CBLAS_TRANSPOSE TransA,  CBLAS_TRANSPOSE TransB,
                                        INTT M,  INTT N,  INTT K,  T alpha,
                                        T *A,  INTT lda,  T *B,  INTT ldb,
                                        T beta, T *C,  INTT ldc);
 /// INTTerface to cblas_*syrk
 template <typename T> void cblas_syrk( CBLAS_ORDER order,
                                        CBLAS_UPLO Uplo,  CBLAS_TRANSPOSE Trans,  INTT N,  INTT K,
                                        T alpha,  T *A,  INTT lda,
                                        T beta, T*C,  INTT ldc);
 /// INTTerface to cblas_*ger
 template <typename T> void cblas_ger( CBLAS_ORDER order,
                                       INTT M,  INTT N,  T alpha,  T *X,  INTT incX,
                                       T* Y,  INTT incY, T*A,  INTT lda);
 /// INTTerface to cblas_*trmm
 template <typename T> void cblas_trmm( CBLAS_ORDER order,
                                        CBLAS_SIDE Side,  CBLAS_UPLO Uplo,
                                        CBLAS_TRANSPOSE TransA,  CBLAS_DIAG Diag,
                                        INTT M,  INTT N,  T alpha,
                                        T*A,  INTT lda,T *B,  INTT ldb);
 
 // interfaces to a few functions from the intel Vector Mathematical Library
 /// INTTerface to v*Sqr
 template <typename T> void vSqrt( INTT n,  T* vecIn, T* vecOut);
 /// INTTerface to v*Sqr
 template <typename T> void vInvSqrt( INTT n,  T* vecIn, T* vecOut);
 /// INTTerface to v*Sqr
 template <typename T> void vSqr( INTT n,  T* vecIn, T* vecOut);
 /// INTTerface to v*Sub
 template <typename T> void vSub( INTT n,  T* vecIn,  T* vecIn2, T* vecOut);
 /// INTTerface to v*Div
 template <typename T> void vDiv( INTT n,  T* vecIn,  T* vecIn2, T* vecOut);
 /// INTTerface to v*Exp
 template <typename T> void vExp( INTT n,  T* vecIn, T* vecOut);
 /// INTTerface to v*Inv
 template <typename T> void vInv( INTT n,  T* vecIn, T* vecOut);
 /// INTTerface to v*Add
 template <typename T> void vAdd( INTT n,  T* vecIn,  T* vecIn2, T* vecOut);
 /// INTTerface to v*Mul
 template <typename T> void vMul( INTT n,  T* vecIn,  T* vecIn2, T* vecOut);
 /// INTTerface to v*Abs
 template <typename T> void vAbs( INTT n,  T* vecIn, T* vecOut);
 
 // interfaces to a few LAPACK functions
 /// interface to *trtri
 template <typename T> void trtri(char& uplo, char& diag,
                                  INTT n, T * a, INTT lda);
 /// interface to *sytri  // call sytrf
 template <typename T> void sytri(char& uplo, INTT n, T* a, INTT lda);
 //, INTT* ipiv,
 //      T* work);
 /// interaface to *lasrt
 template <typename T> void lasrt(char& id,  INTT n, T *d);
 //template <typename T> void lasrt2(char& id,  INTT& n, T *d, int* key);
 template <typename T> void gesvd( char& jobu, char& jobvt, INTT m,
                                   INTT n, T* a, INTT lda, T* s,
                                   T* u, INTT ldu, T* vt, INTT ldvt);
 template <typename T> void syev( char& jobz, char& uplo, INTT n,
                                  T* a, INTT lda, T* w);

#include <cblas_alt_template.tpp>
 
#endif
#endif // CBLAS_ALT_TEMPLATE_H
