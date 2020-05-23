#ifndef MINI_CBLAS_H
#define MINI_CBLAS_H

	#include <arma_base.h>
	#include "cblas.h"

	#ifndef BLAS_extern
	#define BLAS_extern extern
	#endif

	#define INTT arma::blas_int

	/*
	 * Enumerated and derived types
	 */
	#define CBLAS_INDEX size_t  /* this may vary between platforms */

/*
	enum CBLAS_LAYOUT {CblasRowMajor=101, CblasColMajor=102};
	enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113};
	enum CBLAS_UPLO {CblasUpper=121, CblasLower=122};
	enum CBLAS_DIAG {CblasNonUnit=131, CblasUnit=132};
	enum CBLAS_SIDE {CblasLeft=141, CblasRight=142};


	char CBLAS_TRANSPOSE_CHAR[] = {'N', 'T', 'C'};
	char *cblas_transpose(CBLAS_TRANSPOSE TransA)
	{
		switch(TransA)
		{
			case 111:	return &CBLAS_TRANSPOSE_CHAR[0];
			case 112:	return &CBLAS_TRANSPOSE_CHAR[1];
			case 113:	return &CBLAS_TRANSPOSE_CHAR[2];
		}
		return NULL;
	}

	char CBLAS_UPLO_CHAR[] = {'U', 'L'};
	char *cblas_uplo(CBLAS_UPLO Uplo)
	{
		switch(Uplo)
		{
			case 121:	return &CBLAS_UPLO_CHAR[0];
			case 122:	return &CBLAS_UPLO_CHAR[1];
		}
		return NULL;
	}

	char CBLAS_DIAG_CHAR[] = {'N', 'U'};
	char *cblas_diag(CBLAS_DIAG Diag)
	{
		switch(Diag)
		{
			case 131:	return &CBLAS_DIAG_CHAR[0];
			case 132:	return &CBLAS_DIAG_CHAR[1];
		}
		return NULL;
	}

	char CBLAS_SIDE_CHAR[] = {'L', 'R'};
	char *cblas_side(CBLAS_SIDE Side)
	{
		switch(Side)
		{
			case 141:	return &CBLAS_SIDE_CHAR[0];
			case 142:	return &CBLAS_SIDE_CHAR[1];
		}
		return NULL;
	}
*/

/*
	#if !defined(ARMA_BLAS_CAPITALS)  
	  #define dscal  dscal
	  #define dcopy  dcopy
	  #define daxpy  daxpy
	  #define ssymv  ssymv
	  #define dsymv  dsymv
	  #define dger  dger  
	  #define ddot  ddot  
	  #define dgemv dgemv
	  
	#else
	  #define dscal  DSCAL
	  #define dcopy  DCOPY
	  #define daxpy  DAXPY
	  #define ssymv  SSYMV
	  #define dsymv  DSYMV
	  #define dger  DGER  
	  #define ddot  DDOT
	  #define dgemv DGEMV
	  
	#endif
*/

	#ifdef  __cplusplus
	extern "C" {
	#endif
		BLAS_extern void arma_fortran(dscal)(INTT *n,double* a, double *x,INTT *incX);
		BLAS_extern void arma_fortran(dcopy)(INTT *n,double *x,INTT *incX, double *y,INTT *incY);
		BLAS_extern void arma_fortran(daxpy)(INTT *n,double* a, double *x,INTT *incX, double *y,INTT *incY);
		BLAS_extern void arma_fortran(ssymv)(char   *uplo, INTT *n, float *alpha, float *a, INTT *lda, float *x, INTT *incx, float *beta, float *y, INTT *incy);
		BLAS_extern void arma_fortran(dsymv)(char   *uplo, INTT *n, double *alpha, double *a, INTT *lda, double *x, INTT *incx, double *beta, double *y, INTT *incy);
		BLAS_extern void arma_fortran(dger)(INTT *m, INTT *n, double *alpha, double *x, INTT *incx, double *y, INTT *incy, double *a, INTT *lda);	
		//BLAS_extern double arma_fortran(ddot)(const INTT *n, const double *x, const INTT *incx, const double *y, const INTT *incy);
		//BLAS_extern void arma_fortran(dgemv)(const char *trans, const INTT *m, const INTT *n, const double *alpha,const double *a, const INTT *lda, const double *x, const INTT *incx, const double *beta, double *y, const INTT *incy);
	#ifdef  __cplusplus
	}
	#endif


	inline double cblas_dot( INTT n,  double* X,
		   INTT incX,  double* Y, INTT incY) {
	   return cblas_ddot(n,X,incX,Y,incY);
	   //return arma_fortran(ddot)(&n,X,&incX,Y,&incY);
	};

	inline void cblas_copy( INTT n,  double* X, 
		   INTT incX, double* Y,  INTT incY) {
	   cblas_dcopy(n,X,incX,Y,incY);
	   //arma_fortran(dcopy)(&n,X,&incX,Y,&incY);
	};

	inline void cblas_symv( CBLAS_LAYOUT order,
		   CBLAS_UPLO Uplo,  INTT N, 
		   float alpha,  float *A,  INTT lda,  float *X, 
		   INTT incX,  float beta,float *Y,   INTT incY) {
	   cblas_ssymv(order,Uplo,N,alpha,A,lda,X,incX,beta,Y,incY);
	   //arma_fortran(ssymv)(cblas_uplo(Uplo),&N,&alpha,A,&lda,X,&incX,&beta,Y,&incY);
	}

	inline void cblas_symv( CBLAS_LAYOUT order,
		   CBLAS_UPLO Uplo,  INTT N, 
		   double alpha,  double *A,  INTT lda,  double *X, 
		   INTT incX,  double beta,double *Y,   INTT incY) {
	   cblas_dsymv(order,Uplo,N,alpha,A,lda,X,incX,beta,Y,incY);
	   //arma_fortran(dsymv)(cblas_uplo(Uplo),&N,&alpha,A,&lda,X,&incX,&beta,Y,&incY);
	}

	inline void cblas_scal( INTT n,  double a, double* X,
		   INTT incX) {
	   cblas_dscal(n,a,X,incX);
	   //arma_fortran(dscal)(&n,&a,X,&incX);
	};


	inline void cblas_axpy( INTT n,  double a,  double* X, 
		   INTT incX, double* Y,  INTT incY) {
	   cblas_daxpy(n,a,X,incX,Y,incY);
	   //arma_fortran(daxpy)(&n,&a,X,&incX,Y,&incY);
	};

	inline void cblas_gemv( CBLAS_LAYOUT order,
		   CBLAS_TRANSPOSE TransA,  INTT M,  INTT N,
		   double alpha,  double *A,  INTT lda,
		   double *X,  INTT incX,  double beta,
		  double *Y,  INTT incY) {
	   cblas_dgemv(order,TransA,M,N,alpha,A,lda,X,incX,beta,Y,incY);
	   //arma_fortran(dgemv)(cblas_transpose(TransA),&M,&N,&alpha,A,&lda,X,&incX,&beta,Y,&incY);
	};


	inline void cblas_ger( CBLAS_LAYOUT order, 
		   INTT M,  INTT N,  double alpha,  double *X,  INTT incX,
		   double* Y,  INTT incY, double *A,  INTT lda) {
	   cblas_dger(order,M,N,alpha,X,incX,Y,incY,A,lda);
	   //arma_fortran(dger)(&M,&N,&alpha,X,&incX,Y,&incY,A,&lda);
	};
	
#endif
