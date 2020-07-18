#include "ACTIONet.h"

#include "cholmod.h"
#include <Rcpp.h>


typedef struct cholmod_sparse_struct  *CHM_SP ;
cholmod_common chol_c;
/* Need our own CHOLMOD error handler */
void attribute_hidden
irlba_R_cholmod_error (int status, const char *file, int line,
					   const char *message)
{
  if (status < 0)
	Rcpp::stop ("Cholmod error '%s' at file:%s, line %d", message, file, line);
  else
	Rcpp::warning ("Cholmod warning '%s' at file:%s, line %d", message, file, line);
}

namespace ACTIONet {
	// Adopted from the irlba R package
	bool check_sorted_chm(CHM_SP A) {
		int *Ai = (int*)(A->i), *Ap = (int*)(A->p);
		int j, p;

		for (j = 0; j < A->ncol; j++) {
		int p1 = Ap[j], p2 = Ap[j + 1] - 1;
		for (p = p1; p < p2; p++)
			if (Ai[p] >= Ai[p + 1])
			return false;
		}
		return true;
	}
	
	CHM_SP as_cholmod_sparse(CHM_SP ans, sp_mat& A) {
		ans->itype = CHOLMOD_INT;	/* characteristics of the system */
		ans->dtype = CHOLMOD_DOUBLE;
		ans->packed = true;
		
		A.sync();		
		
		/*	
		std::vector<double> x(A.values, A.values + A.n_nonzero ) ;
		std::vector<int> i(A.row_indices, A.row_indices + A.n_nonzero);
		std::vector<int> p(A.col_ptrs, A.col_ptrs + A.n_cols+1 ) ;		

		ans->x = x.data();
		ans->i = i.data();
		ans->p = p.data();
		*/

		ans->i = new int[A.n_nonzero];
		ans->x = new double[A.n_nonzero];
		double *x_ptr = (double *)ans->x;
		int *i_ptr = (int *)ans->i;
		for(int k = 0; k < A.n_nonzero; k++) {
			x_ptr[k] = A.values[k];
			i_ptr[k] = A.row_indices[k];
		}
		
		ans->p = new int[(A.n_cols+1)];
		int *ptr = (int *)ans->p;
		for(int k = 0; k < A.n_cols+1; k++) {
			ptr[k] = A.col_ptrs[k];
		}
		
		ans->nrow = A.n_rows;
		ans->ncol = A.n_cols;		
		ans->nzmax = A.n_nonzero;
	
		ans->xtype = CHOLMOD_REAL;
		ans->stype = 0;
		ans->dtype = 0;
		
		if(!check_sorted_chm(ans)) {
			cholmod_common c; 			
			if (!cholmod_sort(ans, &c))
				printf("Couldn't sort sparse matrix indices!\n");
		}
		ans->sorted = 1;

		return ans;
	}


	void dsdmult (char transpose, int m, int n, void * a, double *b, double *c) {
		int t = transpose == 't' ? 1 : 0;
		CHM_SP cha = (CHM_SP) a;

		cholmod_dense chb;
		chb.nrow = transpose == 't' ? m : n;
		chb.d = chb.nrow;
		chb.ncol = 1;
		chb.nzmax = chb.nrow;
		chb.xtype = cha->xtype;
		chb.dtype = 0;
		chb.x = (void *) b;
		chb.z = (void *) NULL;

		cholmod_dense chc;
		chc.nrow = transpose == 't' ? n : m;
		chc.d = chc.nrow;
		chc.ncol = 1;
		chc.nzmax = chc.nrow;
		chc.xtype = cha->xtype;
		chc.dtype = 0;
		chc.x = (void *) c;
		chc.z = (void *) NULL;

		double one[] = { 1, 0 }, zero[] = { 0, 0};
		cholmod_sdmult(cha, t, one, zero, &chb, &chc, &chol_c);		
	}
	
	
	
	field<mat> IRLB_SVD(sp_mat &A, int dim, int iters = 1000, int seed = 0) {
		
		int m = A.n_rows;
		int n = A.n_cols;
		printf("\t\t* IRLB (sparse) -- A: %d x %d\n", m, n); fflush(stdout);
		
		
		
		cholmod_start (&chol_c);
		chol_c.final_ll = 1;          /* LL' form of simplicial factorization */
		chol_c.error_handler = irlba_R_cholmod_error;
		
		CHM_SP AS = new cholmod_sparse_struct;		
		as_cholmod_sparse(AS, A);
	
		
		double eps = 3e-13;
		//double eps = 2.22e-16;
		double tol = 1e-05, svtol = 1e-5;

		srand(seed);
		std::default_random_engine gen (seed);	
		std::normal_distribution<double> normDist(0.0, 1.0);		
		
		int work = dim + 7;
		int lwork = 7 * work * (1 + work);
		
		double *s = new double[dim];
		double *U = new double[m*work];
		double *V = new double[n*work];


		double *V1 = new double[n * work];  
		double *U1 = new double[m * work];
		double *W = new double[m * work];  
		double *F = new double[n];
		double *B = new double[work * work];
		double *BU = new double[work * work];
		double *BV = new double[work * work];
		double *BS = new double[work];
		double *BW = new double[lwork];
		double *res = new double[work];
		double *T = new double[lwork];
		double *svratio = new double[work];
		
		
		
		
		double d, S, R, R_F, SS;
		double *x;
		int jj, kk;
		int converged;
		int info, j, k = 0;
		int iter = 0;
		double Smax = 0;


		memset (B, 0, work * work * sizeof (double));	
		memset(svratio, 0, work * sizeof(double));

		double alpha = 1, beta = 0;
		int inc = 1;

		
		vec v, y;

	
		// Initialize first column of V
		randN_BM(V, n);
		/*
        for ( int i = 0; i < n; i ++ ) {
            V[i]   = normDist(gen);;
        }   		
		*/
		
		
		/* Main iteration */
		while (iter < iters) {	
				
			j = 0;
			
			/*  Normalize starting vector */
			if (iter == 0) {
				d = cblas_dnrm2(n, V, inc);
				d = 1 / d;
				cblas_dscal(n, d, V, inc);			
			}
			else
				j = k;
			
			// Compute Ax
			x = V + j * n;
			/*v = vec(x, n, true);
			y = A * v;			
			memcpy(W + j * m, y.memptr(), y.n_elem*sizeof(double));		
			*/
			dsdmult ('n', m, n, AS, x, W + j * m);

									  
			if (iter > 0)
				orthog (W, W + j * m, T, m, j, 1);
			
			S = cblas_dnrm2(m, W + j * m, inc);
			SS = 1.0 / S;		
			cblas_dscal(m, SS, W + j * m, inc);

			/* The Lanczos process */
			while (j < work) {			
				/*
				v = vec(W + j * m, m, true);				
				y = At*v;				
				memcpy(F, y.memptr(), y.n_elem*sizeof(double));
				*/
				dsdmult ('t', m, n, AS, W + j * m, F);
				//v = vec(F, A.n_cols, true);				
				//v(span(0, 5)).print("F");

				
				SS = -S;
				cblas_daxpy(n, SS, V + j * n, inc, F, inc);
				orthog (V, F, T, n, j + 1, 1);

				if (j + 1 < work) {
					R_F = cblas_dnrm2(n, F, inc);
					R = 1.0 / R_F;
				  
					if (R_F < eps) {        // near invariant subspace
						//randN_BM(F, n);
						for ( int i = 0; i < n; i ++ ) {
							F[i]   = normDist(gen);;
						}   						

						orthog (V, F, T, n, j + 1, 1);
						R_F = cblas_dnrm2(n, F, inc);
						R = 1.0 / R_F;
						R_F = 0;
					}
				  
					memmove (V + (j + 1) * n, F, n * sizeof (double));
					cblas_dscal(n, R, V + (j + 1) * n, inc);
					B[j * work + j] = S;
					B[(j + 1) * work + j] = R_F;
					
					x = V + (j + 1) * n;
					/*
					v = vec(x, n, true);
					y = A*v;
					memcpy(W + (j + 1) * m, y.memptr(), y.n_elem*sizeof(double));
					*/
				    dsdmult ('n', m, n, AS, x, W + (j + 1) * m);

					/* One step of classical Gram-Schmidt */
					R = -R_F;
					cblas_daxpy(m, R, W + j * m, inc, W + (j + 1) * m, inc);				
					
					/* full re-orthogonalization of W_{j+1} */				
					orthog (W, W + (j + 1) * m, T, m, j + 1, 1);
					S = cblas_dnrm2 (m, W + (j + 1) * m, inc);
					SS = 1.0 / S;
					
					if (S < eps) {
						jj = (j + 1) * m;
						//randN_BM(W+jj, m);
						for ( int i = 0; i < m; i ++ ) {
							W[jj+i]   = normDist(gen);;
						}   						
						
						orthog (W, W + (j + 1) * m, T, m, j + 1, 1);
						S = cblas_dnrm2(m, W + (j + 1) * m, inc);
						SS = 1.0 / S;
						cblas_dscal(m, SS, W + (j + 1) * m, inc);
						S = 0;
					}                
					else
						cblas_dscal(m, SS, W + (j + 1) * m, inc);
				}
				else {
					B[j * work + j] = S;
				}

				
				j++;
			}			
			
			mat tmp(B, work, work, false);			
			mat Umat(BU, work, work, false);
			vec svec(BS, work, false);
			mat Vmat(BV, work, work, false);

			arma::svd(Umat, svec, Vmat, tmp, "dc");
			Vmat = trans(Vmat);
			
			/*
			Umat(span(0, 5), span(0, 5)).print("U1");
			Vmat(span(0, 5), span(0, 5)).print("V1");
			* 
			Umat(span(0, 5), span(0, 5)).print("tmp (after)");
			
			for(int i = 0; i < 20; i++) {
				printf("%d- %f\n", i, BU[i]);
			}	
									
			memmove (BU, B, work * work * sizeof (double));   // Make a working copy of B
			int *BI = (int *) T;
			F77_NAME (dgesdd) ("O", &work, &work, BU, &work, BS, BU, &work, BV, &work, BW, &lwork, BI, &info);			

			mat Umat2(BU, work, work, false);
			mat Vmat2(BV, work, work, false);
			Umat2(span(0, 5), span(0, 5)).print("U2");
			Vmat2(span(0, 5), span(0, 5)).print("V2");
			*/

			
			/*
			for(int i = 0; i < 20; i++) {
				printf("%d, %d- %f\n", iter, i, BU[i]);
			}	
			*/
					
			R_F = cblas_dnrm2(n, F, inc);
			R = 1.0 / R_F;
			cblas_dscal(n, R, F, inc);
			
			/* Force termination after encountering linear dependence */
			if (R_F < eps)
				R_F = 0;

			Smax = 0;
			for (jj = 0; jj < j; ++jj) {
				if (BS[jj] > Smax)
				Smax = BS[jj];
				svratio[jj] = fabs (svratio[jj] - BS[jj]) / BS[jj];
			}
			
			
			for (kk = 0; kk < j; ++kk)
				res[kk] = R_F * BU[kk * work + (j - 1)];
				
			/* Update k to be the number of converged singular values. */
			convtests (j, dim, tol, svtol, Smax, svratio, res, &k, &converged, S);
			if (converged == 1) {
				iter++;
				break;
			}
			
			for (jj = 0; jj < j; ++jj)
				svratio[jj] = BS[jj];

			cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, k, j, alpha, V, n, BV, work, beta, V1, n);
			
			memmove (V, V1, n * k * sizeof (double));
			memmove (V + n * k, F, n * sizeof (double));

			memset (B, 0, work * work * sizeof (double));
			for (jj = 0; jj < k; ++jj) {
				B[jj * work + jj] = BS[jj];
				B[k * work + jj] = res[jj];
			}

			/*   Update the left approximate singular vectors */
			cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, k, j, alpha, W, m, BU, work, beta, U1, m);
			
			memmove (W, U1, m * k * sizeof (double));
			iter++;
		}

		/* Results */
		memmove (s, BS, dim * sizeof (double));        /* Singular values */
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, dim, work, alpha, W, m, BU, work, beta, U, m);
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, dim, work, alpha, V, n, BV, work, beta, V1, n);
		memmove (V, V1, n * dim * sizeof (double));

		field<mat> out(3);			
		out(0) = mat(U, m, dim);
		out(1) = vec(s, dim);
		out(2) = mat(V, n, dim);


		delete [] s;
		delete [] U;
		delete [] V;	
		delete [] V1;
		delete [] U1;
		delete [] W;
		delete [] F;
		delete [] B;
		delete [] BU;
		delete [] BV;
		delete [] BS;
		delete [] BW;
		delete [] res;
		delete [] T;
		delete [] svratio;	
		
		delete [] AS->x;
		delete [] AS->i;
		delete [] AS->p;
		delete AS; 

		cholmod_finish (&chol_c);
		
		if (converged != 1) {
			fprintf(stderr, "IRLB_SVD did NOT converge! Try in creasing the number of iterations\n");
		}		
		
		return(out);
	}
	

	field<mat> IRLB_SVD(mat &A, int dim, int iters = 1000, int seed = 0) {
		printf("\t\t* IRLB (dense) -- A: %d x %d\n", A.n_rows, A.n_cols); fflush(stdout);
		
		double eps = 3e-13;
		//double eps = 2.22e-16;
		double tol = 1e-05, svtol = 1e-5;

		srand(seed);
		
		int m = A.n_rows;
		int n = A.n_cols;
		int work = dim + 7;
		int lwork = 7 * work * (1 + work);
		
		double *s = new double[dim];
		double *U = new double[m*work];
		double *V = new double[n*work];


		double *V1 = new double[n * work];  
		double *U1 = new double[m * work];
		double *W = new double[m * work];  
		double *F = new double[n];
		double *B = new double[work * work];
		double *BU = new double[work * work];
		double *BV = new double[work * work];
		double *BS = new double[work];
		double *BW = new double[lwork];
		double *res = new double[work];
		double *T = new double[lwork];
		double *svratio = new double[work];
		
		
		
		
		double d, S, R, R_F, SS;
		double *x;
		int jj, kk;
		int converged;
		int info, j, k = 0;
		int iter = 0;
		double Smax = 0;


		memset (B, 0, work * work * sizeof (double));	
		memset(svratio, 0, work * sizeof(double));

		double alpha = 1, beta = 0;
		int inc = 1;

		
		vec v, y;

	
		// Initialize first column of V
		randN_BM(V, n);
		
		/* Main iteration */
		while (iter < iters) {			
			j = 0;
			
			/*  Normalize starting vector */
			if (iter == 0) {
				d = cblas_dnrm2(n, V, inc);
				d = 1 / d;
				cblas_dscal(n, d, V, inc);
			}
			else
				j = k;
			
			// Compute Ax
			x = V + j * n;
			v = vec(x, n, true);
			y = A * v;			
			memcpy(W + j * m, y.memptr(), y.n_elem*sizeof(double));		
			
									  
			if (iter > 0)
				orthog (W, W + j * m, T, m, j, 1);
			
			S = cblas_dnrm2(m, W + j * m, inc);
			SS = 1.0 / S;		
			cblas_dscal(m, SS, W + j * m, inc);

			/* The Lanczos process */
			while (j < work) {			
				v = vec(W + j * m, m, true);				
				y = trans(trans(v)*A);				
				memcpy(F, y.memptr(), y.n_elem*sizeof(double));
				
				
				SS = -S;
				cblas_daxpy(n, SS, V + j * n, inc, F, inc);
				orthog (V, F, T, n, j + 1, 1);

				if (j + 1 < work) {
					R_F = cblas_dnrm2(n, F, inc);
					R = 1.0 / R_F;
				  
					if (R_F < eps) {        // near invariant subspace
						randN_BM(F, n);

						orthog (V, F, T, n, j + 1, 1);
						R_F = cblas_dnrm2(n, F, inc);
						R = 1.0 / R_F;
						R_F = 0;
					}
				  
					memmove (V + (j + 1) * n, F, n * sizeof (double));
					cblas_dscal(n, R, V + (j + 1) * n, inc);
					B[j * work + j] = S;
					B[(j + 1) * work + j] = R_F;
					
					x = V + (j + 1) * n;
					v = vec(x, n, true);
					y = A*v;
					memcpy(W + (j + 1) * m, y.memptr(), y.n_elem*sizeof(double));
				  
					/* One step of classical Gram-Schmidt */
					R = -R_F;
					cblas_daxpy(m, R, W + j * m, inc, W + (j + 1) * m, inc);				
					
					/* full re-orthogonalization of W_{j+1} */				
					orthog (W, W + (j + 1) * m, T, m, j + 1, 1);
					S = cblas_dnrm2 (m, W + (j + 1) * m, inc);
					SS = 1.0 / S;
					
					if (S < eps) {
						jj = (j + 1) * m;
						randN_BM(W+jj, m);

						orthog (W, W + (j + 1) * m, T, m, j + 1, 1);
						S = cblas_dnrm2(m, W + (j + 1) * m, inc);
						SS = 1.0 / S;
						cblas_dscal(m, SS, W + (j + 1) * m, inc);
						S = 0;
					}                
					else
						cblas_dscal(m, SS, W + (j + 1) * m, inc);
				}
				else {
					B[j * work + j] = S;
				}


				j++;
			}
			
			
			mat tmp(B, work, work, false);			
			mat Umat(BU, work, work, false);
			vec svec(BS, work, false);
			mat Vmat(BV, work, work, false);

			arma::svd(Umat, svec, Vmat, tmp, "dc");
			Vmat = trans(Vmat);
			
			/*
			Umat(span(0, 5), span(0, 5)).print("U1");
			Vmat(span(0, 5), span(0, 5)).print("V1");
			* 
			Umat(span(0, 5), span(0, 5)).print("tmp (after)");
			
			for(int i = 0; i < 20; i++) {
				printf("%d- %f\n", i, BU[i]);
			}	
									
			memmove (BU, B, work * work * sizeof (double));   // Make a working copy of B
			int *BI = (int *) T;
			F77_NAME (dgesdd) ("O", &work, &work, BU, &work, BS, BU, &work, BV, &work, BW, &lwork, BI, &info);			

			mat Umat2(BU, work, work, false);
			mat Vmat2(BV, work, work, false);
			Umat2(span(0, 5), span(0, 5)).print("U2");
			Vmat2(span(0, 5), span(0, 5)).print("V2");
			*/


			R_F = cblas_dnrm2(n, F, inc);
			R = 1.0 / R_F;
			cblas_dscal(n, R, F, inc);
			
			/* Force termination after encountering linear dependence */
			if (R_F < eps)
				R_F = 0;

			Smax = 0;
			for (jj = 0; jj < j; ++jj) {
				if (BS[jj] > Smax)
				Smax = BS[jj];
				svratio[jj] = fabs (svratio[jj] - BS[jj]) / BS[jj];
			}
			
			
			for (kk = 0; kk < j; ++kk)
				res[kk] = R_F * BU[kk * work + (j - 1)];
				
			/* Update k to be the number of converged singular values. */
			convtests (j, dim, tol, svtol, Smax, svratio, res, &k, &converged, S);
			if (converged == 1) {
				iter++;
				break;
			}
			
			for (jj = 0; jj < j; ++jj)
				svratio[jj] = BS[jj];

			cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, k, j, alpha, V, n, BV, work, beta, V1, n);
			
			memmove (V, V1, n * k * sizeof (double));
			memmove (V + n * k, F, n * sizeof (double));

			memset (B, 0, work * work * sizeof (double));
			for (jj = 0; jj < k; ++jj) {
				B[jj * work + jj] = BS[jj];
				B[k * work + jj] = res[jj];
			}

			/*   Update the left approximate singular vectors */
			cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, k, j, alpha, W, m, BU, work, beta, U1, m);
			
			memmove (W, U1, m * k * sizeof (double));
			iter++;
		}

		/* Results */
		memmove (s, BS, dim * sizeof (double));        /* Singular values */
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, dim, work, alpha, W, m, BU, work, beta, U, m);
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, dim, work, alpha, V, n, BV, work, beta, V1, n);
		memmove (V, V1, n * dim * sizeof (double));

		field<mat> out(3);			
		out(0) = mat(U, m, dim);
		out(1) = vec(s, dim);
		out(2) = mat(V, n, dim);


		delete [] s;
		delete [] U;
		delete [] V;	
		delete [] V1;
		delete [] U1;
		delete [] W;
		delete [] F;
		delete [] B;
		delete [] BU;
		delete [] BV;
		delete [] BS;
		delete [] BW;
		delete [] res;
		delete [] T;
		delete [] svratio;	
		
		if (converged != 1) {
			fprintf(stderr, "IRLB_SVD did NOT converge! Try in creasing the number of iterations\n");
		}		
		
		return(out);
	}
	


	//****************************************************************************************************************************************************************************
	// From: Xu Feng, Yuyang Xie, and Yaohang Li, "Fast Randomzied SVD for Sparse Data," in Proc. the 10th Asian Conference on Machine Learning (ACML), Beijing, China, Nov. 2018.
	//****************************************************************************************************************************************************************************
	field<mat> FengSVD(sp_mat &A, int dim, int iters, int seed = 0) {	
		int s = 5;

		int m = A.n_rows;
		int n = A.n_cols;
		
		printf("\t\t* Feng (sparse) -- A: %d x %d\n", A.n_rows, A.n_cols); fflush(stdout);
				
		vec S;
		mat Q, L, U, V;
		field<mat> SVD_out;
		
		if(m < n) {
			printf("\t\t\tInitializing SVD (mode 1) ... ");
			//arma_rng::set_seed(seed);			
			//Q = randn( n, dim+s );
			//Q = sampleUnif(n, dim+s, 0.0, 1.0, seed);
			Q = randNorm(n, dim+s, seed);
			Q = A*Q;
			if (iters == 0) {
				SVD_out = eigSVD(Q);
				Q = SVD_out(0);					
			}
			else {
				lu(L, U, Q);
				Q = L;
			}
			printf("done\n");
			
			for (int i = 1; i <= iters; i++) {
				printf("\t\t\t\tIter %d/%d ... ", i, iters);				
				
				if (i == iters) {
					SVD_out = eigSVD(A*(trans(A)*Q));
					Q = SVD_out(0);									
				}
				else {
					lu(L, U, A*(trans(A)*Q));
					Q = L;
				}
				printf("done\n");				
			}
			
			SVD_out = eigSVD(trans(A)*Q);
			V = SVD_out(0);
			S = vec(SVD_out(1));
			U = SVD_out(2);
			
			U = Q*fliplr(U.cols(s, dim+s-1));
			V = fliplr(V.cols(s, dim+s-1));
			S = flipud(S(span(s, dim+s-1)));
		}
		else {
			printf("\t\t\tInitializing SVD (mode 2) ... ");				
			// arma_rng::set_seed(seed);
			// Q = randn( m, dim+s ); 
			// Q = sampleUnif(m, dim+s, 0.0, 1.0, seed);
			Q = randNorm(m, dim+s, seed);
			Q = trans(A)*Q;
			if (iters == 0) {
				SVD_out = eigSVD(Q);
				Q = SVD_out(0);					
			}
			else {
				lu(L, U, Q);
				Q = L;
			}
			printf("done\n");
			
			for (int i = 1; i <= iters; i++) {
				printf("\t\t\t\tIter %d/%d ... ", i, iters);				
				
				if (i == iters) {
					SVD_out = eigSVD(trans(A)*(A*Q));
					Q = SVD_out(0);									
				}
				else {
					lu(L, U, trans(A)*(A*Q));
					Q = L;
				}
				printf("done\n");				
				
			}
			
			SVD_out = eigSVD(A*Q);
			U = SVD_out(0);
			S = vec(SVD_out(1));
			V = SVD_out(2);
						
			
			U = fliplr(U.cols(s, dim+s-1));
			V = Q*fliplr(V.cols(s, dim+s-1));
			S = flipud(S(span(s, dim+s-1)));
		}		
	
	
		field<mat> out(3);		
		out(0) = U;
		out(1) = S;
		out(2) = V;
		
		printf("\t\t\tdone\n");	fflush(stdout);			
		
		return(out);
	}


	// From: Xu Feng, Yuyang Xie, and Yaohang Li, "Fast Randomzisped SVD for Sparse Data," in Proc. the 10th Asian Conference on Machine Learning (ACML), Beijing, China, Nov. 2018.
	field<mat> FengSVD(mat &A, int dim, int iters, int seed = 0) {	
		int s = 5;

		int m = A.n_rows;
		int n = A.n_cols;
		
		printf("\t\t* Feng (dense) -- A: %d x %d\n", A.n_rows, A.n_cols); fflush(stdout);
				
		vec S;
		mat Q, L, U, V;
		field<mat> SVD_out;
		
		if(m < n) {
			printf("\t\t\tInitializing SVD (mode 1) ... ");
			//arma_rng::set_seed(seed);			
			//Q = randn( n, dim+s );
			//Q = sampleUnif(n, dim+s, 0.0, 1.0, seed);
			Q = randNorm(n, dim+s, seed);

			Q = A*Q;
			if (iters == 0) {
				SVD_out = eigSVD(Q);
				Q = SVD_out(0);					
			}
			else {
				lu(L, U, Q);
				Q = L;
			}
			printf("done\n");
			
			for (int i = 1; i <= iters; i++) {
				printf("\t\t\t\tIter %d/%d ... ", i, iters);				
				
				if (i == iters) {
					SVD_out = eigSVD(A*(trans(A)*Q));
					Q = SVD_out(0);									
				}
				else {
					lu(L, U, A*(trans(A)*Q));
					Q = L;
				}
				printf("done\n");				
			}
			
			SVD_out = eigSVD(trans(A)*Q);
			V = SVD_out(0);
			S = vec(SVD_out(1));
			U = SVD_out(2);
			
			U = Q*fliplr(U.cols(s, dim+s-1));
			V = fliplr(V.cols(s, dim+s-1));
			S = flipud(S(span(s, dim+s-1)));
		}
		else {
			printf("\t\t\tInitializing SVD (mode 2) ... ");				
			// arma_rng::set_seed(seed);
			// Q = randn( m, dim+s ); 
			//Q = sampleUnif(m, dim+s, 0.0, 1.0, seed);
			Q = randNorm(m, dim+s, seed);
			Q = trans(A)*Q;
			if (iters == 0) {
				SVD_out = eigSVD(Q);
				Q = SVD_out(0);					
			}
			else {
				lu(L, U, Q);
				Q = L;
			}
			printf("done\n");
			
			for (int i = 1; i <= iters; i++) {
				printf("\t\t\t\tIter %d/%d ... ", i, iters);				
				
				if (i == iters) {
					SVD_out = eigSVD(trans(A)*(A*Q));
					Q = SVD_out(0);									
				}
				else {
					lu(L, U, trans(A)*(A*Q));
					Q = L;
				}
				printf("done\n");				
				
			}
			
			SVD_out = eigSVD(A*Q);
			U = SVD_out(0);
			S = vec(SVD_out(1));
			V = SVD_out(2);
						
			
			U = fliplr(U.cols(s, dim+s-1));
			V = Q*fliplr(V.cols(s, dim+s-1));
			S = flipud(S(span(s, dim+s-1)));
		}		
	
	
		field<mat> out(3);		
		out(0) = U;
		out(1) = S;
		out(2) = V;
		
		printf("\t\t\tdone\n");	fflush(stdout);			
		
		return(out);
	}

	//**************************************************************************************************************************************************************************************************
	// From: N Halko, P. G Martinsson, and J. A Tropp. Finding structure with randomness: Probabilistic algorithms for constructing approximate matrix decompositions. Siam Review, 53(2):217-288, 2011.
	//**************************************************************************************************************************************************************************************************
	field<mat> HalkoSVD(sp_mat &A, int dim, int iters, int seed = 0) {	
		field<mat> results(3);

		int m = A.n_rows;
		int n = A.n_cols;
		int l = dim + 2;		
				
		vec s;
		mat R, Q;
		mat U, V, X;
		
		printf("\t\t* Halko (sparse) -- A: %d x %d\n", A.n_rows, A.n_cols); fflush(stdout);
		
		if (m < n) {
			//R = stats::runif<arma::mat>(l, m, -1.0, 1.0, seed);
			// R = sampleUnif(l, m, -1.0, 1.0, 0);
			R = randNorm(l, m, seed);
			
			sp_mat At = A.t();
			Q = At*R.t(); 
		}
		else {
			//R = stats::runif<arma::mat>(n, l, -1.0, 1.0, seed);
			//R = sampleUnif(n, l, -1.0, 1.0, 0);
			R = randNorm(n, l, seed);
			
			Q = A*R;
		}				


		// Form a matrix Q whose columns constitute a well-conditioned basis for the columns of the earlier Q.			
		gram_schmidt(Q);
		//Q = orth(Q);
		
		if (m < n) {
			// Conduct normalized power iterations.
			for(int it = 1; it <= iters; it++) {
				printf("\t\tIteration %d\n", it);
				
				
				Q = A*Q; 
				gram_schmidt(Q);
				//Q = orth(Q);

				Q = A.t()*Q; 
				gram_schmidt(Q);								
				//Q = orth(Q);
			}

			X = mat(A*Q);
			printf("\t\tReduced SVD ... ");
			svd_econ( U, s, V,  X);
			printf("done\n");
			V = Q*V;		
		}				
		else {
			// Conduct normalized power iterations.
			for(int it = 1; it <= iters; it++) {
				printf("\t\tIteration %d\n", it);
				
				Q = A.t()*Q;
				gram_schmidt(Q);
				//Q = orth(Q);
					
				Q = A*Q; // Apply A to a random matrix, obtaining Q.
				gram_schmidt(Q);			
				//Q = orth(Q);
			}
			
			// SVD Q' applied to the centered A to obtain approximations to the singular values and right singular vectors of the A;
			
			X = mat(Q.t()*A);
			svd_econ( U, s, V,  X);
			U = Q*U;		
		}

		U.shed_cols(dim, dim+1);
		s = s(span(0, dim-1));
		V.shed_cols(dim, dim+1);
		
		results(0) = U;
		results(1) = s;
		results(2) = V;
		
		return results;	
	}



	field<mat> HalkoSVD(mat &A, int dim, int iters, int seed = 0) {	
		field<mat> results(3);

		int m = A.n_rows;
		int n = A.n_cols;
		int l = dim + 2;		
				
		vec s;
		mat R, Q;
		mat U, V, X;
		
		printf("\t\t* Halko (dense) -- A: %d x %d\n", A.n_rows, A.n_cols); fflush(stdout);
		
		if (m < n) {
			//R = stats::runif<arma::mat>(l, m, -1.0, 1.0, seed);
			//R = sampleUnif(l, m, -1.0, 1.0, 0);
			R = randNorm(l, m, seed);
			
			mat At = A.t();
			Q = At*R.t(); 
		}
		else {
			//R = stats::runif<arma::mat>(n, l, -1.0, 1.0, seed);
			//R = sampleUnif(n, l, -1.0, 1.0, 0);
			R = randNorm(n, l, seed);
			
			Q = A*R;
		}				


		// Form a matrix Q whose columns constitute a well-conditioned basis for the columns of the earlier Q.			
		gram_schmidt(Q);
		//Q = orth(Q);
		
		if (m < n) {
			// Conduct normalized power iterations.
			for(int it = 1; it <= iters; it++) {
				printf("\t\tIteration %d\n", it);
				
				
				Q = A*Q; 
				gram_schmidt(Q);
				//Q = orth(Q);

				Q = A.t()*Q; 
				gram_schmidt(Q);								
				//Q = orth(Q);
			}

			X = mat(A*Q);
			printf("\t\tReduced SVD ... ");
			svd_econ( U, s, V,  X);
			printf("done\n");
			V = Q*V;		
		}				
		else {
			// Conduct normalized power iterations.
			for(int it = 1; it <= iters; it++) {
				printf("\t\tIteration %d\n", it);
				
				Q = A.t()*Q;
				gram_schmidt(Q);
				//Q = orth(Q);
					
				Q = A*Q; // Apply A to a random matrix, obtaining Q.
				gram_schmidt(Q);			
				//Q = orth(Q);
			}
			
			// SVD Q' applied to the centered A to obtain approximations to the singular values and right singular vectors of the A;
			
			X = mat(Q.t()*A);
			svd_econ( U, s, V,  X);
			U = Q*U;		
		}

		U.shed_cols(dim, dim+1);
		s = s(span(0, dim-1));
		V.shed_cols(dim, dim+1);
		
		results(0) = U;
		results(1) = s;
		results(2) = V;
		
		return results;	
	}		
}
