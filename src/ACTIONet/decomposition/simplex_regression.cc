#include <ACTIONet.h>
//#include <mini_cblas.h>
#include <cassert>

// Re-implemented from: Fast and Robust Archetypal Analysis for Representation Learning
namespace ACTIONet {
		
	/* **************************
	 * Active-Set Method with direct inversion, with update(matrix inversion lemma)
	 * **************************/

	vec activeSet_arma(arma::mat& M, arma::vec& b, double lambda2 = double(1e-5), double epsilon = double(1e-5)) {
	   int m = M.n_rows;
	   int p = M.n_cols;
	   int L = min(m,p)+1;
	   
	   mat Mt = trans(M);

	   arma::vec c = -Mt*b;
	   
	   
	   double lam2sq =lambda2*lambda2;

	   arma::vec x(p);
	   double* pr_M = M.memptr();
	   // constraint matrix
	   arma::vec A = ones(L);	   
	   double* pr_A = A.memptr();
	   
	   // Non-Active Constraints Set
	   arma::ivec NASet(L); NASet.ones(); NASet = -NASet;
	   arma::ivec NAMask(p); NAMask.zeros();
	   
	   int na;
	   arma::vec xRed(L);
	   arma::vec cRed(L);
	   arma::mat MRed(m,L);
	   double* pr_MRed = MRed.memptr();
	   arma::mat GRed(L,L);
	   double* pr_GRed = GRed.memptr();
	   arma::mat GRedinv(L, L);
	   double* pr_GRedinv = GRedinv.memptr();
	   
	   // cold-start
	  x.zeros();
	  x[0] = double(1.0);
	  // Non-Active Constraints Set
	  NASet[0] = 0;
	  NAMask[0] = 1;
	  na = 1;
	  xRed[0] = x[0];
	  cRed[0] = c[0];
	  cblas_dcopy(m, pr_M, 1, pr_MRed, 1);

	  // BLAS GRed = MRedT * MRed + lam2sq (na = 1 for now)
	  double coeff = cblas_ddot(m,pr_MRed,1,pr_MRed,1) + lam2sq;
	  GRed(0,0) = coeff;
	  GRedinv(0,0)= double(1.0) / GRed(0,0);


	   arma::vec Gplus = Mt*(M*x) + lam2sq*x + c;
	   double* pr_Gplus = Gplus.memptr();
	   
	   

	   arma::vec gRed(L);
	   double* pr_gRed = gRed.memptr();
	   arma::vec MRedxRed(m);

	   arma::vec GinvA(L);
	   double* pr_GinvA = GinvA.memptr();
	   arma::vec Ginvg(L);
	   double* pr_Ginvg = Ginvg.memptr();
	   arma::vec PRed(L);
	   double* pr_PRed = PRed.memptr();
	   arma::vec MRedPRed(m);
	   double* pr_MRedPRed = MRedPRed.memptr();
	   arma::vec UB(L);
	   double* pr_UB = UB.memptr();
	   arma::vec UAiB(L);
	   double* pr_UAiB = UAiB.memptr();
	   // main loop active set
	   int iter = 0;
	   while(iter <= 100*p) {
		  ++iter;
		  // update of na, NASet, NAMask, xRed, cRed, Gplus, GRedinv, already done
		  // now update MRed, gRed, GinvA, Ginvg  (no need to update GRed)
		  // MRed
		  for(int i = 0; i < na; ++i) {
			 // BLAS copy first columns of M to MRed
			 cblas_dcopy(m, pr_M + m*NASet[i], 1, pr_MRed + m*i , 1);
		  }
		  // gRed
		  for(int i = 0; i<na; ++i) {
			 gRed[i] = Gplus[NASet[i]];
		  }

		  // GinvA
		  // BLAS GinvA = GRedinv * A (ARed == A)
		  cblas_dsymv(CblasColMajor,CblasUpper,na,double(1.0),pr_GRedinv,L,pr_A,1,double(),pr_GinvA,1);
		  // Ginvg
		  // BLAS Ginvg = GRedinv * gRed
		  cblas_dsymv(CblasColMajor,CblasUpper,na,double(1.0),pr_GRedinv,L,pr_gRed,1,double(),pr_Ginvg,1);
		  double sGinvg = double();
		  double sGinvA = double();
		  for(int i = 0; i< na; ++i) {
			 sGinvg += Ginvg[i];
			 sGinvA += GinvA[i];
		  }
		  double lambdaS = sGinvg / sGinvA;
		  // BLAS PRed = GinvA * lambdaS - Ginvg
		  cblas_dcopy(na, pr_GinvA, 1, pr_PRed, 1);
		  cblas_dscal(na, lambdaS, pr_PRed, 1);
		  cblas_daxpy(na, double(-1.0), pr_Ginvg, 1, pr_PRed, 1);

		  double maxPRed = abs(PRed[0]);
		  for(int i = 0; i< na; ++i) {
			 if(abs(PRed[i])>maxPRed)
				maxPRed = abs(PRed[i]);
		  }
		  if(maxPRed < 1e-10) {
			 // P = 0, no advance possible
			 bool isOpt = true;
			 double lamMin = -epsilon;
			 int indexMin = -1;
			 for(int i = 0; i <p; ++i) {
				if(!NAMask[i] && Gplus[i]-lambdaS < lamMin) {
				   isOpt = false;
				   lamMin = Gplus[i]-lambdaS;
				   indexMin = i;
				}
			 }

			 if(isOpt) {
				// Got the optimal, STOP!
				return(x);
			 } else {
				// Add one constraint
				NAMask[indexMin] = 1;
				NASet[na] = indexMin;
				xRed[na] = x[indexMin];
				cRed[na] = c[indexMin];
				// Gplus inchange

				// update GRedinv
				// BLAS UB = MRed.double * M[:, indexMin]
				cblas_dgemv(CblasColMajor,CblasTrans,m,na,double(1.0),pr_MRed,m,pr_M+indexMin*m,1,double(),pr_UB,1);
				// BLAS UC = M[:,indexMin].double* M[:, indexMin]
				double UC = cblas_ddot(m,pr_M+indexMin*m,1, pr_M+indexMin*m,1) + lam2sq;
				// BLAS UAiB = GRedinv * UB
				cblas_dsymv(CblasColMajor,CblasUpper,na,double(1.0),pr_GRedinv,L,pr_UB,1,double(),pr_UAiB,1);
				double USi = double(1.0)/(UC - cblas_ddot(na,pr_UB,1,pr_UAiB,1));
				// GRedinv (restricted) += USi * UAiB*UAiB
				//replace cblas_syr(CblasColMajor,CblasUpper,na,USi, pr_UAiB, 1, pr_GRedinv, L);  
				cblas_dger(CblasColMajor,na,na,USi,pr_UAiB,1,pr_UAiB,1,pr_GRedinv,L);
				// copy -UAiB*USi, -UAiB.double*USi, USi to GRedinv
				cblas_dcopy(na, pr_UAiB, 1, pr_GRedinv+na*L, 1);
				cblas_dscal(na, -USi, pr_GRedinv+na*L,1);
				cblas_dcopy(na, pr_UAiB, 1, pr_GRedinv+na, L);
				cblas_dscal(na, -USi, pr_GRedinv+na,L);
				GRedinv(na,na) = USi;

				na += 1;
				assert(na <= L);
			 }
		  } else {
			 // P != 0, can advance
			 int indexMin = -1;
			 double alphaMin = double(1.0);
			 for(int i = 0; i < na; ++i) {
				if(PRed[i] < 0 && -xRed[i]/PRed[i] < alphaMin) {
				   indexMin = i;
				   alphaMin = -xRed[i]/PRed[i];
				}
			 }
			 // update x and Gplus
			 cblas_dscal(na, min(double(1.0), alphaMin), pr_PRed, 1);
			 for(int i = 0; i< na; ++i) {
				x[NASet[i]] += PRed[i];
				xRed[i] = x[NASet[i]];
				// BLAS Gplus += M.double * M[:, NASet[i]] * cAdv
				//cblas_dgemv(CblasColMajor,CblasTrans,m,p,cAdv,pr_M,m, pr_M+NASet[i]*m,1,double(1.0),pr_Gplus,1);
				// BLAS Gplus
				Gplus[NASet[i]] += PRed[i]*lam2sq;
			 }
			 // Gplus += M.double * MRed * (scaled PRed)
			 cblas_dgemv(CblasColMajor,CblasNoTrans,m,na,double(1.0),pr_MRed,m,pr_PRed,1,double(),pr_MRedPRed,1);    
			 cblas_dgemv(CblasColMajor,CblasTrans,m,p,double(1.0),pr_M,m,pr_MRedPRed,1,double(1.0),pr_Gplus,1);


			 // delete one constraint or not?
			 if(indexMin != -1) {
				// give true 0
				// x[NASet[indexMin]] = double();
				// delete one constraint
				NAMask[NASet[indexMin]] = 0;
				// downdate remove this -1;
				na -= 1;
				for(int i = indexMin; i<na; ++i) {
				   NASet[i] = NASet[i+1];
				   xRed[i] = xRed[i+1];
				   cRed[i] = cRed[i+1];
				}
				NASet[na] = -1;
				xRed[na] = double();
				cRed[na] = double();
				// PRed also
				PRed[na] = double();

				// downdate GRedinv
				double UCi = double(1.0)/GRedinv(indexMin, indexMin); 
				// BLAS UB = GRedinv[ALL\indexMin,indexMin]
				cblas_dcopy(na+1, pr_GRedinv+indexMin*L, 1, pr_UB, 1);
				for(int i = indexMin; i<na; ++i)
				   UB[i] = UB[i+1];
				UB[na] = double();
				// get (GRedinv translated)
				// column first
				for(int i = indexMin; i<na; ++i)
				   cblas_dcopy(na+1, pr_GRedinv+(i+1)*L,1,pr_GRedinv+i*L,1);
				// row then
				for(int i = indexMin; i<na; ++i)
				   cblas_dcopy(na+1, pr_GRedinv+i+1, L, pr_GRedinv+i, L);

				// BLAS GRedinv = (GRedinv translated) - UB*UB.double*UCi 
				//replace cblas_syr(CblasColMajor,CblasUpper,na,-UCi, pr_UB, 1, pr_GRedinv, L);  
				cblas_dger(CblasColMajor,na,na,-UCi,pr_UB,1,pr_UB,1,pr_GRedinv,L);
			 }
		  }
	   }
	   return(x);
	}

	/// Active-Set Method with direct inversion, with update(matrix inversion lemma)
	/// Memorize M.double* M + lam2sq = G
	vec activeSetS_arma(arma::mat& M, arma::vec& b, arma::mat& G, double lambda2 = 1e-5, double epsilon = 1e-5) {
	   int m = M.n_rows;
	   int p = M.n_cols;
	   int L = min(m,p)+1;
	   double lam2sq =lambda2*lambda2;
	   double* pr_G = G.memptr();
	   
	   mat Mt = arma::trans(M);
	   
	   arma::vec c = -Mt*b;

	   arma::vec x(p);
	   double* pr_M = M.memptr();
	   // constraint matrix
	   arma::vec A = ones(L);
	   double* pr_A = A.memptr();;
	   
	   // Non-Active Constraints Set
	   arma::ivec NASet(L); NASet.ones(); NASet = -NASet;
	   arma::ivec NAMask(p); NAMask.zeros();

	   int na;
	   arma::vec xRed(L);
	   arma::vec cRed(L);
	   arma::mat MRed(m,L);
	   double* pr_MRed = MRed.memptr();
	   arma::mat GRed(L,L);
	   arma::mat GRedinv(L, L);
	   double* pr_GRedinv = GRedinv.memptr();
	   arma::mat MTMRed(p,L);
	   double* pr_MTMRed = MTMRed.memptr();
	   
	  x.zeros();
	  x[0] = double(1.0);
	  // Non-Active Constraints Set
	  NASet[0] = 0;
	  NAMask[0] = 1;
	  na = 1;
	  xRed[0] = x[0];
	  cRed[0] = c[0];
	  cblas_dcopy(m, pr_M, 1, pr_MRed, 1);

	  // BLAS GRed = MRedT * MRed + lam2sq (na = 1 for now)
	  double coeff = cblas_ddot(m,pr_MRed,1,pr_MRed,1) + lam2sq;
	  GRed(0,0) = coeff;
	  GRedinv(0,0)= double(1.0) / GRed(0,0);
	  cblas_dcopy(p, pr_G, 1, pr_MTMRed, 1);


	   arma::vec Gplus = G*x + c;
	   double* pr_Gplus = Gplus.memptr();

	   arma::vec gRed(L);
	   double* pr_gRed = gRed.memptr();
	   arma::vec MRedxRed(m);

	   arma::vec GinvA(L);
	   double* pr_GinvA = GinvA.memptr();
	   arma::vec Ginvg(L);
	   double* pr_Ginvg = Ginvg.memptr();
	   arma::vec PRed(L);
	   double* pr_PRed = PRed.memptr();
	   arma::vec UB(L);
	   double* pr_UB = UB.memptr();
	   arma::vec UAiB(L);
	   double* pr_UAiB = UAiB.memptr();
	   // main loop active set
	   int iter = 0;
	   while(iter <= 100*p) {
		  ++iter;
		  // update of na, NASet, NAMask, xRed, cRed, Gplus, GRedinv, already done
		  // now update gRed, GinvA, Ginvg  (no need to update GRed)

		  // gRed
		  // gRed = Gplus[NASet]
		  for(int i = 0; i<na; ++i) {
			 gRed[i] = Gplus[NASet[i]];
		  }

		  // GinvA
		  // BLAS GinvA = GRedinv * A (ARed == A)
		  cblas_dsymv(CblasColMajor,CblasUpper,na,double(1.0),pr_GRedinv,L,pr_A,1,double(),pr_GinvA,1);
		  // Ginvg
		  // BLAS Ginvg = GRedinv * gRed
		  cblas_dsymv(CblasColMajor,CblasUpper,na,double(1.0),pr_GRedinv,L,pr_gRed,1,double(),pr_Ginvg,1);
		  double sGinvg = double();
		  double sGinvA = double();
		  for(int i = 0; i< na; ++i) {
			 sGinvg += Ginvg[i];
			 sGinvA += GinvA[i];
		  }
		  double lambdaS = sGinvg / sGinvA;
		  // BLAS PRed = GinvA * lambdaS - Ginvg
		  cblas_dcopy(na, pr_GinvA, 1, pr_PRed, 1);
		  cblas_dscal(na, lambdaS, pr_PRed, 1);
		  cblas_daxpy(na, double(-1.0), pr_Ginvg, 1, pr_PRed, 1);

		  double maxPRed = abs(PRed[0]);
		  for(int i = 0; i< na; ++i) {
			 if(abs(PRed[i])>maxPRed)
				maxPRed = abs(PRed[i]);
		  }
		  if(maxPRed < 1e-10) {
			 // P = 0, no advance possible
			 bool isOpt = true;
			 double lamMin = -epsilon;
			 int indexMin = -1;
			 for(int i = 0; i <p; ++i) {
				if(!NAMask[i] && Gplus[i]-lambdaS < lamMin) {
				   isOpt = false;
				   lamMin = Gplus[i]-lambdaS;
				   indexMin = i;
				}
			 }

			 if(isOpt) {
				// Got the optimal, STOP!				
				return(x);
			 } else {
				// Add one constraint
				NAMask[indexMin] = 1;
				NASet[na] = indexMin;
				xRed[na] = x[indexMin];
				cRed[na] = c[indexMin];
				// Gplus inchange
				// update MTMRed
				cblas_dcopy(p, pr_G+indexMin*p, 1, pr_MTMRed+na*p, 1);

				// update GRedinv
				// BLAS UB = MRed.double * M[:, indexMin]
				// Use G instead here
				for(int i = 0; i< na; ++i) {
				   UB[i] = G(NASet[i], indexMin);
				}
				// BLAS UC = M[:,indexMin].double* M[:, indexMin]
				double UC = G(indexMin, indexMin);
				// BLAS UAiB = GRedinv * UB
				cblas_dsymv(CblasColMajor,CblasUpper,na,double(1.0),pr_GRedinv,L,pr_UB,1,double(),pr_UAiB,1);
				double USi = double(1.0)/(UC - cblas_ddot(na,pr_UB,1,pr_UAiB,1));
				// GRedinv (restricted) += USi * UAiB*UAiB
				//replace cblas_syr(CblasColMajor,CblasUpper,na,USi, pr_UAiB, 1, pr_GRedinv, L);  
				cblas_dger(CblasColMajor,na,na,USi,pr_UAiB,1,pr_UAiB,1,pr_GRedinv,L);
				// copy -UAiB*USi, -UAiB.double*USi, USi to GRedinv
				cblas_dcopy(na, pr_UAiB, 1, pr_GRedinv+na*L, 1);
				cblas_dscal(na, -USi, pr_GRedinv+na*L,1);
				cblas_dcopy(na, pr_UAiB, 1, pr_GRedinv+na, L);
				cblas_dscal(na, -USi, pr_GRedinv+na,L);
				GRedinv(na,na) = USi;

				na += 1;
				assert(na <= L);
			 }
		  } else {
			 // P != 0, can advance
			 int indexMin = -1;
			 double alphaMin = double(1.0);
			 for(int i = 0; i < na; ++i) {
				if(PRed[i] < 0 && -xRed[i]/PRed[i] < alphaMin) {
				   indexMin = i;
				   alphaMin = -xRed[i]/PRed[i];
				}
			 }
			 // update x and Gplus
			 cblas_dscal(na, min(double(1.0), alphaMin), pr_PRed, 1);
			 for(int i = 0; i< na; ++i) {
				x[NASet[i]] += PRed[i];
				xRed[i] = x[NASet[i]];
			 }
			 // Gplus += MTMRed * (scaled PRed)
			 cblas_dgemv(CblasColMajor,CblasNoTrans,p,na,double(1.0),pr_MTMRed,p,pr_PRed,1,double(1.0),pr_Gplus,1);


			 // delete one constraint or not?
			 if(indexMin != -1) {
				// give true 0
				// x[NASet[indexMin]] = double();
				// delete one constraint
				NAMask[NASet[indexMin]] = 0;
				// downdate remove this -1;
				na -= 1;
				for(int i = indexMin; i<na; ++i) {
				   NASet[i] = NASet[i+1];
				   xRed[i] = xRed[i+1];
				   cRed[i] = cRed[i+1];
				}
				NASet[na] = -1;
				xRed[na] = double();
				cRed[na] = double();
				// PRed also
				PRed[na] = double();

				// downdate MTMRed
				for(int i = indexMin; i<na; ++i)
				   cblas_dcopy(p, pr_MTMRed+(i+1)*p, 1, pr_MTMRed+i*p, 1);

				// downdate GRedinv
				double UCi = double(1.0)/GRedinv(indexMin, indexMin); 
				// BLAS UB = GRedinv[ALL\indexMin,indexMin]
				cblas_dcopy(na+1, pr_GRedinv+indexMin*L, 1, pr_UB, 1);
				for(int i = indexMin; i<na; ++i)
				   UB[i] = UB[i+1];
				UB[na] = double();
				// get (GRedinv translated)
				// column first
				for(int i = indexMin; i<na; ++i)
				   cblas_dcopy(na+1, pr_GRedinv+(i+1)*L,1,pr_GRedinv+i*L,1);
				// row then
				for(int i = indexMin; i<na; ++i)
				   cblas_dcopy(na+1, pr_GRedinv+i+1, L, pr_GRedinv+i, L);

				// BLAS GRedinv = (GRedinv translated) - UB*UB.double*UCi 
				//replace cblas_syr(CblasColMajor,CblasUpper,na,-UCi, pr_UB, 1, pr_GRedinv, L);  
				cblas_dger(CblasColMajor,na,na,-UCi,pr_UB,1,pr_UB,1,pr_GRedinv,L);
			 }
		  }
	   }
	   return(x);
	}
	
	
	
	// min(|| AX - B ||) s.t. simplex constraint
	mat run_simplex_regression(mat &A, mat &B) {
		double lambda2 = 1e-5, epsilon = 1e-5;
		
		mat X = zeros(A.n_cols, B.n_cols);		
		if(A.n_cols < 1000) {
			double lam2sq =lambda2*lambda2;
			mat G = trans(A)*A + lam2sq;
			for(int i = 0; i < B.n_cols; i++) {
				vec b = B.col(i);
				X.col(i) = activeSetS_arma( A,  b, G, lambda2, epsilon);
			}			
		} else {
			for(int i = 0; i < B.n_cols; i++) {
				vec b = B.col(i);
				X.col(i) = activeSet_arma( A,  b, lambda2, epsilon);
			}
		}
		
		
		X = clamp(X, 0, 1);
		X = normalise(X, 1);

		return(X);
	}	
	
}

