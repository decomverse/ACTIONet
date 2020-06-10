#include "arch.h"
#define NEW_VERSION

namespace ACTIONet {
	void arch_dense(SPAMS_Matrix<double>& X, SPAMS_Matrix<double>& Z0, SPAMS_Matrix<double>& Z,  SPAMS_Matrix<double>& A, SPAMS_Matrix<double>& B, int I1, int I2, double lambda2, double epsilon, bool computeXtX) {	
	   const int m = X.m();
	   const int n = X.n();
	   const int p = Z0.n();
	   Z.copy(Z0);
	   SPAMS_Matrix<double> AlphaT(p,n);
	   SPAMS_Matrix<double> BetaT(n,p);
	   double RSS = -1.0;
	   Vector<double> refColZ;
	   Vector<double> copRowAlphaT;
	   Vector<double> refColBetaT;
	   SPAMS_Matrix<double> matRSS(m,n);
	   Vector<double> vBarre(m);
	   Vector<double> norms;
	   cout.precision(8);


	   for(int t=0; t<I2; ++t) {
		   printf("%d\n", t);
		  SPAMS_Matrix<double> G;
		  if (computeXtX) {
			 Z.XtX(G);
			 G.addDiag(lambda2*lambda2);
		  }
		  // step 1: fix Z to compute Alpha
		  for(int i=0; i<n; ++i) {
			 Vector<double> refColX;
			 Vector<double> refColAlphaT;
			 X.refCol(i,refColX);
			 AlphaT.refCol(i, refColAlphaT);
			 if (computeXtX) {
				activeSetS<double>(Z,refColX, refColAlphaT, G, lambda2, epsilon);
			 } else {
				activeSet<double>(Z,refColX, refColAlphaT, lambda2, epsilon);
			 }
		  }
		  // step 2: fix Alpha, fix all but one to compute Zi
	#ifdef NEW_VERSION
		  // new version
		  Vector<double> refColX;
		  Vector<double> tmp;
		  matRSS.copy(X);
		  Z.mult(AlphaT, matRSS, false, false, double(-1.0), double(1.0));
		  for(int l=0; l<p; ++l) {
			 AlphaT.copyRow(l, copRowAlphaT);
			 double sumAsq =  copRowAlphaT.nrm2sq();
			 Z.refCol(l, refColZ);
			 // matRSS = X- Z*AlphaT
			 if(sumAsq < double(10e-8)) {
				// singular
				matRSS.norm_2_cols(norms);
				int k = norms.max();
				X.refCol(k, refColX);
				refColZ.copy(refColX);
			 } else {
				//matRSS.rank1Update(refColZ, copRowAlphaT);
				matRSS.mult(copRowAlphaT, vBarre, 1/sumAsq, double());
				vBarre.add(refColZ);
				tmp.copy(refColZ);
				// least square to get Beta
				BetaT.refCol(l, refColBetaT);
				activeSet<double>(X, vBarre, refColBetaT, lambda2, epsilon);
				X.mult(refColBetaT, refColZ);
				tmp.sub(refColZ);
				matRSS.rank1Update(tmp, copRowAlphaT);
			 }
		  }
	#else
		  // end new version

		  Vector<double> refColX;
		  for(int l=0; l<p; ++l) {
			 AlphaT.copyRow(l, copRowAlphaT);
			 double sumAsq =  copRowAlphaT.nrm2sq();
			 Z.refCol(l, refColZ);
			 // matRSS = X- Z*AlphaT
			 matRSS.copy(X);
			 Z.mult(AlphaT, matRSS, false, false, double(-1.0), double(1.0));
			 if(sumAsq < double(10e-8)) {
				// singular
				matRSS.norm_2_cols(norms);
				int k = norms.max();
				X.refCol(k, refColX);
				refColZ.copy(refColX);
			 } else {
				matRSS.rank1Update(refColZ, copRowAlphaT);
				matRSS.mult(copRowAlphaT, vBarre, 1/sumAsq, double());
				// least square to get Beta
				BetaT.refCol(l, refColBetaT);
				activeSet<double>(X, vBarre, refColBetaT, lambda2, epsilon);
				X.mult(refColBetaT, refColZ);
			 }
		  }

		  matRSS.copy(X);
		  Z.mult(AlphaT, matRSS, false, false, double(-1.0), double(1.0));
	#endif
		  RSS = matRSS.normFsq();
		  
	#ifdef DEBUG      
		  cout << "RSS AS = " << RSS << endl;
	#endif      
		  flush(cout);
	   }
		  
	   memcpy(A._X, AlphaT._X, p*n*sizeof(double));
	   memcpy(B._X, BetaT._X, n*p*sizeof(double));
	}

	template <typename T>
	void arch(const SPAMS_Matrix<T>& X, const SPAMS_Matrix<T>& Z0, SPAMS_Matrix<T>& Z,  SpSPAMS_Matrix<T>& A, SpSPAMS_Matrix<T>& B, const int I1, const int I2, const T lambda2, const T epsilon, const bool computeXtX) {
	   const int m = X.m();
	   const int n = X.n();
	   const int p = Z0.n();
	   Z.copy(Z0);
	   SPAMS_Matrix<T> AlphaT(p,n);
	   SPAMS_Matrix<T> BetaT(n,p);
	   T RSS = -1.0;
	   Vector<T> refColZ;
	   Vector<T> copRowAlphaT;
	   Vector<T> refColBetaT;
	   SPAMS_Matrix<T> matRSS(m,n);
	   Vector<T> vBarre(m);
	   Vector<T> norms;
	   cout.precision(8);

	   for(int t=0; t<I1; ++t) {
		  // step 1: fix Z to compute Alpha
		  for(int i=0; i<n; ++i) {
			 Vector<T> refColX;
			 Vector<T> refColAlphaT;
			 X.refCol(i,refColX);
			 AlphaT.refCol(i, refColAlphaT);
			 gpFISTAFor(Z,refColX, refColAlphaT, T(1.0), T(1.0/0.7), 50, true);
		  }
		  // step 2: fix Alpha, fix all but one to compute Zi
		  Vector<T> refColX;
		  for(int l=0; l<p; ++l) {
			 AlphaT.copyRow(l, copRowAlphaT);
			 T sumAsq =  copRowAlphaT.nrm2sq();
			 Z.refCol(l, refColZ);
			 // matRSS = X- Z*AlphaT
			 matRSS.copy(X);
			 Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
			 if(sumAsq < T(10e-8)) {
				// singular
				matRSS.norm_2_cols(norms);
				int k = norms.max();
				X.refCol(k, refColX);
				refColZ.copy(refColX);
			 } else {
				matRSS.rank1Update(refColZ, copRowAlphaT);
				matRSS.mult(copRowAlphaT, vBarre, 1/sumAsq, T());
				// least square to get Beta
				BetaT.refCol(l, refColBetaT);
				gpFISTAFor(X, vBarre, refColBetaT, T(1.0), T(1.0/0.7), 50, true);
				X.mult(refColBetaT, refColZ);
			 }
		  }

		  matRSS.copy(X);
		  Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
		  RSS = matRSS.normFsq();
		  cout << "RSS FISTA = " << RSS << endl;
		  flush(cout);
	   }  

	   for(int t=0; t<I2; ++t) {
		  SPAMS_Matrix<T> G;
		  if (computeXtX) {
			 Z.XtX(G);
			 G.addDiag(lambda2*lambda2);
		  }
		  // step 1: fix Z to compute Alpha
		  for(int i=0; i<n; ++i) {
			 Vector<T> refColX;
			 Vector<T> refColAlphaT;
			 X.refCol(i,refColX);
			 AlphaT.refCol(i, refColAlphaT);
			 if (computeXtX) {
				activeSetS<T>(Z,refColX, refColAlphaT, G, lambda2, epsilon);
			 } else {
				activeSet<T>(Z,refColX, refColAlphaT, lambda2, epsilon);
			 }
		  }
		  // step 2: fix Alpha, fix all but one to compute Zi
	#ifdef NEW_VERSION
		  // new version
		  Vector<T> refColX;
		  Vector<T> tmp;
		  matRSS.copy(X);
		  Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
		  for(int l=0; l<p; ++l) {
			 AlphaT.copyRow(l, copRowAlphaT);
			 T sumAsq =  copRowAlphaT.nrm2sq();
			 Z.refCol(l, refColZ);
			 // matRSS = X- Z*AlphaT
			 if(sumAsq < T(10e-8)) {
				// singular
				matRSS.norm_2_cols(norms);
				int k = norms.max();
				X.refCol(k, refColX);
				refColZ.copy(refColX);
			 } else {
				//matRSS.rank1Update(refColZ, copRowAlphaT);
				matRSS.mult(copRowAlphaT, vBarre, 1/sumAsq, T());
				vBarre.add(refColZ);
				tmp.copy(refColZ);
				// least square to get Beta
				BetaT.refCol(l, refColBetaT);
				activeSet<T>(X, vBarre, refColBetaT, lambda2, epsilon);
				X.mult(refColBetaT, refColZ);
				tmp.sub(refColZ);
				matRSS.rank1Update(tmp, copRowAlphaT);
			 }
		  }
	#else
		  // end new version

		  Vector<T> refColX;
		  for(int l=0; l<p; ++l) {
			 AlphaT.copyRow(l, copRowAlphaT);
			 T sumAsq =  copRowAlphaT.nrm2sq();
			 Z.refCol(l, refColZ);
			 // matRSS = X- Z*AlphaT
			 matRSS.copy(X);
			 Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
			 if(sumAsq < T(10e-8)) {
				// singular
				matRSS.norm_2_cols(norms);
				int k = norms.max();
				X.refCol(k, refColX);
				refColZ.copy(refColX);
			 } else {
				matRSS.rank1Update(refColZ, copRowAlphaT);
				matRSS.mult(copRowAlphaT, vBarre, 1/sumAsq, T());
				// least square to get Beta
				BetaT.refCol(l, refColBetaT);
				activeSet<T>(X, vBarre, refColBetaT, lambda2, epsilon);
				X.mult(refColBetaT, refColZ);
			 }
		  }

		  matRSS.copy(X);
		  Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
	#endif
		  RSS = matRSS.normFsq();
		  cout << "RSS AS = " << RSS << endl;
		  flush(cout);
	   }
	   AlphaT.toSparse(A);
	   BetaT.toSparse(B);
	}

	template <typename T>
	void archRobust(const SPAMS_Matrix<T>& X, const SPAMS_Matrix<T>& Z0, SPAMS_Matrix<T>& Z,  SpSPAMS_Matrix<T>& A, SpSPAMS_Matrix<T>& B, const int I1, const int I2, const T lambda2, const T epsilon, const T epsilon2, const bool computeXtX) {
	   const int m = X.m();
	   const int n = X.n();
	   const int p = Z0.n();
	   Z.copy(Z0);
	   SPAMS_Matrix<T> AlphaT(p,n);
	   SPAMS_Matrix<T> BetaT(n,p);

	   T RSN = -1.0;
	   Vector<T> refColZ;
	   Vector<T> copRowAlphaT;
	   Vector<T> refColBetaT;
	   SPAMS_Matrix<T> matRSS(m,n);
	   Vector<T> vBarre(m);
	   Vector<T> norms;
	   cout.precision(8);

	   for(int t=0; t<I1; ++t) {
		  // step 1: fix Z to compute Alpha
		  for(int i=0; i<n; ++i) {
			 Vector<T> refColX;
			 Vector<T> refColAlphaT;
			 X.refCol(i,refColX);
			 AlphaT.refCol(i, refColAlphaT);
			 gpFISTAFor(Z, refColX, refColAlphaT, T(1.0), T(1.0/0.7), 10, true);
		  }
		  // update scale factors
		  matRSS.copy(X);
		  Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
		  matRSS.norm_2_cols(norms);
		  norms.thrsmax(epsilon2);
		  norms.Sqrt();
		  Vector<T> refColX;
		  // step 2: fix Alpha, fix all but one to compute Zi
		  for(int l=0; l<p; ++l) {
			 Z.refCol(l, refColZ);

			 AlphaT.copyRow(l, copRowAlphaT);
			 copRowAlphaT.div(norms);
			 T sumAsq =  copRowAlphaT.nrm2sq();

			 matRSS.copy(X);
			 Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));

			 if(sumAsq < T(10e-8)) {
				// singular
				matRSS.norm_2_cols(norms);
				int k = norms.max();
				X.refCol(k, refColX);
				refColZ.copy(refColX);
			 } else {
				// absorbe the weights by rowAlphaT
				copRowAlphaT.div(norms);
				matRSS.mult(copRowAlphaT, vBarre, 1/sumAsq, T());
				vBarre.add(refColZ);
				// least square to get Beta
				BetaT.refCol(l, refColBetaT);
				gpFISTAFor(X, vBarre, refColBetaT, T(1.0), T(1.0/0.7), 10, true); 
				X.mult(refColBetaT, refColZ);
			 }
		  }

		  matRSS.copy(X);
		  Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
		  matRSS.norm_2_cols(norms);
		  for (int i=0; i<norms.n(); ++i)
			 if (norms[i] <= epsilon2)
				norms[i]=norms[i]*norms[i]/(2*epsilon2) + epsilon2/2;
		  RSN = norms.sum();
		  cout << "RSN FISTA= " << RSN << endl;
		  flush(cout);
	   }

	   for(int t=0; t<I2; ++t) {
		  SPAMS_Matrix<T> G;
		  if (computeXtX) {
			 Z.XtX(G);
			 G.addDiag(lambda2*lambda2);
		  }
		  // step 1: fix Z to compute Alpha
		  for(int i=0; i<n; ++i) {
			 Vector<T> refColX;
			 Vector<T> refColAlphaT;
			 X.refCol(i,refColX);
			 AlphaT.refCol(i, refColAlphaT);
			 if (computeXtX) {
				activeSetS<T>(Z,refColX, refColAlphaT, G, lambda2, epsilon);
			 } else {
				activeSet<T>(Z,refColX, refColAlphaT, lambda2, epsilon);
			 }
		  }
		  // update scale factors
	#ifndef NEW_VERSION
		  matRSS.copy(X);
		  Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
		  matRSS.norm_2_cols(norms);
		  norms.thrsmax(epsilon2);
		  norms.Sqrt();
		  // step 2: fix Alpha, fix all but one to compute Zi
		  Vector<T> refColX;
		  for(int l=0; l<p; ++l) {
			 Z.refCol(l, refColZ);

			 AlphaT.copyRow(l, copRowAlphaT);
			 copRowAlphaT.div(norms);
			 T sumAsq =  copRowAlphaT.nrm2sq();

			 matRSS.copy(X);
			 Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));

			 if(sumAsq < T(10e-8)) {
				// singular
				matRSS.norm_2_cols(norms);
				int k = norms.max();
				X.refCol(k, refColX);
				refColZ.copy(refColX);
			 } else {
				// absorbe the weights by rowAlphaT
				copRowAlphaT.div(norms);
				matRSS.mult(copRowAlphaT, vBarre, 1/sumAsq, T());
				vBarre.add(refColZ);
				// least square to get Beta
				BetaT.refCol(l, refColBetaT);
				activeSet<T>(X, vBarre, refColBetaT, lambda2, epsilon);
				X.mult(refColBetaT, refColZ);
			 }
		  }

		  matRSS.copy(X);
		  Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
	#else
		  /// new version
		  Vector<T> refColX;
		  Vector<T> tmp;
		  Vector<T> tmp2;
		  matRSS.copy(X);
		  Z.mult(AlphaT, matRSS, false, false, T(-1.0), T(1.0));
		  matRSS.norm_2_cols(norms);
		  norms.thrsmax(epsilon2);
		  norms.Sqrt();
		  for(int l=0; l<p; ++l) {
			 Z.refCol(l, refColZ);

			 AlphaT.copyRow(l, copRowAlphaT);
			 tmp2.copy(copRowAlphaT);
			 copRowAlphaT.div(norms);
			 T sumAsq =  copRowAlphaT.nrm2sq();

			 if(sumAsq < T(10e-8)) {
				// singular
				matRSS.norm_2_cols(tmp);
				int k = tmp.max();
				X.refCol(k, refColX);
				refColZ.copy(refColX);
			 } else {
				// absorbe the weights by rowAlphaT
				copRowAlphaT.div(norms);
				matRSS.mult(copRowAlphaT, vBarre, 1/sumAsq, T());
				vBarre.add(refColZ);
				tmp.copy(refColZ);
				// least square to get Beta
				BetaT.refCol(l, refColBetaT);
				activeSet<T>(X, vBarre, refColBetaT, lambda2, epsilon);
				X.mult(refColBetaT, refColZ);
				tmp.sub(refColZ);
				matRSS.rank1Update(tmp,tmp2);
			 }
		  }
	#endif
		  /// end new version
		  
		  matRSS.norm_2_cols(norms);
		  for (int i=0; i<norms.n(); ++i)
			 if (norms[i] <= epsilon2)
				norms[i]=norms[i]*norms[i]/(2*epsilon2) + epsilon2/2;
		  RSN = norms.sum();
		  cout << "RSN AS= " << RSN << endl;
		  flush(cout);
	   }
	   AlphaT.toSparse(A);
	   BetaT.toSparse(B);
	}

	template <typename T>
	void archetypalAnalysis(const SPAMS_Matrix<T>& X, const SPAMS_Matrix<T>& Z0, SPAMS_Matrix<T>& Z, SpSPAMS_Matrix<T>& A, SpSPAMS_Matrix<T>& B, const bool robust, const T epsilon2, const bool computeXtX, const int stepsFISTA, const int stepsAS, const int numThreads) {
	   init_omp(numThreads);
	   const T epsilon = 1e-5;
	   const T lambda2 = 1e-5;
	   if (!robust) {
		  arch(X, Z0, Z, A, B, stepsFISTA, stepsAS, epsilon,lambda2,computeXtX);
	   } else {
		  archRobust(X, Z0, Z, A, B, stepsFISTA, stepsAS, epsilon,lambda2,epsilon2,computeXtX);
	   }
	}

	template <typename T>
	void archetypalAnalysis(const SPAMS_Matrix<T>& X, SPAMS_Matrix<T>& Z, SpSPAMS_Matrix<T>& A, SpSPAMS_Matrix<T>& B, const bool robust, const T epsilon2, const bool computeXtX, const int stepsFISTA, const int stepsAS, const bool randominit, const int numThreads) {

	   const int m = X.m();
	   const int n = X.n();
	   const int p = Z.n();
	   SPAMS_Matrix<T> Z0(m,p);
	   Vector<T> refColZ0;
	   Vector<T> refColX;
	   if(!randominit) {
		  for(int i=0; i<p; i++) {
			 X.refCol(i%n, refColX);
			 Z0.refCol(i%n, refColZ0);
			 refColZ0.copy(refColX);
		  }
	   } else {
		  srandom(0);
		  for(int i=0; i<p; i++) {
			 int k = random() % n;
			 X.refCol(k, refColX);
			 Z0.refCol(i, refColZ0);
			 refColZ0.copy(refColX);
		  }
	   }
	   archetypalAnalysis(X, Z0, Z, A, B, robust, epsilon2, computeXtX, stepsFISTA, stepsAS,numThreads);
	}

	template <typename T>
	void decompSimplex(const SPAMS_Matrix<T>& X, const SPAMS_Matrix<T>& Z, SpSPAMS_Matrix<T>& alpha, const bool computeZtZ, const int numThreads) {
	   init_omp(numThreads);
	   const int n = X.n();
	   const int p = Z.n();
	   SPAMS_Matrix<T> AlphaT(p,n);
	   int i;
	   if(computeZtZ) {
		  SPAMS_Matrix<T> G;
		  Z.XtX(G);
		  T lambda2 = 1e-5;
		  G.addDiag(lambda2*lambda2);
		  for(i=0; i<n; ++i) {
			 Vector<T> refColX;
			 Vector<T> refColAlphaT;
			 X.refCol(i,refColX);
			 AlphaT.refCol(i, refColAlphaT);
			 activeSetS(Z,refColX, refColAlphaT, G);
		  }
		  AlphaT.toSparse(alpha);
	   } else {
		  for(i=0; i<n; ++i) {
			 Vector<T> refColX;
			 Vector<T> refColAlphaT;
			 X.refCol(i,refColX);
			 AlphaT.refCol(i, refColAlphaT);
			 activeSet(Z,refColX, refColAlphaT);
		  }
		  AlphaT.toSparse(alpha);
	   }
	}


	void decompSimplex_dense (const SPAMS_Matrix<double>& X, const SPAMS_Matrix<double>& Z, SPAMS_Matrix<double>& alpha, const bool computeZtZ, const int numThreads) {
	   init_omp(numThreads);
	   const int n = X.n();
	   const int p = Z.n();
	   SPAMS_Matrix<double> AlphaT(p,n);
	   int i;
	   if(computeZtZ) {
		  SPAMS_Matrix<double> G;
		  Z.XtX(G);
		  double lambda2 = 1e-5;
		  G.addDiag(lambda2*lambda2);
		  for(i=0; i<n; ++i) {
			 Vector<double> refColX;
			 Vector<double> refColAlphaT;
			 X.refCol(i,refColX);
			 AlphaT.refCol(i, refColAlphaT);
			 activeSetS(Z,refColX, refColAlphaT, G);
		  }
	   } else {
		  for(i=0; i<n; ++i) {
			 Vector<double> refColX;
			 Vector<double> refColAlphaT;
			 X.refCol(i,refColX);
			 AlphaT.refCol(i, refColAlphaT);
			 activeSet(Z,refColX, refColAlphaT);
		  }
	   }
	   
	   memcpy(alpha._X, AlphaT._X, n*p*sizeof(double));
	}


	void AA (double *A_ptr, int A_rows, int A_cols, double *W0_ptr, int W0_cols, double *C_ptr, double *H_ptr, int stepsAS = 50) {
		SPAMS_Matrix<double> C_spam;		
		C_spam._X = C_ptr;
		C_spam._externAlloc = true;	
		
		SPAMS_Matrix<double> H_spam;		
		H_spam._X = H_ptr;
		H_spam._externAlloc = true;

		SPAMS_Matrix<double> A_spam(A_ptr, A_rows, A_cols);
		SPAMS_Matrix<double> W0_spam(W0_ptr, A_rows, W0_cols);

		SPAMS_Matrix<double> W_spam(A_rows, W0_cols);
		
		double lambda2 = 1e-5;
		double epsilon = 1e-5;
		int stepsFISTA = 3;
			
		arch_dense(A_spam, W0_spam, W_spam, H_spam, C_spam, stepsFISTA, stepsAS, epsilon,lambda2, false); //arch_dense(A, W0, Z, A, B....)
	}


	void simplexRegression(double *A_ptr, int A_cols, double *B_ptr, int B_rows, int B_cols, double *X_ptr) { // min(|| AX - B ||) s.t. simplex constraint
		double lambda2 = 1e-5;
		double epsilon = 1e-5;
		
		SPAMS_Matrix<double> Z(A_ptr, B_rows, A_cols); 
		SPAMS_Matrix<double> X(B_ptr, B_rows, B_cols);	
		SPAMS_Matrix<double> AlphaT(X_ptr, A_cols, B_cols);
					
		Vector<double> refColX;
		Vector<double> refColAlphaT;	
		for(register int i = 0; i < B_cols; i++) {
			X.refCol(i, refColX);
			AlphaT.refCol(i, refColAlphaT);

			activeSet<double>(Z, refColX, refColAlphaT, lambda2, epsilon);
		}			
	}
}
