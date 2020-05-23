#include "ACTIONet.h"


namespace ACTIONet {
	mat sampleUnif(int l, int m, double a, double b, int seed) {
		std::default_random_engine gen (seed);	
		std::uniform_real_distribution<double> unif(a, b);
		
		
		mat R(l, m);
		for (int j = 0; j < m; j++) {
			for(int i = 0; i < l; i++) {
				R(i, j) = unif(gen);
			}
		}
		return R;
	}

	void gram_schmidt(mat& A) {
		for(uword i = 0; i < A.n_cols; ++i) {
			for(uword j = 0; j < i; ++j) {
				double r = dot(A.col(i), A.col(j));
				A.col(i) -= r * A.col(j);
			}
			
			double col_norm = norm(A.col(i), 2);
			
			if(col_norm < 1E-4) {
				for(uword k = i; k < A.n_cols; ++k)
					A.col(k).zeros();

				return;
			}
			A.col(i) /= col_norm;
		}
	}
	
	field<mat> eigSVD(mat A) {
		int n = A.n_cols;
		mat B = trans(A)*A;
				
		vec d;
		mat V;
		eig_sym( d, V, B );		
		d = sqrt(d);
		
		// Compute U
		sp_mat S(n, n);
		S.diag() = 1 / d; 
		mat U = (S*trans(V))*trans(A);
		U = trans(U);
		
		field<mat> out(3);
		
		out(0) = U;
		out(1) = d;
		out(2) = V;
		
		return(out);
	}	
	
	mat zscore(mat A) {	
		rowvec mu = mean(A, 0);
		rowvec sigma = stddev(A, 0);

		
		for(int j = 0; j < A.n_cols; j++) {
			A.col(j) = (A.col(j) - mu(j)) / sigma(j);
		}
		
		return A;
	}	
	
	
}
