#include "ACTIONet.h"

// From: N Halko, P. G Martinsson, and J. A Tropp. Finding structure with randomness: Probabilistic algorithms for constructing approximate matrix decompositions. Siam Review, 53(2):217-288, 2011.
namespace ACTIONet {
	field<mat> HalkoSVD(sp_mat &A, int dim, int iters, int seed = 0) {	
		field<mat> results(3);

		int m = A.n_rows;
		int n = A.n_cols;
		int l = dim + 2;		
				
		vec s;
		mat R, Q;
		mat U, V, X;
		
		printf("\t\tRunning randomized SVD. Matrix size: %d x %d\n", A.n_rows, A.n_cols); fflush(stdout);
		
		if (m < n) {
			//R = stats::runif<arma::mat>(l, m, -1.0, 1.0, seed);
			R = sampleUnif(l, m, -1.0, 1.0, 0);

			sp_mat At = A.t();
			Q = At*R.t(); 
		}
		else {
			//R = stats::runif<arma::mat>(n, l, -1.0, 1.0, seed);
			R = sampleUnif(n, l, -1.0, 1.0, 0);
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
		
		printf("\t\tRunning randomized SVD. Matrix size: %d x %d\n", A.n_rows, A.n_cols); fflush(stdout);
		
		if (m < n) {
			//R = stats::runif<arma::mat>(l, m, -1.0, 1.0, seed);
			R = sampleUnif(l, m, -1.0, 1.0, 0);

			mat At = A.t();
			Q = At*R.t(); 
		}
		else {
			//R = stats::runif<arma::mat>(n, l, -1.0, 1.0, seed);
			R = sampleUnif(n, l, -1.0, 1.0, 0);
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
