#include "ACTIONet.h"

namespace ACTIONet {
	
	//****************************************************************************************************************************************************************************
	// From: Xu Feng, Yuyang Xie, and Yaohang Li, "Fast Randomzied SVD for Sparse Data," in Proc. the 10th Asian Conference on Machine Learning (ACML), Beijing, China, Nov. 2018.
	//****************************************************************************************************************************************************************************
	field<mat> FengSVD(sp_mat &A, int dim, int iters, int seed = 0) {	
		int s = 5;

		int m = A.n_rows;
		int n = A.n_cols;
		
		printf("\t\tRunning randomized SVD. Matrix size: %d x %d (# iters = %d)\n", m, n, iters); fflush(stdout);
				
		vec S;
		mat Q, L, U, V;
		field<mat> SVD_out;
		
		if(m < n) {
			printf("\t\t\tInitializing SVD (mode 1) ... ");
			//arma_rng::set_seed(seed);			
			//Q = randn( n, dim+s );
			Q = sampleUnif(n, dim+s, 0.0, 1.0, seed);
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
			Q = sampleUnif(m, dim+s, 0.0, 1.0, seed);
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
		
		printf("\t\tRunning randomized SVD (full matrix version). Matrix size: %d x %d (# iters = %d)\n", m, n, iters); fflush(stdout);
				
		vec S;
		mat Q, L, U, V;
		field<mat> SVD_out;
		
		if(m < n) {
			printf("\t\t\tInitializing SVD (mode 1) ... ");
			//arma_rng::set_seed(seed);			
			//Q = randn( n, dim+s );
			Q = sampleUnif(n, dim+s, 0.0, 1.0, seed);
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
			Q = sampleUnif(m, dim+s, 0.0, 1.0, seed);
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
