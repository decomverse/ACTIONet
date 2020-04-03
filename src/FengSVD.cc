#include "ACTIONet.h"

// From: Xu Feng, Yuyang Xie, and Yaohang Li, "Fast Randomzied SVD for Sparse Data," in Proc. the 10th Asian Conference on Machine Learning (ACML), Beijing, China, Nov. 2018.
namespace ACTIONet {
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
}
