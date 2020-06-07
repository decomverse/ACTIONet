#include "ACTIONet.h"


namespace ACTIONet {
	field<mat> PCA2SVD(sp_mat &S, field<mat> PCA_results) {	
		int n = S.n_rows;
				
		printf("\tRetriving the PCA results\n"); fflush(stdout);
		mat U = PCA_results(0);
		vec s = PCA_results(1);
		mat V = PCA_results(2);
		
		int dim = U.n_cols;


		mat A = ones(S.n_rows, 1);
		mat B = mat(trans(mean(S, 0)));
		
	
		printf("\tPCA => SVD ..."); fflush(stdout);
		vec s_prime;
		mat U_prime, V_prime;
			
		mat M = U.t()*A; 
		mat A_ortho_proj = A - U*M;   
		mat P = A_ortho_proj;// = orth(A_ortho_proj);
		gram_schmidt(P);		
		mat R_P = P.t()*A_ortho_proj;
		
		
		mat N = V.t()*B; 
		mat B_ortho_proj = B - V*N; 
		mat Q = B_ortho_proj; //orth(B_ortho_proj); 
		gram_schmidt(Q);
		mat R_Q = Q.t()*B_ortho_proj;	
		
		mat K1 = zeros(s.n_elem+A.n_cols, s.n_elem+A.n_cols);
		for(int i = 0; i < s.n_elem; i++) {
			K1(i, i) = s(i);
		}

		mat K2 = join_vert(M, R_P)*trans(join_vert(N, R_Q));

		
		mat K = K1 + K2;

		svd( U_prime, s_prime, V_prime, K );
		
		mat U_updated = join_horiz(U, P)*U_prime;
		mat V_updated = join_horiz(V, Q)*V_prime;
		printf("done.\n"); fflush(stdout);
		
		field<mat> SVD_results(3);
		SVD_results(0) = U_updated.cols(0, dim-1);
		SVD_results(1) = s_prime(span(0, dim-1));
		SVD_results(2) = V_updated.cols(0, dim-1);
		
		return SVD_results;					
	}	


	field<mat> PCA2SVD(mat &S, field<mat> PCA_results) {	
		int n = S.n_rows;
				
		printf("\tRetriving the PCA results\n"); fflush(stdout);
		mat U = PCA_results(0);
		vec s = PCA_results(1);
		mat V = PCA_results(2);
		
		int dim = U.n_cols;


		mat A = ones(S.n_rows, 1);
		mat B = mat(trans(mean(S, 0)));
		
	
		printf("\tPCA => SVD ..."); fflush(stdout);
		vec s_prime;
		mat U_prime, V_prime;
			
		mat M = U.t()*A; 
		mat A_ortho_proj = A - U*M;   
		mat P = A_ortho_proj;// = orth(A_ortho_proj);
		gram_schmidt(P);		
		mat R_P = P.t()*A_ortho_proj;
		
		
		mat N = V.t()*B; 
		mat B_ortho_proj = B - V*N; 
		mat Q = B_ortho_proj; //orth(B_ortho_proj); 
		gram_schmidt(Q);
		mat R_Q = Q.t()*B_ortho_proj;	
		
		mat K1 = zeros(s.n_elem+A.n_cols, s.n_elem+A.n_cols);
		for(int i = 0; i < s.n_elem; i++) {
			K1(i, i) = s(i);
		}

		mat K2 = join_vert(M, R_P)*trans(join_vert(N, R_Q));

		
		mat K = K1 + K2;

		svd( U_prime, s_prime, V_prime, K );
		
		mat U_updated = join_horiz(U, P)*U_prime;
		mat V_updated = join_horiz(V, Q)*V_prime;
		printf("done.\n"); fflush(stdout);
		
		field<mat> SVD_results(3);
		SVD_results(0) = U_updated.cols(0, dim-1);
		SVD_results(1) = s_prime(span(0, dim-1));
		SVD_results(2) = V_updated.cols(0, dim-1);
		
		return SVD_results;					
	}	




	field<mat> SVD2PCA(sp_mat &S, field<mat> SVD_results) {	
		int n = S.n_rows;
				
		printf("\tRetriving the PCA results\n"); fflush(stdout);
		mat U = SVD_results(0);
		vec s = SVD_results(1);
		mat V = SVD_results(2);
		
		int dim = U.n_cols;


		mat A = ones(S.n_rows, 1);
		mat B = -mat(trans(mean(S, 0)));
		
	
		printf("\tSVD => PCA ..."); fflush(stdout);
		vec s_prime;
		mat U_prime, V_prime;
			
		mat M = U.t()*A; 
		mat A_ortho_proj = A - U*M;   
		mat P = A_ortho_proj;// = orth(A_ortho_proj);
		gram_schmidt(P);		
		mat R_P = P.t()*A_ortho_proj;
		
		
		mat N = V.t()*B; 
		mat B_ortho_proj = B - V*N; 
		mat Q = B_ortho_proj; //orth(B_ortho_proj); 
		gram_schmidt(Q);
		mat R_Q = Q.t()*B_ortho_proj;	
		
		mat K1 = zeros(s.n_elem+A.n_cols, s.n_elem+A.n_cols);
		for(int i = 0; i < s.n_elem; i++) {
			K1(i, i) = s(i);
		}

		mat K2 = join_vert(M, R_P)*trans(join_vert(N, R_Q));

		
		mat K = K1 + K2;

		svd( U_prime, s_prime, V_prime, K );
		
		mat U_updated = join_horiz(U, P)*U_prime;
		mat V_updated = join_horiz(V, Q)*V_prime;
		printf("done.\n"); fflush(stdout);
		
		field<mat> PCA_results(3);
		PCA_results(0) = U_updated.cols(0, dim-1);
		PCA_results(1) = s_prime(span(0, dim-1));
		PCA_results(2) = V_updated.cols(0, dim-1);
		
		return PCA_results;					
	}	



	field<mat> SVD2PCA(mat &S, field<mat> SVD_results) {	
		int n = S.n_rows;
				
		printf("\tRetriving the PCA results\n"); fflush(stdout);
		mat U = SVD_results(0);
		vec s = SVD_results(1);
		mat V = SVD_results(2);
		
		int dim = U.n_cols;


		mat A = ones(S.n_rows, 1);
		mat B = -mat(trans(mean(S, 0)));
		
	
		printf("\tSVD => PCA ..."); fflush(stdout);
		vec s_prime;
		mat U_prime, V_prime;
			
		mat M = U.t()*A; 
		mat A_ortho_proj = A - U*M;   
		mat P = A_ortho_proj;// = orth(A_ortho_proj);
		gram_schmidt(P);		
		mat R_P = P.t()*A_ortho_proj;
		
		
		mat N = V.t()*B; 
		mat B_ortho_proj = B - V*N; 
		mat Q = B_ortho_proj; //orth(B_ortho_proj); 
		gram_schmidt(Q);
		mat R_Q = Q.t()*B_ortho_proj;	
		
		mat K1 = zeros(s.n_elem+A.n_cols, s.n_elem+A.n_cols);
		for(int i = 0; i < s.n_elem; i++) {
			K1(i, i) = s(i);
		}

		mat K2 = join_vert(M, R_P)*trans(join_vert(N, R_Q));

		
		mat K = K1 + K2;

		svd( U_prime, s_prime, V_prime, K );
		
		mat U_updated = join_horiz(U, P)*U_prime;
		mat V_updated = join_horiz(V, Q)*V_prime;
		printf("done.\n"); fflush(stdout);
		
		field<mat> PCA_results(3);
		PCA_results(0) = U_updated.cols(0, dim-1);
		PCA_results(1) = s_prime(span(0, dim-1));
		PCA_results(2) = V_updated.cols(0, dim-1);
		
		return PCA_results;					
	}	

	
	ReducedKernel ACTION_reduction(sp_mat &S, int dim, int iter = 5, int seed = 0, int SVD_algorithm = HALKO_ALG, bool prenormalize = false) {			
		int n = S.n_rows;
		if(prenormalize)
			S = normalise(S, 2);    

		printf("\tComputing reduced ACTION kernel. Input matrix size: %d x %d\n", S.n_rows, S.n_cols); fflush(stdout);
		
		// Update 1: Orthogonalize columns w.r.t. background (mean)
		vec mu = vec(mean(S, 1));
		vec v = mu / norm(mu, 2);
		vec a1 = v;
		vec b1 = -trans(S)*v;
		
		// Update 2: Center columns of orthogonalized matrix before performing SVD
		vec c = vec(trans(mean(S, 0)));
		double a1_mean = mean(a1);
		vec a2 = ones(S.n_rows);
		vec b2 = -(a1_mean*b1 + c);

		mat A = join_rows(a1, a2);
		mat B = join_rows(b1, b2);
		
		printf("\tPerform SVD on the original matrix\n"); fflush(stdout);
		vec s;
		mat U, V;
		field<mat> SVD_results;
		switch(SVD_algorithm) {
			case FULL_SVD:
				svd_econ(U, s, V, mat(S));			
				SVD_results(0) = U;
				SVD_results(1) = s;
				SVD_results(2) = V;
				break;			
			case IRLB_ALG:
				SVD_results = IRLB_SVD(S, dim, iter, seed);				
				break;
			case HALKO_ALG:
				SVD_results = HalkoSVD(S, dim, iter, seed);				
				break;
			case FENG_ALG:
				SVD_results = FengSVD(S, dim, iter, seed);				
				break;
			default:
				fprintf(stderr, "Unknown SVD algorithm chosen (%d). Switching to Halko.\n", SVD_algorithm);
				SVD_results = HalkoSVD(S, dim, iter, seed);				
				break;
		}

		U = SVD_results(0);
		s = SVD_results(1);
		V = SVD_results(2);
		
		printf("\tUpdate SVD ..."); fflush(stdout);
		vec s_prime;
		mat U_prime, V_prime;
			
		mat M = U.t()*A; 
		mat A_ortho_proj = A - U*M;   
		mat P = A_ortho_proj;// = orth(A_ortho_proj);
		gram_schmidt(P);		
		mat R_P = P.t()*A_ortho_proj;
		
		
		mat N = V.t()*B; 
		mat B_ortho_proj = B - V*N; 
		mat Q = B_ortho_proj; //orth(B_ortho_proj); 
		gram_schmidt(Q);
		mat R_Q = Q.t()*B_ortho_proj;	
		
		mat K1 = zeros(s.n_elem+A.n_cols, s.n_elem+A.n_cols);
		for(int i = 0; i < s.n_elem; i++) {
			K1(i, i) = s(i);
		}

		mat K2 = join_vert(M, R_P)*trans(join_vert(N, R_Q));

		
		mat K = K1 + K2;

		svd( U_prime, s_prime, V_prime, K );
		
		mat U_updated = join_horiz(U, P)*U_prime;
		mat V_updated = join_horiz(V, Q)*V_prime;
		printf("done.\n"); fflush(stdout);
		
		
		ReducedKernel output;		
		output.S_r = trans(V_updated.cols(0, dim-1));
		for(int i = 0; i < dim; i++) {
			output.S_r.row(i) *= s_prime(i);
		}
		output.V = U_updated.cols(0, dim-1);	
		
		vec sigma_sq = square(s_prime);

		output.lambda = sigma_sq / (n-1);

		vec explained_var = cumsum(sigma_sq) / sum(sigma_sq);
		output.exp_var = explained_var(span(0, dim-1));
		
		return output;	
	}

	ReducedKernel ACTION_reduction(mat &S, int dim, int iter = 5, int seed = 0, int SVD_algorithm = HALKO_ALG, bool prenormalize = false) {			
		int n = S.n_rows;
		if(prenormalize)
			S = normalise(S, 2);    

		printf("\tComputing reduced ACTION kernel. Input matrix size: %d x %d\n", S.n_rows, S.n_cols); fflush(stdout);
		
		// Update 1: Orthogonalize columns w.r.t. background (mean)
		vec mu = vec(mean(S, 1));
		vec v = mu / norm(mu, 2);
		vec a1 = v;
		vec b1 = -trans(S)*v;
		
		// Update 2: Center columns of orthogonalized matrix before performing SVD
		vec c = vec(trans(mean(S, 0)));
		double a1_mean = mean(a1);
		vec a2 = ones(S.n_rows);
		vec b2 = -(a1_mean*b1 + c);

		mat A = join_rows(a1, a2);
		mat B = join_rows(b1, b2);
		
		printf("\tPerform SVD on the original matrix\n"); fflush(stdout);
		vec s;
		mat U, V;
		field<mat> SVD_results;
		switch(SVD_algorithm) {
			case FULL_SVD:
				svd_econ(U, s, V, S);			
				SVD_results(0) = U;
				SVD_results(1) = s;
				SVD_results(2) = V;
				break;			
			case IRLB_ALG:
				SVD_results = IRLB_SVD(S, dim, iter, seed);				
				break;
			case HALKO_ALG:
				SVD_results = HalkoSVD(S, dim, iter, seed);				
				break;
			case FENG_ALG:
				SVD_results = FengSVD(S, dim, iter, seed);				
				break;
			default:
				fprintf(stderr, "Unknown SVD algorithm chosen (%d). Switching to Halko.\n", SVD_algorithm);
				SVD_results = HalkoSVD(S, dim, iter, seed);				
				break;
		}

		U = SVD_results(0);
		s = SVD_results(1);
		V = SVD_results(2);

		
		printf("\tUpdate SVD ..."); fflush(stdout);
		vec s_prime;
		mat U_prime, V_prime;
			
		mat M = U.t()*A; 
		mat A_ortho_proj = A - U*M;   
		mat P = A_ortho_proj;// = orth(A_ortho_proj);
		gram_schmidt(P);		
		mat R_P = P.t()*A_ortho_proj;
		
		
		mat N = V.t()*B; 
		mat B_ortho_proj = B - V*N; 
		mat Q = B_ortho_proj; //orth(B_ortho_proj); 
		gram_schmidt(Q);
		mat R_Q = Q.t()*B_ortho_proj;	
		
		mat K1 = zeros(s.n_elem+A.n_cols, s.n_elem+A.n_cols);
		for(int i = 0; i < s.n_elem; i++) {
			K1(i, i) = s(i);
		}

		mat K2 = join_vert(M, R_P)*trans(join_vert(N, R_Q));

		
		mat K = K1 + K2;

		svd( U_prime, s_prime, V_prime, K );
		
		mat U_updated = join_horiz(U, P)*U_prime;
		mat V_updated = join_horiz(V, Q)*V_prime;
		printf("done.\n"); fflush(stdout);
		
		ReducedKernel output;		
		output.S_r = trans(V_updated.cols(0, dim-1));
		for(int i = 0; i < dim; i++) {
			output.S_r.row(i) *= s_prime(i);
		}
		output.V = U_updated.cols(0, dim-1);	
		
		vec sigma_sq = square(s_prime);

		output.lambda = sigma_sq / (n-1);

		vec explained_var = cumsum(sigma_sq) / sum(sigma_sq);
		output.exp_var = explained_var(span(0, dim-1));
		
		return output;	
	}
	
	ReducedKernel reduce_kernel(sp_mat &S, int dim, int iter = 5, int seed = 0, int reduction_algorithm = ACTIONRED_ALG, int SVD_algorithm = HALKO_ALG, bool prenormalize = false) {
		ReducedKernel output = ACTION_reduction(S, dim, iter, seed, SVD_algorithm, prenormalize);
		
		return(output);
	}

	ReducedKernel reduce_kernel(mat &S, int dim, int iter = 5, int seed = 0, int reduction_algorithm = ACTIONRED_ALG, int SVD_algorithm = HALKO_ALG, bool prenormalize = false) {
		ReducedKernel output = ACTION_reduction(S, dim, iter, seed, SVD_algorithm, prenormalize);
		
		return(output);	
	}

	
	ReducedKernel SVD2ACTIONred(sp_mat &S, field<mat> SVD_results) {	
		int n = S.n_rows;
				
		printf("\tRetriving the original SVD results\n"); fflush(stdout);
		mat U = SVD_results(0);
		vec s = SVD_results(1);
		mat V = SVD_results(2);
		
		int dim = U.n_cols;
		
		// Update 1: Orthogonalize columns w.r.t. background (mean)
		vec mu = vec(mean(S, 1));
		vec v = mu / norm(mu, 2);
		vec a1 = v;
		vec b1 = -trans(S)*v;
		
		// Update 2: Center columns of orthogonalized matrix before performing SVD
		vec c = vec(trans(mean(S, 0)));
		double a1_mean = mean(a1);
		vec a2 = ones(S.n_rows);
		vec b2 = -(a1_mean*b1 + c);

		mat A = join_rows(a1, a2);
		mat B = join_rows(b1, b2);
		
	
		printf("\tUpdate SVD ..."); fflush(stdout);
		vec s_prime;
		mat U_prime, V_prime;
			
		mat M = U.t()*A; 
		mat A_ortho_proj = A - U*M;   
		mat P = A_ortho_proj;// = orth(A_ortho_proj);
		gram_schmidt(P);		
		mat R_P = P.t()*A_ortho_proj;
		
		
		mat N = V.t()*B; 
		mat B_ortho_proj = B - V*N; 
		mat Q = B_ortho_proj; //orth(B_ortho_proj); 
		gram_schmidt(Q);
		mat R_Q = Q.t()*B_ortho_proj;	
		
		mat K1 = zeros(s.n_elem+A.n_cols, s.n_elem+A.n_cols);
		for(int i = 0; i < s.n_elem; i++) {
			K1(i, i) = s(i);
		}

		mat K2 = join_vert(M, R_P)*trans(join_vert(N, R_Q));

		
		mat K = K1 + K2;

		svd( U_prime, s_prime, V_prime, K );
		
		mat U_updated = join_horiz(U, P)*U_prime;
		mat V_updated = join_horiz(V, Q)*V_prime;
		printf("done.\n"); fflush(stdout);
		
		
		ReducedKernel output;		
		output.S_r = trans(V_updated.cols(0, dim-1));
		for(int i = 0; i < dim; i++) {
			output.S_r.row(i) *= s_prime(i);
		}
		output.V = U_updated.cols(0, dim-1);	
		
		vec sigma_sq = square(s_prime);

		output.lambda = sigma_sq / (n-1);

		vec explained_var = cumsum(sigma_sq) / sum(sigma_sq);
		output.exp_var = explained_var(span(0, dim-1));
		
		return output;					
	}

	ReducedKernel SVD2ACTIONred(mat &S, field<mat> SVD_results) {			
		int n = S.n_rows;
		printf("\tComputing reduced ACTION kernel from prior SVD. Input matrix size: %d x %d\n", S.n_rows, S.n_cols); fflush(stdout);

				
		printf("\tRetriving the original SVD results\n"); fflush(stdout);
		mat U = SVD_results(0);
		vec s = SVD_results(1);
		mat V = SVD_results(2);
		printf("done.\n"); fflush(stdout);
				
		int dim = U.n_cols;
		
		// Update 1: Orthogonalize columns w.r.t. background (mean)
		vec mu = vec(mean(S, 1));
		vec v = mu / norm(mu, 2);
		vec a1 = v;
		vec b1 = -trans(S)*v;
		
		// Update 2: Center columns of orthogonalized matrix before performing SVD
		vec c = vec(trans(mean(S, 0)));
		double a1_mean = mean(a1);
		vec a2 = ones(S.n_rows);
		vec b2 = -(a1_mean*b1 + c);

		mat A = join_rows(a1, a2);
		mat B = join_rows(b1, b2);
		
	
		printf("\tUpdate SVD ..."); fflush(stdout);
		vec s_prime;
		mat U_prime, V_prime;
			
		mat M = U.t()*A; 
		mat A_ortho_proj = A - U*M;   
		mat P = A_ortho_proj;// = orth(A_ortho_proj);
		gram_schmidt(P);		
		mat R_P = P.t()*A_ortho_proj;
		
		
		mat N = V.t()*B; 
		mat B_ortho_proj = B - V*N; 
		mat Q = B_ortho_proj; //orth(B_ortho_proj); 
		gram_schmidt(Q);
		mat R_Q = Q.t()*B_ortho_proj;	
		
		mat K1 = zeros(s.n_elem+A.n_cols, s.n_elem+A.n_cols);
		for(int i = 0; i < s.n_elem; i++) {
			K1(i, i) = s(i);
		}

		mat K2 = join_vert(M, R_P)*trans(join_vert(N, R_Q));

		
		mat K = K1 + K2;

		svd( U_prime, s_prime, V_prime, K );
		
		mat U_updated = join_horiz(U, P)*U_prime;
		mat V_updated = join_horiz(V, Q)*V_prime;
		printf("done.\n"); fflush(stdout);
		
		
		ReducedKernel output;		
		output.S_r = trans(V_updated.cols(0, dim-1));
		for(int i = 0; i < dim; i++) {
			output.S_r.row(i) *= s_prime(i);
		}
		output.V = U_updated.cols(0, dim-1);	
		
		vec sigma_sq = square(s_prime);

		output.lambda = sigma_sq / (n-1);

		vec explained_var = cumsum(sigma_sq) / sum(sigma_sq);
		output.exp_var = explained_var(span(0, dim-1));
		
		return output;	
	}


	ReducedKernel PCA2ACTIONred(sp_mat &S, field<mat> PCA_results) {	
				
		printf("\tReverting column-centering ... "); fflush(stdout);
		field<mat> SVD_results = PCA2SVD(S, PCA_results);
		printf("done reverting\n");
		
		ReducedKernel output = SVD2ACTIONred(S, SVD_results);
		return output;					
	}
	

	ReducedKernel PCA2ACTIONred(mat &S, field<mat> PCA_results) {	
				
		printf("\tReverting column-centering ... "); fflush(stdout);
		field<mat> SVD_results = PCA2SVD(S, PCA_results);
		printf("done reverting\n");
		
		ReducedKernel output = SVD2ACTIONred(S, SVD_results);
		return output;					
	}	



}
