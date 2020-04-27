#include <ACTIONet.h>

namespace ACTIONet {
	field<mat> compute_feature_specificity(sp_mat &S, mat &H) {			
		field<mat> res(3);

		sp_mat Sb = spones(S);
		
		printf("Compute row stats");
		vec row_means = vec(mean(S, 1));
		vec row_p = vec(mean(Sb, 1));
		vec alpha = row_means / row_p;
		printf("done\n");
						
		
		printf("Computing column stats ... ");
		rowvec col_p = rowvec(mean(Sb, 0));
		double rho = mean(col_p);
		rowvec beta = col_p/rho; // Relative density compared to the overall density
		
		// Mix the two weights together
		mat Gamma = H;
		vec a(H.n_rows);
		for(int i = 0; i < H.n_rows; i++) {
			Gamma.row(i) %= beta;
			a(i) = max(Gamma.row(i));
		}		
		printf("done\n");
		
		
		printf("Computing observation statistics ... ");
		mat Obs = S * trans(H);
		printf("done\n");
		
		
    
		printf("Computing expectation statistics ... ");
		mat Exp = row_means * trans(sum(Gamma, 1));
		mat Nu  = row_means * trans(sum(square(Gamma), 1));
		printf("done\n");
		
		
		printf("Computing significance ... ");
		
		mat Lambda = Obs - Exp;

		
		mat logPvals_lower = square(Lambda)/(2 * Nu);

		uvec uidx = find(Lambda>=0);
		logPvals_lower(uidx) = zeros(uidx.n_elem);
		
		mat Lambda_scaled = Lambda;
		for(int j = 0; j < Lambda_scaled.n_cols; j++) {
			Lambda_scaled.col(j) *= (a(j) / 3);
		}
		mat logPvals_upper = square(Lambda)/(2*(Nu + Lambda_scaled));
		
		uvec lidx = find(Lambda<=0);
		logPvals_upper(lidx) = zeros(lidx.n_elem);

		logPvals_lower /= log(10);
		logPvals_upper /= log(10);
    

		res(0) = Obs;
		res(1) = logPvals_upper;
		res(2) = logPvals_lower;
		
		return(res);
	}	
	
	field<mat> compute_feature_specificity(mat &S, mat &H) {			
		field<mat> res(3);
		
		mat Sb = S;		
		uvec nnz_idx = find(Sb > 0);
		(Sb(nnz_idx)).ones();
		
		printf("Compute row stats");
		vec row_means = vec(mean(S, 1));
		vec row_p = vec(mean(Sb, 1));
		vec alpha = row_means / row_p;
		printf("done\n");
						
		
		printf("Computing column stats ... ");
		rowvec col_p = mean(Sb, 0);
		double rho = mean(col_p);
		rowvec beta = col_p/rho; // Relative density compared to the overall density
		
		// Mix the two weights together
		mat Gamma = H;
		vec a(H.n_rows);
		for(int i = 0; i < H.n_rows; i++) {
			Gamma.row(i) %= beta;
			a(i) = max(Gamma.row(i));
		}		
		printf("done\n");
		
		
		printf("Computing observation statistics ... ");
		mat Obs = S * trans(H);
		printf("done\n");
		
		
    
		printf("Computing expectation statistics ... ");
		mat Exp = row_means * trans(sum(Gamma, 1));
		mat Nu  = row_means * trans(sum(square(Gamma), 1));
		printf("done\n");
		
		
		printf("Computing significance ... ");
		
		mat Lambda = Obs - Exp;

		
		mat logPvals_lower = square(Lambda)/(2 * Nu);
		uvec uidx = find(Lambda>=0);
		logPvals_lower(uidx) = zeros(uidx.n_elem);
		
				
		mat Lambda_scaled = Lambda;
		for(int j = 0; j < Lambda_scaled.n_cols; j++) {
			Lambda_scaled.col(j) *= (a(j) / 3);
		}
		mat logPvals_upper = square(Lambda)/(2*(Nu + Lambda_scaled));
		uvec lidx = find(Lambda<=0);
		logPvals_upper(lidx) = zeros(lidx.n_elem);
		
		
		logPvals_lower /= log(10);
		logPvals_upper /= log(10);
    

		res(0) = Obs;
		res(1) = logPvals_upper;
		res(2) = logPvals_lower;

		
		return(res);
	}	
	
	
	field<mat> compute_feature_specificity(sp_mat &S, uvec sample_assignments) {			
		mat H(max(sample_assignments), S.n_cols);
		
		for(int i = 1; i <= max(sample_assignments); i++) {
			vec v = zeros(S.n_cols);
			uvec idx = find(sample_assignments == i);
			v(idx) = ones(idx.n_elem);
			H.row(i-1) = trans(v);
		}
		
		field<mat> res = compute_feature_specificity(S, H);
		
		return(res);		
	}	
	
	
	field<mat> compute_feature_specificity(mat &S, uvec sample_assignments) {			
		mat H(max(sample_assignments), S.n_cols);
		
		for(int i = 1; i <= max(sample_assignments); i++) {
			vec v = zeros(S.n_cols);
			uvec idx = find(sample_assignments == i);
			v(idx) = ones(idx.n_elem);
			H.row(i-1) = trans(v);
		}
		
		field<mat> res = compute_feature_specificity(S, H);
		
		return(res);		
	}		
	
}
