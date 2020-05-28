#include <ACTIONet.h>

namespace ACTIONet {
	field<mat> compute_feature_specificity_bin(sp_mat &Sb, mat &H) {			
		field<mat> res(3);

		mat Ht = trans(H);
		Ht.each_col( [](vec& h){ double mu = mean(h); h /= (mu == 0)?1:mu; } ); // For numerical stability


		printf("Compute stats ... ");
		//vec row_p = vec(mean(Sb, 1));
		//rowvec col_p = rowvec(mean(Sb, 0));
		
		// Heuristic optimization! Shall add parallel for later on
		vec row_p = zeros(Sb.n_rows);		
		vec col_p = zeros(Sb.n_cols);				
		sp_mat::const_iterator it     = Sb.begin();
		sp_mat::const_iterator it_end = Sb.end();
		for(; it != it_end; ++it) {			
			row_p[it.row()]++;
			col_p[it.col()]++;			
		}		
		row_p /= Sb.n_cols;
		col_p /= Sb.n_rows;		
		printf("done\n");
			
				
		printf("Computing observation statistics ... ");
		sp_mat Ht_sparse = sp_mat(Ht); // Shall we sparsify first?
		mat Obs = mat(Sb * Ht_sparse);
		printf("done\n");		
		
    
		printf("Computing expectation statistics ... ");
		double rho = mean(col_p);
		vec beta = col_p/rho; // Relative density compared to the overall density
		mat Gamma = Ht;
		vec a(H.n_rows);
		for(int i = 0; i < H.n_rows; i++) {
			Gamma.col(i) %= beta;
			a(i) = max(Gamma.col(i));
		}		
		
		mat Exp = row_p * sum(Gamma, 0);
		mat Nu  = row_p * sum(square(Gamma), 0);
		printf("done\n");
		

		
		printf("Computing significance ... ");		
		mat Lambda = Obs - Exp;

		
		mat logPvals_lower = square(Lambda)/(2 * Nu);
		uvec uidx = find(Lambda>=0);
		logPvals_lower(uidx) = zeros(uidx.n_elem);
		logPvals_lower.replace(datum::nan, 0);  // replace each NaN with 0
		
				
		mat Lambda_scaled = Lambda;
		for(int j = 0; j < Lambda_scaled.n_cols; j++) {
			Lambda_scaled.col(j) *= (a(j) / 3);
		}
		mat logPvals_upper = square(Lambda)/(2*(Nu + Lambda_scaled));
		uvec lidx = find(Lambda<=0);
		logPvals_upper(lidx) = zeros(lidx.n_elem);
		logPvals_upper.replace(datum::nan, 0);  // replace each NaN with 0
		
		
		logPvals_lower /= log(10);
		logPvals_upper /= log(10);
		printf("done\n");
    

		res(0) = Obs / Ht.n_rows;
		res(1) = logPvals_upper;
		res(2) = logPvals_lower;
		
		return(res);
	}	



	field<mat> compute_feature_specificity(sp_mat &S, mat &H) {			
		field<mat> res(3);

		mat Ht = trans(H);
		Ht.each_col( [](vec& h){ double mu = mean(h); h /= (mu == 0)?1:mu; } ); // For numerical stability
		
		// make sure all values are positive
		double min_val = S.min();
		S.for_each( [min_val](mat::elem_type& val) { val -= min_val; } );

		printf("Compute stats ... ");
		//sp_mat Sb = spones(S);		
		
		// Heuristic optimization! Shall add parallel for later on
		vec row_p = zeros(S.n_rows);		
		vec col_p = zeros(S.n_cols);
		vec row_factor = zeros(S.n_rows);
						
		sp_mat::const_iterator it     = S.begin();
		sp_mat::const_iterator it_end = S.end();
		for(; it != it_end; ++it) {			
			col_p[it.col()]++;			
			row_p[it.row()]++;
			row_factor[it.row()] += (*it);
		}		
		row_factor /= row_p;		
		row_p /= S.n_cols;
		col_p /= S.n_rows;		
		
		
		printf("done\n");
			
				
		printf("Computing observation statistics ... ");
		sp_mat Ht_sparse = sp_mat(Ht); // Shall we sparsify first?
		mat Obs = mat(S * Ht_sparse);
		printf("done\n");		
		
    
		printf("Computing expectation statistics ... ");
		double rho = mean(col_p);
		vec beta = col_p/rho; // Relative density compared to the overall density
		mat Gamma = Ht;
		vec a(H.n_rows);
		for(int i = 0; i < H.n_rows; i++) {
			Gamma.col(i) %= beta;
			a(i) = max(Gamma.col(i));
		}		
		
		mat Exp = (row_p % row_factor) * sum(Gamma, 0);
		mat Nu  = (row_p % square(row_factor)) * sum(square(Gamma), 0);		
		mat A = (row_factor * trans(a));
		printf("done\n");
		

		
		printf("Computing significance ... ");		
		mat Lambda = Obs - Exp;

		
		mat logPvals_lower = square(Lambda)/(2 * Nu);
		uvec uidx = find(Lambda>=0);
		logPvals_lower(uidx) = zeros(uidx.n_elem);
		logPvals_lower.replace(datum::nan, 0);  // replace each NaN with 0


		mat logPvals_upper = square(Lambda)/(2*(Nu + (Lambda%A/3)));
		uvec lidx = find(Lambda<=0);
		logPvals_upper(lidx) = zeros(lidx.n_elem);
		logPvals_lower.replace(datum::nan, 0);  // replace each NaN with 0
		
		
		logPvals_lower /= log(10);
		logPvals_upper /= log(10);
		printf("done\n");
    

		res(0) = Obs / Ht.n_rows;
		res(1) = logPvals_upper;
		res(2) = logPvals_lower;
		
		return(res);
	}	


	field<mat> compute_feature_specificity(mat &S, mat &H) {			
		field<mat> res(3);

		// make sure all values are positive
		double min_val = S.min();
		S.for_each( [min_val](mat::elem_type& val) { val -= min_val; } );

		mat Sb = S;		
		uvec nnz_idx = find(Sb > 0);
		(Sb(nnz_idx)).ones();

		mat Ht = trans(H);
		Ht.each_col( [](vec& h){ double mu = mean(h); h /= (mu == 0)?1:mu; } ); // For numerical stability

		printf("Compute stats ... ");
		vec row_p = vec(sum(Sb, 1));
		vec row_factor = vec(sum(S, 1)) / row_p; // mean of nonzero elements
		row_p /= Sb.n_cols;		
		vec col_p = vec(trans(mean(Sb, 0)));
		printf("done\n");
			
				
		printf("Computing observation statistics ... ");
		mat Obs = mat(S * Ht);
		printf("done\n");		
		
    
		printf("Computing expectation statistics ... ");
		double rho = mean(col_p);
		vec beta = col_p/rho; // Relative density compared to the overall density
		mat Gamma = Ht;
		vec a(H.n_rows);
		for(int i = 0; i < H.n_rows; i++) {
			Gamma.col(i) %= beta;
			a(i) = max(Gamma.col(i));
		}		
		
		mat Exp = (row_p % row_factor) * sum(Gamma, 0);
		mat Nu  = (row_p % square(row_factor)) * sum(square(Gamma), 0);		
		mat A = (row_factor * trans(a));
		printf("done\n");
		

		
		printf("Computing significance ... ");		
		mat Lambda = Obs - Exp;

		
		mat logPvals_lower = square(Lambda)/(2 * Nu);
		uvec uidx = find(Lambda>=0);
		logPvals_lower(uidx) = zeros(uidx.n_elem);
		logPvals_lower.replace(datum::nan, 0);  // replace each NaN with 0


		mat logPvals_upper = square(Lambda)/(2*(Nu + (Lambda%A/3)));
		uvec lidx = find(Lambda<=0);
		logPvals_upper(lidx) = zeros(lidx.n_elem);
		logPvals_lower.replace(datum::nan, 0);  // replace each NaN with 0
		
		
		
		logPvals_lower /= log(10);
		logPvals_upper /= log(10);
		printf("done\n");
    

		res(0) = Obs / Ht.n_rows;
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


