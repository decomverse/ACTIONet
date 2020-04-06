#include "ACTIONet.h"


namespace ACTIONet {

	unification_results unify_archetypes(sp_mat &G, mat &S_r, mat &C_stacked, mat &H_stacked) {
		
		mat W_r = S_r * C_stacked;
		W_r = normalise(W_r, 1); 
		
		int k = std::min(W_r.n_rows, W_r.n_cols);		
		SPA_results SPA_out = run_SPA(W_r, k);
		
		double x_sum = sum(SPA_out.column_norms);
		double x2_sum = sum(square(SPA_out.column_norms));
		int nnz = round( (x_sum*x_sum) / x2_sum );
		
		mat subW_r = W_r.cols(SPA_out.selected_columns(span(0, nnz-1)));
		mat arch_membership = run_simplex_regression(subW_r, W_r);		
		uvec archetype_groups = trans(index_max( arch_membership, 0 ));
		
		mat subC;
		mat C_unified(C_stacked.n_rows, nnz);
		for(int i = 0; i < nnz; i++) {
			uvec idx = find(archetype_groups == i);
			if(idx.n_elem == 1) {
				subC = C_stacked.col(idx(0));
				C_unified.col(i) = subC;
			} else {				
				subC = C_stacked.cols(idx);
				C_unified.col(i) = mean(subC, 1);				
			}
		}
		mat W_unified = S_r * C_unified;
		mat H_unified = run_simplex_regression(W_unified, S_r);
		
		uvec sample_assignments = trans(index_max( H_unified, 0 ));
		
		unification_results output;
		
		output.archetype_groups = archetype_groups;
		output.C_unified = C_unified;
		output.H_unified = H_unified;
		output.sample_assignments = sample_assignments;

		return(output);
	}
}
