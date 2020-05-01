#include "ACTIONet.h"


namespace ACTIONet {
	
	field<vec> run_HDBSCAN(mat &X, int minPoints = 5, int minClusterSize = 5) {
		Hdbscan hdbscan(X);
		hdbscan.execute(minPoints, minClusterSize, "Euclidean");
		
		vec labels(X.n_rows);
		vec membershipProbabilities(X.n_rows);
		vec outlierScores(X.n_rows);
		
		for(int i = 0; i < X.n_rows; i++) {
			labels[i] = hdbscan.labels_[i];
			membershipProbabilities[i] = hdbscan.membershipProbabilities_[i];
			outlierScores[hdbscan.outlierScores_[i].id] = hdbscan.outlierScores_[i].score;
		}
		
		field<vec> out(3);
		out(0) = labels;
		out(1) = membershipProbabilities;
		out(2) = outlierScores;
		
		return(out);
	}

	unification_results unify_archetypes(sp_mat &G, mat &S_r, mat &C_stacked, mat &H_stacked, int minPoints = 5, int minClusterSize = 5, double outlier_threshold = 0.0) {
		unification_results output;
		

		mat W_r = S_r * C_stacked;
		mat X = trans(W_r);

	
		field<vec> res = run_HDBSCAN(X, minPoints, minClusterSize);

		vec clusters = res(0);
		uvec selected_archs = find((0 < res(0)) && (res(2) <= outlier_threshold));
		clusters = clusters(selected_archs);

		output.archetype_groups = res(0);
		output.selected_archetypes = selected_archs;		

		vec uc = sort(unique(clusters));
		
		mat subC;
		mat C_unified(C_stacked.n_rows, uc.n_elem);
		for(int i = 0; i < uc.n_elem; i++) {
			
			uvec ii = find(clusters == uc(i));
			uvec idx = selected_archs(ii);
			
			subC = C_stacked.cols(idx);
			if(idx.n_elem == 1) {
				C_unified.col(i) = subC;
			} else {				
				C_unified.col(i) = mean(subC, 1);				
			}
		}
		C_unified = normalise(C_unified, 1, 0);
		
		mat W_unified = S_r * C_unified;
		mat H_unified = run_simplex_regression(W_unified, S_r);
		
		uvec sample_assignments = trans(index_max( H_unified, 0 ));
		
		output.C_unified = C_unified;
		output.H_unified = H_unified;
		output.sample_assignments = sample_assignments;

		return(output);
	}
}
