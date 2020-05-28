#include "ACTIONet.h"


namespace ACTIONet {
	
	double Kappa(double p, double q) {
		double a = 0.0, b = 0.0;
		if( (1e-300 < p) & (1e-300 < q) ){
			a = p * log(p/q);
		}
		if( (p < (1-1e-300)) & (q < (1-1e-300)) ) {
			b = (1-p)*log((1-p)/(1-q));
		} 		
		
		double k = a + b;
		return(k);			
	}
		
	double log_HGT_tail(int population_size, int success_count, int sample_size, int observed_success) {
		if ( observed_success == 0) 
			return(0);
			
		double success_rate = (double)success_count/population_size;		
		double expected_success = sample_size * success_rate;
		double delta = (observed_success/expected_success) - 1.0;
		if(delta < 0) {
			return(0);
		}
		
		double log_tail_bound = sample_size * Kappa((1.0 + delta) * success_rate, success_rate);
		
		return(log_tail_bound);
	}
	
	double assess_overlap(uvec i1, uvec i2, int population_size) {
		int success_count = i1.n_elem;
		int sample_size = i2.n_elem;
		
		uvec shared = intersect(i1, i2);
		int observed_success = shared.n_elem;
		
		double log_pval = log_HGT_tail(population_size, success_count, sample_size, observed_success);
		
		return(log_pval);		
	}
	
	mat compute_overlap_matrix(mat C) {
		int N = C.n_cols;
		mat O = zeros(N, N);
		
		vector<uvec> indices(N);
		for(int i = 0; i < N; i++) {
			uvec idx = find(C.col(i) > 0);
			indices[i] = idx;
		}
		
		for(int i = 0; i < N; i++) {
			uvec i1 = indices[i];
			for(int j = i+1; j < N; j++) {
				uvec i2 = indices[j];
				
				O(i, j) = O(j, i) = assess_overlap(i1, i2, N);				
			}
		}
		
		return(O);
	}
	
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

	//unification_results unify_archetypes(sp_mat &G, mat &S_r, mat &archetypes, mat &C_stacked, mat &H_stacked, int minPoints = 10, int minClusterSize = 10, double outlier_threshold = 0.0, int reduced_dim = 50) {
	unification_results unify_archetypes(mat &S_r, mat &C_stacked, mat &H_stacked, double min_overlap = 10.0, double resolution = 1.0) {
		unification_results output;

		mat A_full = compute_overlap_matrix(C_stacked);
		
		
		uvec idx = find(A_full < min_overlap);
		if(0 < idx.n_elem)
			A_full(idx).zeros();
			
		
		sp_mat A(A_full);
		
		
		uvec initial_clusters(A.n_rows);
		for (int i = 0; i < A.n_rows; i++) initial_clusters(i) = i;

		vec clusters = unsigned_cluster(A, resolution, initial_clusters, 0);
		
		uvec cn = compute_core_number(A);
		
		
		uvec selected_archs = find(1 < cn);

		output.archetype_groups = clusters;
		output.selected_archetypes = selected_archs;		
		
		clusters = clusters(selected_archs);

	
		/*
		ReducedKernel reduction = reduce_kernel(archetypes, reduced_dim, 5, 0, 1, 1);

		mat X = trans(reduction.S_r);
	
		field<vec> res = run_HDBSCAN(X, minPoints, minClusterSize);

		vec clusters = res(0);
		uvec selected_archs = find((0 < res(0)) && (res(2) <= outlier_threshold));
		clusters = clusters(selected_archs);

		output.archetype_groups = res(0);
		uvec filtered_archs = find((0 >= res(0)) && (res(2) > outlier_threshold));
		output.archetype_groups(filtered_archs).zeros();
		*/
		

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
		
		uvec assigned_archetypes = trans(index_max( H_unified, 0 ));
		
		output.C_unified = C_unified;
		output.H_unified = H_unified;
		output.assigned_archetypes = assigned_archetypes;
		

		return(output);
	}
}
