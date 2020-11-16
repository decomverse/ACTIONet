#include "ACTIONet.h"
#include "dagConstruct.h"

vector<double> Corrector::vals;


namespace ACTIONet {
	mat unsigned_cluster_batch(sp_mat A, vec resolutions, uvec initial_clusters = uvec(), int seed = 0);
	mat signed_cluster_batch(sp_mat A, vec resolutions, uvec initial_clusters = uvec(), int seed = 0);
	
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
	
	
	mat NetEnh(mat A) {		
		A.diag().zeros();
		mat P = normalise(A, 1, 1);
		
		mat D = diagmat(sqrt(sum(P) + 1e-16));
		mat W = P * D;		
		mat P2 = W * trans(W);		
		//P2.diag().zeros();
		
		return(P2);
	}

	unification_results unify_archetypes(mat &S_r, mat &C_stacked, mat &H_stacked, double z_threshold = -1.0, double cor_threshold = 0.95) {
		printf("Unify archetypes (%d archs)\n", H_stacked.n_rows);
														
		unification_results output;

		mat C_norm = normalise(C_stacked, 1, 0);
		mat Ht_norm = normalise(trans(H_stacked), 1, 0);					
		mat S_r_arch = S_r * Ht_norm;
		
		
		// Filter noisy archetypes		
		vec x = trans(sqrt(sum(square(S_r_arch))));
		vec z = robust_zscore(x);
		S_r_arch = S_r_arch.cols(find(z > z_threshold));
		C_norm = C_norm.cols(find(z > z_threshold));
		
		
		// Normalize wrt norm-1
		S_r_arch = normalise(S_r_arch, 1, 0);
		
				
		// Run SPA to order columns
		int dim = min(100, (int)min(S_r_arch.n_cols, S_r_arch.n_rows));		
		SPA_results res = run_SPA(S_r_arch, dim);
		uvec selected_columns = res.selected_columns;
		
		
		
		// Prune selected columns in order
		mat CC = cor(S_r_arch);
		vec is_selected = zeros(S_r_arch.n_cols);
		is_selected(selected_columns(0)) = 1;		
		for(int i = 1; i < dim; i++) {			
			int j = selected_columns[i];
			vec v = CC.col(j);
			vec cc_vals = v(find(is_selected == 1));
			double mm = max(cc_vals);
			
			if(cor_threshold < mm) {
				continue;
			}
			is_selected(j) = 1;
		}		
		uvec idx = find(is_selected == 1);		
		
		mat C_unified = C_norm.cols(idx);
		mat W_unified = S_r_arch.cols(idx);
				
		
		field<mat> res_AA = run_AA(S_r_arch, W_unified, 50, 1e-16);
		mat C_AA = res_AA(0);	
		mat H_AA = res_AA(1);		
						
		W_unified = S_r_arch * C_AA;
		C_unified = C_norm * C_AA;
		
		
		mat H_unified = run_simplex_regression(W_unified, S_r, false);			


		uvec assigned_archetypes = trans(index_max( H_unified, 0 ));
		
		output.C_unified = C_unified;
		output.H_unified = H_unified;
		output.assigned_archetypes = assigned_archetypes;
		
		//output.dag_adj = dag_adj;
		//output.dag_node_annotations = dag_node_annotations;

		return(output);
	}
	
}
