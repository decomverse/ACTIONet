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

	unification_results unify_archetypes(sp_mat& G,
										mat &S_r,
										mat &C_stacked,
										double alpha = 0.85,
										int core_threshold = 1,
										double sim_threshold = 0.5,
										int thread_no = 0) {
		printf("Unify archetypes (%d archs, alpha = %f, core_threshold = %d, sim_threshold = %f)\n", C_stacked.n_cols, alpha, core_threshold, sim_threshold);
														
		unification_results output;
		
		// Smooth archetypes using ACTIONet
		sp_mat X0 = sp_mat(normalise(C_stacked, 1, 0));
		mat C_filtered = compute_network_diffusion(G, X0, thread_no, alpha, 5);
		
		
		
		// Build archetype-archetype network from smoothed archetype footprints
		double M = 16, ef_construction = 200, ef = 50, density = 1.0;
		sp_mat arch_G = build_ACTIONet(C_filtered, density, thread_no, M, ef_construction, ef, true);


		// Prune unreliable archetypes
		uvec core_num = compute_core_number(arch_G);
		uvec selected_archetypes = find(core_threshold < core_num);
		
		arch_G = sp_mat(mat(arch_G)(selected_archetypes, selected_archetypes));


		// Transform similarity matrix
		vec resolutions = regspace(0.1, 0.1, 5);
		mat X = trans(unsigned_cluster_batch(arch_G, resolutions));			
		mat Sim = trans(X) * X;
		Sim = Sim / max(max(Sim));

/*		
		output.C_unified = mat(Sim);
		output.H_unified = mat(arch_G);
		output.selected_archetypes = selected_archetypes;
		return(output);
*/		
		// Reduced profile of each archetype
		C_filtered = C_filtered.cols(selected_archetypes);
		mat S_r_arch = S_r * C_filtered;
		
		
		// Prioritize archetypes
		int dim = min(100, (int)min(X.n_cols, X.n_rows));		
		SPA_results res = run_SPA(X, dim);		
		uvec selected_columns = res.selected_columns;		
		uvec row_idx = find(1e-10 < res.column_norms);		
		selected_columns = selected_columns(row_idx);
		
		dim = selected_columns.n_elem;	

						
		vec is_selected = zeros(S_r_arch.n_cols);
		is_selected(selected_columns(0)) = 1;		
		for(int i = 1; i < dim; i++) {			
			int j = selected_columns[i];
			vec v = Sim.col(j);
			vec cc_vals = v(find(is_selected == 1));
			double mm = max(cc_vals);
			
			if(sim_threshold < mm) {
				continue;
			}
			is_selected(j) = 1;
		}		
		uvec idx = find(is_selected == 1);		
		
		mat C_unified = C_filtered.cols(idx);
		mat W_unified = S_r_arch.cols(idx);
				
				
		mat H_unified = run_simplex_regression(W_unified, S_r, false);			


		uvec assigned_archetypes = trans(index_max( H_unified, 0 ));
		
		output.C_unified = C_unified;
		output.H_unified = H_unified;
		output.assigned_archetypes = assigned_archetypes;
		output.selected_archetypes = selected_archetypes(idx);		
		
		return(output);
	}
	
}
