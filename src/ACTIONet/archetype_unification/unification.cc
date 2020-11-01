#include "ACTIONet.h"
#include "dagConstruct.h"

vector<double> Corrector::vals;


namespace ACTIONet {
	mat unsigned_cluster_batch(sp_mat A, vec resolutions, uvec initial_clusters = uvec(), int seed = 0);
	
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

	unification_results unify_archetypes(mat &S_r, mat &C_stacked, mat &H_stacked, double sensitivity = 1.0, int normalization_type = 1, double edge_threshold = 0.5) {
		printf("Unify archetypes: sensitivity = %.2f (%d archs)\n", sensitivity, H_stacked.n_rows);
								

		unification_results output;

		mat C_norm = normalise(C_stacked, 1, 0);
		//mat W_stacked = S_r * C_norm;
		mat arch_H_stacked = H_stacked * C_norm;
		mat arch_S_r = S_r * C_norm;
		
		mat Wg = trans(arch_S_r) * arch_S_r;		
		Wg = Wg / max(max(abs(Wg)));
		mat Hg = computeFullSim(arch_H_stacked, 1);
		mat G = Wg % Hg;
		
		G.transform( [edge_threshold](double val) { return (val < edge_threshold? 0:val); } );

		if(normalization_type == 1) {
			G = normalise(G, 1, 0);
		}
		else if(normalization_type == 3) {
			G = NetEnh(G);						
		}
		
		int dim = std::min((int)G.n_cols, 100);
		SPA_results res = run_SPA(G, dim);
		
		uvec selected_columns = res.selected_columns;

	
		printf("Computing lead independent vertex set\n"); fflush(stdout);	
		mat subG = Hg(selected_columns, selected_columns);
		subG.diag().zeros();
		subG = normalise(subG, 1, 0);
			
		vec is_selected = zeros(selected_columns.n_elem);
		for(int i = 0; i < selected_columns.n_elem; i++) {
			vec N = subG.col(i);
			vec neighbors_contrib = vec( trans(N) * is_selected );
			double w = neighbors_contrib(0);
			if(w <= sensitivity) {
				//printf("pick %d (N. contrib = %f)\n", i, w);
				is_selected(i) = 1;
			} else {
				//printf("drop %d (N. contri = %f)\n", i, w);
			}
		}
		printf("Selected %d states\n", (int)sum(is_selected));

		uvec idx = find(is_selected == 1);
		uvec selected_archs = selected_columns(idx);

		mat G_ext = subG;
		G_ext.diag() = 5.0 * ones(G_ext.n_rows);
		G_ext = normalise(G_ext, 1, 0);
		
		mat C_avg = C_norm.cols(selected_columns) * G_ext;
		mat C_unified = C_avg.cols(idx);


		//mat C_unified = C_norm.cols(selected_columns) * subG.cols(idx);
		
		/*
		mat G_ext = G;
		G_ext.diag() = 5.0 * ones(G_ext.n_rows);
		G_ext = normalise(G_ext, 1, 0);
		
		mat C_avg = C_norm * G_ext;
		mat C_unified = C_avg.cols(selected_archs);
		*/
		
		//mat C_unified = C_norm.cols(selected_archs);


/*
		field<vec> res = run_HDBSCAN(G, 10, 5);
		
		vec cl = res(0);
		vec cl_p = res(1);
		
		int d = arma::max(cl);
		mat M = zeros(d, cl.n_elem);		
		for(int ii = 0; ii < cl.n_elem; ii++) {
			if(cl(ii) > 0) {
				M(cl(ii)-1, ii) = cl_p(ii);
			}
		}
		mat C_unified = C_norm * trans(M);

		

		/*
		mat A = H_stacked * C_unified;
		mat H_unified = run_simplex_regression(A, H_stacked, false);			
		*/

/*
		sp_mat A = sp_mat(G);
		uvec initial_clusters_uvec(A.n_rows);
		for (int i = 0; i < A.n_rows; i++) initial_clusters_uvec(i) = i;
		
		vec cl = unsigned_cluster(A, sensitivity, initial_clusters_uvec, 0);
		
		int d = arma::max(cl);
		mat M = zeros(d+1, cl.n_elem);		
		for(int ii = 0; ii < cl.n_elem; ii++) {
				M(cl(ii), ii) = 1;
		}
		vec cc = sum(M, 1);
		cc.print("cc");
		
		uvec idx = find(cc >= 5);
		
		mat C_unified = C_norm * trans(M.rows(idx));
*/
		
		C_unified = normalise(C_unified, 1, 0);					
		
		mat W_unified = S_r * C_unified;
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
