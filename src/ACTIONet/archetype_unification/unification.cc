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

	unification_results unify_archetypes(mat &S_r, mat &C_stacked, mat &H_stacked, double min_edge_weight = 0.5, int min_coreness = 2, double resolution = 1.0, int min_repeat = 2, int thread_no = 0, double alpha = 0.05, double beta = 0.5, double outlier_threshold = 0.5, int minPoints = 5, int minClusterSize = 5, double cond_threshold = 10, int normalization_type = 0, bool preprocess_adj = true, bool reduce_G = false, int method_type = 3) {
		printf("Unify archetypes: resolution = %.2f, min. weight = %.2f, alpha = %.2f, beta = %.2f\n", resolution, min_edge_weight, alpha, beta);
		
		
		
		unification_results output;

		mat C_norm = normalise(C_stacked, 1, 0);
		//mat W_stacked = S_r * C_norm;
		mat arch_H_stacked = H_stacked * C_norm;
		
		
		mat G = computeFullSim(arch_H_stacked, thread_no);
		if(preprocess_adj == true) {
			vec resolutions = regspace(0.1, 0.1, 10);
			mat M_total = unsigned_cluster_batch(sp_mat(G), resolutions);
			G = M_total * trans(M_total);
		}
		
		if(normalization_type == 1) {
			G = NetEnh(G);
		} else if(normalization_type == 2) {
			vec cs = sum(G, 1);
			for(int i = 0; i < G.n_rows; i++) {
				for(int j = i+1; j < G.n_cols; j++) {
					double w = sqrt(cs(i)*cs(j));
					G(i, j) /= (w == 0?1:w);
				}
			}		
		} else if(normalization_type == 3) {
			G = normalise(G, 1, 0);					
		}

		int dim = G.n_cols;
		mat G_red;
		if(reduce_G == true) {
			int dim = std::min((int)G.n_cols, 50);
			
			field<mat> reduction = reduce_kernel(G, dim, 500, 0, IRLB_ALG, false);
					
			vec sigma = reduction(1).col(0);
			G_red = reduction(2);
			for(int i = 0; i < G_red.n_cols; i++) {
				G_red.col(i) *= sigma(i);
			}	
			G_red = trans(G_red);
		} else {
			G_red = (G);
		}
		printf("GRed: %d x %d\n", G_red.n_rows, G_red.n_cols);
		
		
		mat C_unified;
		if(method_type == 1) {
			printf("Running SPA (G:%d x %d, sim = %d)\n", G_red.n_rows, G_red.n_cols, dim); fflush(stdout);
			SPA_results res = run_SPA(G_red, dim);
			
			uvec selected_columns = res.selected_columns;
			
			printf("Computing condition numbers\n"); fflush(stdout);
			vec cn = zeros(dim-1);
			for(int i = 1; i < dim; i++) {
				cn(i-1) = cond(G_red.cols(selected_columns(span(0, i))));
				//cn(i-1) = cond(W_stacked.cols(selected_columns(span(0, i))));				
			}
			
			vec x = log(cn); //1/sqrt(cn);
			
			int state_no = 0;
			if(cond_threshold == -1) {
				double x1 = sum(x);
				double x2 = sum(square(x));
				int nnz = ceil(x1*x1 / x2);			
				state_no = x.n_elem - nnz + 1;
			} else if(cond_threshold == -2) {
				double x1 = sum(square(x));
				double x2 = sum(square(square(x)));
				int nnz = ceil(x1*x1 / x2);			
				state_no = x.n_elem - nnz + 1;
			} else {
				state_no = sum(cn <= cond_threshold);
			}			
			printf("Selected %d states\n", state_no);
							
			mat Aa = G_red.cols(selected_columns(span(0, state_no-1)));
			mat cc = run_simplex_regression(Aa, G_red, false);
			uvec idx = find(cc < 0.5);
			cc(idx).zeros();
			
			mat w = trans(cc);
			w = normalise(w, 1, 0);		

			C_unified = C_norm * w;
		} else if(method_type == 2) {
			printf("Running SPA (G:%d x %d, sim = %d)\n", G_red.n_rows, G_red.n_cols, dim); fflush(stdout);
			SPA_results res = run_SPA(G_red, dim);
			
			uvec selected_columns = res.selected_columns;
			
			printf("Computing condition numbers\n"); fflush(stdout);
			vec cn = zeros(dim-1);
			for(int i = 1; i < dim; i++) {
				cn(i-1) = cond(G_red.cols(selected_columns(span(0, i))));
				//cn(i-1) = cond(W_stacked.cols(selected_columns(span(0, i))));				
			}
			
			vec x = log(cn); //1/sqrt(cn);
			
			int state_no = 0;
			if(cond_threshold == -1) {
				double x1 = sum(x);
				double x2 = sum(square(x));
				int nnz = ceil(x1*x1 / x2);			
				state_no = x.n_elem - nnz + 1;
			} else if(cond_threshold == -2) {
				double x1 = sum(square(x));
				double x2 = sum(square(square(x)));
				int nnz = ceil(x1*x1 / x2);			
				state_no = x.n_elem - nnz + 1;
			} else {
				state_no = sum(cn <= cond_threshold);
			}			
			printf("Selected %d states\n", state_no);
							
					
			C_unified = C_norm.cols(selected_columns(span(0, state_no-1)));				
		} else if(method_type == 3) {			
			field<vec> hd_out = run_HDBSCAN(G_red, minPoints, minClusterSize);
			
			uvec selected_archs = find(hd_out(2) <= outlier_threshold); 
			vec clusters = hd_out(0);
			clusters = clusters(selected_archs);
			
			vec uc = sort(unique(clusters));				
			mat subC;
			C_unified = zeros(C_stacked.n_rows, uc.n_elem);
			vec trivial_mask = zeros(uc.n_elem);		
			for(int i = 0; i < uc.n_elem; i++) {			
				if(uc(i) == 0) {
					trivial_mask(i) = 1;
					continue;
				}
				
				uvec ii = find(clusters == uc(i));
				uvec idx = selected_archs(ii);
				
				if(min_repeat <= idx.n_elem) {
					subC = C_norm.cols(idx);
					C_unified.col(i) = mean(subC, 1);	
				} else {
					trivial_mask(i) = 1;
				}
				
			}
			C_unified = C_unified.cols(find(trivial_mask == 0));
					
			C_unified = normalise(C_unified, 1, 0);					
		} else if(method_type == 4) {
			G.diag().zeros();
			
			uvec initial_clusters(G.n_rows);
			for (int i = 0; i < G_red.n_rows; i++) initial_clusters(i) = i;					
			vec clusters = unsigned_cluster(sp_mat(G), resolution, initial_clusters, 0);
			
			uvec selected_archs = find(clusters >= 0);
			clusters = clusters(selected_archs);
			
			
			vec uc = sort(unique(clusters));				
			mat subC;
			C_unified = zeros(C_stacked.n_rows, uc.n_elem);
			vec trivial_mask = zeros(uc.n_elem);		
			for(int i = 0; i < uc.n_elem; i++) {			
				if(uc(i) == 0) {
					trivial_mask(i) = 1;
					continue;
				}
				
				uvec ii = find(clusters == uc(i));
				uvec idx = selected_archs(ii);
				
				if(min_repeat <= idx.n_elem) {
					subC = C_norm.cols(idx);
					C_unified.col(i) = mean(subC, 1);	
				} else {
					trivial_mask(i) = 1;
				}
				
			}
			C_unified = C_unified.cols(find(trivial_mask == 0));
					
			C_unified = normalise(C_unified, 1, 0);			
		}
		mat W_unified = S_r * C_unified;
		mat H_unified = run_simplex_regression(W_unified, S_r, false);			
		

		uvec assigned_archetypes = trans(index_max( H_unified, 0 ));
		
		output.C_unified = C_unified;
		output.H_unified = H_unified;
		output.assigned_archetypes = assigned_archetypes;
		
		//output.dag_adj = dag_adj;
		//output.dag_node_annotations = dag_node_annotations;
		//output.selected_archetypes = selected_archs;		

		return(output);
	}
}
