#include "ACTIONet.h"
#include "dagConstruct.h"

vector<double> Corrector::vals;

	template<class Function>
	inline void ParallelFor(size_t start, size_t end, size_t thread_no, Function fn) {
		if (thread_no <= 0) {
			thread_no = std::thread::hardware_concurrency();
		}

		if (thread_no == 1) {
			for (size_t id = start; id < end; id++) {
				fn(id, 0);
			}
		} else {
			std::vector<std::thread> threads;
			std::atomic<size_t> current(start);

			// keep track of exceptions in threads
			// https://stackoverflow.com/a/32428427/1713196
			std::exception_ptr lastException = nullptr;
			std::mutex lastExceptMutex;

			for (size_t threadId = 0; threadId < thread_no; ++threadId) {
				threads.push_back(std::thread([&, threadId] {
					while (true) {
						size_t id = current.fetch_add(1);

						if ((id >= end)) {
							break;
						}

						try {
							fn(id, threadId);
						} catch (...) {
							std::unique_lock<std::mutex> lastExcepLock(lastExceptMutex);
							lastException = std::current_exception();
							/*
							 * This will work even when current is the largest value that
							 * size_t can fit, because fetch_add returns the previous value
							 * before the increment (what will result in overflow
							 * and produce 0 instead of current + 1).
							 */
							current = end;
							break;
						}
					}
				}));
			}
			for (auto &thread : threads) {
				thread.join();
			}
			if (lastException) {
				std::rethrow_exception(lastException);
			}
		}

	}



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
										double alpha = 0.99,
										double outlier_threshold = -1.65,
										double sim_threshold = 3.0,
										int thread_no = 0) {
		printf("Unify archetypes (%d archs, alpha = %f, outlier_threshold = %f, sim_threshold = %f)\n", C_stacked.n_cols, alpha, outlier_threshold, sim_threshold);
														
		unification_results output;
		
		// Smooth archetypes using ACTIONet
		C_stacked = normalise(C_stacked, 1, 0);
		sp_mat X0 = sp_mat(C_stacked);		
		mat C_imputed = compute_network_diffusion(G, X0, thread_no, alpha, 5);

		// Compute similarity matrix		
		/*
		mat logC = log(C_imputed);
		mat Z = zeros(size(logC));
		for(int j = 0; j < logC.n_cols; j++) {
			vec v = logC.col(j);
			uvec idx = find( (v != datum::nan) && (v != -datum::inf) && (v != +datum::inf));
			double mu = mean(v(idx));
			double sigma = stddev(v(idx));
			vec z = (v-mu)/sigma;
			uvec mask_idx = find( (v == datum::nan) || (v == -datum::inf) || (v == +datum::inf));
			z(mask_idx).zeros();
			Z.col(j) = z;
		}
		
		vec mu = mean(Z, 1);		
		mat R = Z - ( (mu * (trans(mu) * Z)) / sum(square(mu)));
		*/
	/*	
		mat Z = zscore(C_imputed);
		
		//mat Sim = trans(Z) * Z / (double)Z.n_rows;

		
		vec v = trans(sum(square(Z)));			
		mat denom = sqrt(v * trans(v));
		denom = sum(sum(G)) * denom / (double)G.n_rows;
		
		mat Sim = (trans(Z) * G * Z) / denom;

		Sim(find(denom == 0)).zeros();
*/


		mat Sim = cor(C_imputed);
	
/*		
		vec resolutions = regspace(0.1, 0.05, 5);
		mat M_total = signed_cluster_batch(sp_mat(Sim), resolutions);
		Sim = (M_total * trans(M_total));
		
		Sim /= max(max(abs(Sim)));
*/		
		
/*		
		output.C_unified = Sim;
		output.H_unified = Z;		
		
		return(output);
*/

		// Prune unreliable archetypes		
/*
		vec zz = zscore(Sim.diag());
		uvec selected_archetypes = find( (zz > outlier_z_threshold) && (Sim.diag() != 0) );
*/
		Sim.transform( [sim_threshold](double val) { return (val < sim_threshold?0:val); } );
		sp_mat Sim_sp = sp_mat(Sim);
		uvec core_num = compute_core_number(Sim_sp);
		
		uvec selected_archetypes = find(outlier_threshold <= core_num);				
		printf("%d selected archs\n", selected_archetypes.n_elem);		


		// Subset archetypes and compute reduced profile of each archetype
		//M_total = M_total.cols(selected_archetypes);
		
		C_imputed = C_imputed.cols(selected_archetypes);
		mat S_r_arch = S_r * C_imputed;
		
		Sim = Sim(selected_archetypes, selected_archetypes);
		
		// Prioritize archetypes
		int dim = min(100, (int)min(C_imputed.n_cols, C_imputed.n_rows));		
		SPA_results res = run_SPA(C_imputed, dim);		
		
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
						
			//printf("%d- %d %f\n", i+1, j+1, mm);
			if(sim_threshold < mm) {
				continue;
			}
			is_selected(j) = 1;

		}		
		uvec idx = find(is_selected == 1);		
		
		printf("%d unfied archetypes\n", idx.n_elem);
		
		mat C_unified = C_imputed.cols(idx);
		mat W_unified = S_r_arch.cols(idx);
				
				
		mat H_unified = run_simplex_regression(W_unified, S_r, false);			
		uvec assigned_archetypes = trans(index_max( H_unified, 0 ));


		// Construct Ontology!		
		graph_undirected inputNetwork(Sim);
		DAGraph ontology;
		ontology.setTerminalName("archetype");
		nodeDistanceObject nodeDistances;
		
		double threshold =  0.05;
		double density = 0.5;
  
		dagConstruct::constructDAG(inputNetwork, ontology, nodeDistances, threshold, density);
		
		vector<vector<int>> dag_nodes;
		vector<int> dag_nodes_type;
		
		for(int k = 0; k < Sim.n_rows; k++) {
			vector<int> v;
			v.push_back(k);
			dag_nodes.push_back(v);
			dag_nodes_type.push_back(0);
		}
				
		for(map< pair<unsigned,unsigned>, string >::iterator edgesIt = ontology.edgesBegin(); edgesIt != ontology.edgesEnd(); ++edgesIt) {
			unsigned ii = edgesIt->first.first;
			unsigned jj = edgesIt->first.second;

			vector<int> v;			
			int tt;
			if(edgesIt->second == "archetype") {
				v.push_back(jj);
				tt = 0;
			} else { // Internal node
				for(int kk = 0; kk < dag_nodes[jj].size(); kk++) {
					v.push_back(dag_nodes[jj][kk]);
				}
				tt = 1;
			}
			
			if(ii >= dag_nodes.size()) {
				dag_nodes.push_back(v);
				dag_nodes_type.push_back(tt);
			} else { // merge
				for(int kk = 0; kk < v.size(); kk++) {
					dag_nodes[ii].push_back(v[kk]);
				}	
				
				if(edgesIt->second != "archetype")
					dag_nodes_type[ii] = 1;
			}
			
			//cout << ontology.getName(edgesIt->first.first) << "\t" << ontology.getName(edgesIt->first.second) << "\t" << edgesIt->second << "\t" << ontology.getWeight(edgesIt->first.first) << endl;
		}
		
		// Get internal adjacency matrix of DAGs
		int dag_node_counts = dag_nodes.size();
		mat dag_adj = zeros(dag_node_counts, dag_node_counts);	
		vec dag_node_annotations = zeros(dag_node_counts);		
		for(map< pair<unsigned,unsigned>, string >::iterator edgesIt = ontology.edgesBegin(); edgesIt != ontology.edgesEnd(); ++edgesIt) {
			unsigned ii = edgesIt->first.first;
			unsigned jj = edgesIt->first.second;
			double w = ontology.getWeight(edgesIt->first.first);

			if(edgesIt->second == "archetype") {
				dag_node_annotations(ii) = 1;
			} else {
				dag_node_annotations(ii) = 2;
			}
			
			dag_adj(ii, jj) = 1;
			
		}
		
		output.dag_adj = dag_adj;
		output.dag_node_annotations = dag_node_annotations;
		
		output.C_unified = C_unified;
		output.H_unified = H_unified;
		output.assigned_archetypes = assigned_archetypes;
		output.selected_archetypes = selected_archetypes(idx);		
		
		return(output);
	}
	
}
