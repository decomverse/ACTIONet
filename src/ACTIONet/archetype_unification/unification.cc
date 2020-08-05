#include "ACTIONet.h"
#include "dagConstruct.h"

vector<double> Corrector::vals;

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
	
	
	mat NetEnh(mat A) {		
		A.diag().zeros();
		mat P = normalise(A, 1, 1);
		
		mat D = diagmat(sqrt(sum(P) + 1e-16));
		mat W = P * D;		
		mat P2 = W * trans(W);		
		//P2.diag().zeros();
		
		return(P2);
	}

	unification_results unify_archetypes(mat &S_r, mat &C_stacked, mat &H_stacked, double min_edge_weight = 0.5, int min_coreness = 2, double resolution = 1.0, int min_repeat = 2, int thread_no = 0, double alpha = 0.05, double beta = 0.5) {
		printf("Unify archetypes: resolution = %.2f, min. weight = %.2f, alpha = %.2f, beta = %.2f\n", resolution, min_edge_weight, alpha, beta);
		
		unification_results output;

		mat C_norm = normalise(C_stacked, 1, 0);
		mat arch_H_stacked = H_stacked * C_norm;
		mat G = computeFullSim(arch_H_stacked, thread_no);

		//sp_mat G = build_ACTIONet(arch_H_stacked, density, thread_no, M, ef_construction, ef, mutual_edges_only);
		
		// Find outlier archetypes
		//G = clamp(G, min_edge_weight, 1); 		
		G.transform( [min_edge_weight](double val) { return (val < min_edge_weight? 0:val); } );

		sp_mat Gs = sp_mat(G);
		uvec cn = ACTIONet::compute_core_number(Gs);
		uvec selected_archs = find(min_coreness <= cn); // Basically ignore for now!

		
		//mat G_enh = NetEnh(mat(G));		
		sp_mat subG = sp_mat(G(selected_archs, selected_archs));
		
		uvec initial_clusters(subG.n_rows);
		for (int i = 0; i < subG.n_rows; i++) initial_clusters(i) = i;		
		
		
		mat M_total;	
		vec  clusters = conv_to<vec>::from(initial_clusters);
		for(double res = 10; res >= 0.1; res -= 0.1) {
			initial_clusters = conv_to<uvec>::from(clusters);
			clusters = unsigned_cluster(subG, res, initial_clusters, 0);
			mat M = zeros(clusters.n_elem, max(clusters)+1);
			for(int i = 0; i < clusters.n_elem; i++) {
				M(i, clusters(i)) = 1;
			}
			if(M_total.n_elem == 0) {
				M_total = M;
			} else {
				M_total = join_horiz(M_total, M);
			}
		}			
		mat new_sim = M_total * trans(M_total);
		new_sim.diag().zeros();

		new_sim = NetEnh(new_sim);

		
		//vec clusters = unsigned_cluster(subG, resolution, initial_clusters, 0);
		
		graph_undirected inputNetwork(new_sim);
		DAGraph ontology;
		ontology.setTerminalName("archetype");
		nodeDistanceObject nodeDistances;

		double threshold = alpha;
		double density = beta;
  
		dagConstruct::constructDAG(inputNetwork, ontology, nodeDistances, threshold, density);
		
		vector<vector<int>> dag_nodes;
		vector<int> dag_nodes_type;
		
		for(int k = 0; k < new_sim.n_rows; k++) {
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
		
/*		
		for(int k = 0; k < dag_nodes.size()-1; k++) {
			if(dag_nodes[k].size() > 1) {
				printf("%d- ", k);
				for(int kk = 0; kk < dag_nodes[k].size(); kk++) {
					printf("%d ", dag_nodes[k][kk]);
				}
				printf("\n");
			}
		}
*/

		
		mat C_unified = zeros(C_stacked.n_rows, dag_nodes.size() - 1 - new_sim.n_rows);
		vec node_type = zeros(dag_nodes.size() - 1 - new_sim.n_rows);
		for(int k = new_sim.n_rows; k < dag_nodes.size()-1; k++) {			
			//printf("%d (%d)- size = %d, type = %d\n", k-new_sim.n_rows + 1, k, dag_nodes[k].size(), dag_nodes_type[k]);
			if(dag_nodes[k].size() <= min_repeat) {
				node_type(k-new_sim.n_rows) = 1;			
				continue;
			}

			
			node_type(k-new_sim.n_rows)	= dag_nodes_type[k];			
			for(int kk = 0; kk < dag_nodes[k].size(); kk++) {
				int idx = selected_archs(dag_nodes[k][kk]);				
				C_unified.col(k-new_sim.n_rows) += C_norm.col(idx);				
			}
			C_unified.col(k-new_sim.n_rows) /= dag_nodes[k].size();
		}
		C_unified = C_unified.cols(find(node_type == 0));
		
/*

		clusters = unsigned_cluster(sp_mat(new_sim), resolution, initial_clusters, 0);
		vec full_clusters = -ones(G.n_rows);
		vec uc = sort(unique(clusters));		
		
		mat subC;
		mat C_unified(C_stacked.n_rows, uc.n_elem);
		vec trivial_mask = zeros(uc.n_elem);		
		for(int i = 0; i < uc.n_elem; i++) {			
			uvec ii = find(clusters == uc(i));
			uvec idx = selected_archs(ii);
			
			if(min_repeat <= idx.n_elem) {
				subC = C_norm.cols(idx);
				C_unified.col(i) = mean(subC, 1);	
				full_clusters(idx) = i*ones(idx.n_elem);			
			} else {
				trivial_mask(i) = 1;
			}
			
		}
		C_unified = C_unified.cols(find(trivial_mask == 0));
*/	
				
		C_unified = normalise(C_unified, 1, 0);
		
		mat W_unified = S_r * C_unified;
		mat H_unified = run_simplex_regression(W_unified, S_r, false);
		
		uvec assigned_archetypes = trans(index_max( H_unified, 0 ));


		
		output.C_unified = C_unified;
		output.H_unified = H_unified;
		output.assigned_archetypes = assigned_archetypes;
		
		output.dag_adj = dag_adj;
		output.dag_node_annotations = dag_node_annotations;
		output.selected_archetypes = selected_archs;		

		return(output);
	}
}
