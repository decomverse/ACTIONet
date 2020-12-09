#include <ACTIONet.h>

namespace ACTIONet {
	multilevel_archetypal_decomposition prune_archetypes(field<mat> C_trace, field<mat> H_trace, double min_specificity_z_threshold, int min_cells = 3) {
		mat C_stacked;
		mat H_stacked;
		int depth = H_trace.size();

		multilevel_archetypal_decomposition results;

		// Vector contains an element for k==0, this have to -1
		Rprintf("Joining the trace of C& H matrices (depth = %d) ... ", depth-1); //fflush(stdout);
		// Group H and C matrices for different values of k (#archs) into joint matrix
		for(int k = 0; k < depth; k++) {
			if(H_trace[k].n_rows == 0)
				continue;

			if(H_stacked.n_elem == 0) {
				C_stacked = C_trace[k];
				H_stacked = H_trace[k];

			}
			else {
				C_stacked = join_rows(C_stacked, C_trace[k]);
				H_stacked = join_cols(H_stacked, H_trace[k]);
			}
		}
		int total_archs = H_stacked.n_rows;
		Rprintf("done\n");



		Rprintf("Pruning non-specific archetypes (based on transitivity) ... ");
		// Construct backbone
		mat backbone = cor(trans(H_stacked));
		backbone.diag().zeros();
		backbone.transform( [](double val) { return (val< 0?0:val); } );

		vec pruned = zeros(total_archs);

		// Barrat weighted transitivity: formulation from "Clustering Coefficients for Weighted Networks" (Kalna)
		vec transitivity = zeros(total_archs);
		vec s = sum(backbone, 1); // strength of nodes
		vec d = vec(sum(spones(sp_mat(backbone)), 1));
		for(int k = 0; k < total_archs; k++) {

			double sum = 0;
			for(int i = 0; i < total_archs; i++) {
				double w_ki = backbone(k, i);
				for(int j = 0; j < total_archs; j++) {
					double w_kj = backbone(k, j);

					double mean_weight = (w_ki + w_kj) / 2.0;
					double triangle_mask = backbone(i, k)*backbone(k, j)*backbone(j, i) > 0?1:0;

					sum += mean_weight*triangle_mask;
				}
			}
			transitivity(k) = sum / (s(k)*(d(k)-1));
		}

		vec transitivity_z = zscore(transitivity);
		uvec nonspecific_idx = find(transitivity_z < min_specificity_z_threshold);
		pruned(nonspecific_idx).ones();
		Rprintf("done (%d archs pruned)\n", nonspecific_idx.n_elem);


		// Find landmark cells, i.e., closest cells to each multi-level archetype (its projection on to the cell space, ish)
		Rprintf("Removing unreliable archetypes (based on the landmark cells) ... ");
		double epsilon = 1e-6;
		int bad_archs = 0;
		vec landmark_cells = -ones(total_archs);
		for(int i = 0; i < total_archs; i++) {
			vec h = trans(H_stacked.row(i));
			vec c = C_stacked.col(i);

			uvec h_landmarks = find( (max(h) - h) < epsilon );
			uvec c_landmarks = find( 0 < c );
			uvec common_landmarks = intersect(h_landmarks, c_landmarks);

			if(0 < common_landmarks.n_elem) { // They don't agree on any samples!
				landmark_cells(i) = common_landmarks(index_max(c(common_landmarks)));
			}
			else { // Potentially noisy archetype
				pruned(i) = 1;
				bad_archs++;
			}
		}
		Rprintf("done (further %d archs removed)\n", bad_archs); fflush(stdout);


		uvec idx = find(C_stacked > 0);
		mat C_bin = C_stacked;
		C_bin(idx).ones();
		uvec trivial_idx = find(sum(C_bin) < min_cells);
		pruned(trivial_idx).ones();

		Rprintf("Found (and removed) %d trivial archetypes\n", trivial_idx.n_elem); fflush(stdout);


		uvec selected_archs = find(pruned == 0);
		results.selected_archs = selected_archs;
		results.C_stacked = C_stacked.cols(selected_archs);
		results.H_stacked = H_stacked.rows(selected_archs);


		return(results);
	}
}
