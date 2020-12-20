#include <ACTIONet.h>

namespace ACTIONet {
	mat compute_pseudo_bulk_per_cluster(sp_mat &S, arma::Col<unsigned long long> sample_assignments) {	
		
		mat pb = zeros(S.n_rows, max(sample_assignments));

		sp_mat::const_iterator it     = S.begin();
		sp_mat::const_iterator it_end = S.end();
		for(; it != it_end; ++it) {
			int i = it.row();
			int j = sample_assignments[it.col()]-1;
			pb(i, j) += (*it);
		}
		
		for(int j = 0; j < pb.n_cols; j++) {
			uvec idx = find(sample_assignments == (j+1));
			pb.col(j) /= idx.n_elem;
		}
		
		return(pb);
	}	

	mat compute_pseudo_bulk_per_cluster(mat &S, arma::Col<unsigned long long> sample_assignments) {			
		mat pb = zeros(S.n_rows, max(sample_assignments));
		
		for(int j = 0; j < pb.n_cols; j++) {
			uvec idx = find(sample_assignments == (j+1));
			if(idx.n_elem == 0)
				continue;
				
			if(idx.n_elem > 1) {
				mat subS = S.cols(idx);
				pb.col(j) = mean(subS, 1);						
			} else {
				pb.col(j) = S.col(idx(0));
			}
		}
		
		return(pb);
	}		
	
	
	
	mat compute_pseudo_bulk_per_archetype(sp_mat &S, mat& H) {	
		mat H_norm = trans(H);
		
		vec col_means = mean(H, 0);
		for(int i = 0; i < H_norm.n_cols; i++) {
			H_norm.row(i) /= col_means(i);
		}
		
		mat pb = mat(S * H_norm);
		
		return(pb);
	}	

	mat compute_pseudo_bulk_per_archetype(mat &S, mat& H) {			
		mat H_norm = trans(H);
		
		vec col_means = mean(H, 0);
		for(int i = 0; i < H_norm.n_cols; i++) {
			H_norm.row(i) /= col_means(i);
		}
		
		mat pb = mat(S * H_norm);
		
		return(pb);
	}		
		
	
	
	
	field<mat> compute_pseudo_bulk_per_ind(sp_mat &S, arma::Col<unsigned long long> sample_assignments, arma::Col<unsigned long long> individuals) {	
		
		field<mat> pbs(max(sample_assignments));
		for(int k = 0; k < max(sample_assignments); k++) {
			pbs(k) = zeros(S.n_rows, max(individuals));
		}

		sp_mat::const_iterator it     = S.begin();
		sp_mat::const_iterator it_end = S.end();
		for(; it != it_end; ++it) {
			int i = it.row();
			int j = individuals[it.col()]-1;
			int k = sample_assignments[it.col()]-1;
						
			pbs(k)(i, j) += (*it);
		}
		
		for(int j = 0; j < max(individuals); j++) {
			for(int k = 0; k < max(sample_assignments); k++) {				
				uvec idx = intersect( find( (sample_assignments == (k+1)) ), find( (individuals == (j+1)) ) );
				
				pbs(k).col(j) /= idx.n_elem;
			}
		}
		
		return(pbs);
	}	

	field<mat> compute_pseudo_bulk_per_ind(mat &S, arma::Col<unsigned long long> sample_assignments, arma::Col<unsigned long long> individuals) {	
		
		field<mat> pbs(max(sample_assignments));
		for(int k = 0; k < max(sample_assignments); k++) {
			pbs(k) = zeros(S.n_rows, max(individuals));
		}

		for(int j = 0; j < max(individuals); j++) {
			for(int k = 0; k < max(sample_assignments); k++) {				
				uvec idx = intersect( find( (sample_assignments == (k+1)) ), find( (individuals == (j+1)) ) );
				mat subS = S.cols(idx);
				pbs(k).col(j) = mean(subS, 1);
			}
		}
		
		return(pbs);
	}			
	
}
