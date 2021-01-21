#include <ACTIONet.h>

namespace ACTIONet {
	
mat assess_label_enrichment(sp_mat &H, sp_mat &M) {
	mat Obs = mat(M * H);
	
	vec p = vec(mean(M, 1));
	mat Exp = (p * sum(H));
	
	mat Lambda = Obs - Exp;
	
	mat Nu = (p * sum(square(H)));
	vec a = vec(trans(max(H, 0)));

	mat Lambda_scaled = Lambda;
	for (int j = 0; j < Lambda_scaled.n_cols; j++) {
		Lambda_scaled.col(j) *= (a(j) / 3);
	}
	mat logPvals_upper = square(Lambda) / (2 * (Nu + Lambda_scaled));
	uvec lidx = find(Lambda <= 0);
	logPvals_upper(lidx) = zeros(lidx.n_elem);
	logPvals_upper.replace(datum::nan, 0);  // replace each NaN with 0


	return logPvals_upper;
}

sp_mat one_hot_encoding(vec labels) {
	int n = labels.n_elem;
	
	vec vals = unique(labels);
	uvec idx = find(0 <= vals);
	vals = vals(idx);
	
	int k = vals.n_elem;
	sp_mat M(k, n);
	for(int i = 0; i < k; i++) {
		uvec idx = find(labels == vals(i));
		for(int j = 0; j < idx.n_elem; j++) {
			M(i, idx(j)) = 1;
		}
	}
	
	return(M);
}

vec LPA(sp_mat &G, vec labels, double lambda = 0, int iters = 3, double sig_threshold = 3, uvec fixed_labels = uvec()) {

	int n = G.n_rows;
	
	vec updated_labels = labels;	
	
	sp_mat H = G;
	H.diag().ones();
	H.diag() *= lambda; // add "inertia"
	H = n*normalise(H, 1, 0); // column-normalize to n
	
	
	for(int it = 0; it < iters; it++) {
		vec vals = unique(updated_labels);
		uvec idx = find(0 <= vals);
		vals = vals(idx);
				
		sp_mat M = one_hot_encoding(updated_labels);

		
		mat logPvals = assess_label_enrichment(H, M);

		
		vec  max_sig = trans(max(logPvals, 0));
		vec new_labels = vals(trans(index_max(logPvals, 0)));
		
		// Only update vertices with significant enrichment in their neighborhood
		uvec sig_idx = find(sig_threshold < max_sig);
		updated_labels(sig_idx) = new_labels(sig_idx);

		// revert-back the vertices that need to be fixed
		updated_labels(fixed_labels) = labels(fixed_labels);		
	}
	
	return(updated_labels);
	
}

/*
vec Kappa(vec p, vec q) {
    vec kl = ones(size(p));
	
	vec a = p * log(p / q);	
	uvec zero_idx = find(p == 0);
	a(zero_idx).zeros();
    
	vec b = (1.0 - p) * log((1.0 - p)/(1.0 - q));
	uvec one_idx = find(p == 1);
	b(one_idx).zeros();	
    
    vec k = a + b;
    
    return(k);
}

HGT_tail(vec success_counts, vec sample_size, observed.success, int population_size) {
    if (sum(success_count) == 0) {
		return(zeros(size(success_counts)));
	} 
	
    vec success_rate = success_count/population_size;
    vec expected_success = sample.size % success.rate;
    
    delta = (observed_success/expected.success) - 1
    
    log.tail_bound = sample.size * Kappa((1 + delta) * success.rate, 
        success.rate)
    log.tail_bound[delta < 0] = 0
    log.tail_bound[is.na(log.tail_bound)] = 0
    return(log.tail_bound)
}
*/

}
