#include <ACTIONet.h>

namespace ACTIONet {
	mat assess_enrichment(mat &scores, mat &associations, int L) {			
		field<mat> res(3);

		if(scores.n_rows != associations.n_rows) {
			fprintf(stderr, "Number of rows in scores and association matrices should both match the number of features\n");
		}

		
		mat sorted_scores = sort(scores, "descend");
		mat top_sorted_scores = sorted_scores.rows(span(0, L-1));


		vec p_success = vec(trans(mean(associations, 0)));
		vec A_sum = vec(trans(sum(scores, 0)));
		vec A_sum_sq = vec(trans(sum(square(scores), 0)));

		vec f = zeros(L);
		for(int i = 0; i < L; i++) {
			f(i) = (i+1.0)/(double)scores.n_rows;
		}
		vec f_sq = square(f);
		
		//vec max_L =min(associations.n_rows*ones(associations.n_rows, 1), 3*sum(associations, 0));
		
		mat logPvals = zeros(scores.n_cols, associations.n_cols);
		for(int i = 0; i < scores.n_cols; i++) {	
			
			vec x = scores.col(i);
			uvec perm = sort_index(x, "descend");
			perm = perm(span(0, L-1));
			
			
			vec a = top_sorted_scores.col(i);					
			mat X = associations.rows(perm);
			
			vec a_cs = cumsum(a);
			vec a_cs_sq = cumsum(square(a));
			double a_max = max(a);
			for(int j = 0; j < associations.n_cols; j++) {
				vec x = X.col(j);						
				vec ax = a % x;
				vec o = cumsum(ax);
				
				
				vec e = a_cs * p_success(j);
				vec nu = a_cs_sq * p_success(j);
				
				//vec e = (A_sum(i) * p_success(j)) * f;
				//vec nu = (A_sum_sq(i) * p_success(j)) * f_sq;
				
				vec lambda = o - e;
				
				vec logpvals = square(lambda) / (2.0 * (nu + (a_max*lambda/3.0)));
				uvec idx = find(lambda <= 0);
				logpvals(idx) = zeros(idx.n_elem);
				
				logPvals(i, j) = max(logpvals);
			}
		}
		
		logPvals /= log(10);
		
		return(logPvals);
	}		
	
}
