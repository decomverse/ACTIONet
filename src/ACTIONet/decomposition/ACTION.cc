#include "ACTIONet.h"

namespace ACTIONet {
	// Solves the standard Archetypal Analysis (AA) problem
	field<mat> run_AA(mat &A, mat &W0, int max_it = 50, double min_delta = 0.01) {	
		
		int sample_no = A.n_cols;
		int d = A.n_rows; // input dimension
		int k = W0.n_cols; // AA components
		   
		
		mat C = zeros(sample_no, k);
		mat H = zeros(k, sample_no);
		
		mat W = W0;
		vec c(sample_no);	

		//printf("(New) %d- %d\n", k, max_it);

		for (int it = 0; it < max_it; it++) {		
			H = run_simplex_regression(W, A);
			
			mat C_old = C;
			mat R = A - W*H;	
			for(int i = 0; i < k; i++) {
				vec w = W.col(i);
				rowvec h = H.row(i);
				vec ht = trans(h);
				
				double norm_sq = arma::dot(h, h);
				if(norm_sq < double(10e-8)) {
					printf("%d- Ooops @<%d, %d>\n",k, it, i);
					
					// singular
					int max_res_idx = index_max(rowvec(sum(square(R), 0)));
					W.col(i) = A.col(max_res_idx);
					c.zeros();
					c(max_res_idx) = 1;
					C.col(i) = c;					
				} else {				 				
					vec b = (1.0 / norm_sq) *R*ht + w;
					c = run_simplex_regression(A, b);		
					R += (w - A*c)*h; // Rank-1 update
					W.col(i) = A*c;
					C.col(i) = c;			
				}
			}				
			double delta = arma::max(rowvec(sum(abs(C - C_old)))) / 2.0;

			//double RSS = norm(R, "fro"); RSS *= RSS;
			//printf("\t<%d, %d>- RSS = %.3e, delta = %.3e\n", l, it, RSS, delta);
			
			if(delta < min_delta)
				break;
		}
		

		C = clamp(C, 0, 1);
		C = normalise(C, 1);
		H = clamp(H, 0, 1);
		H = normalise(H, 1);
		
		field<mat> decomposition(2,1);
		decomposition(0) = C;
		decomposition(1) = H;
		
		return decomposition;
	}

	template<class Function>
	inline void ParallelFor(size_t start, size_t end, size_t numThreads, Function fn) {
		if (numThreads <= 0) {
			numThreads = std::thread::hardware_concurrency();
		}

		if (numThreads == 1) {
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

			for (size_t threadId = 0; threadId < numThreads; ++threadId) {
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
	
	// Solves separable NMF problem
	SPA_results run_SPA(mat A, int k) {	
		
		SPA_results res;
		
		int n = A.n_cols;
		uvec K(k); // selected columns from A
			
	 
		rowvec normM = sum(A % A, 0); 
		rowvec normM1 = normM;
		
		mat U(A.n_rows, k);
		
		vec norm_trace = zeros(k);
		double eps = 1e-6;
		for (int i = 0; i < k; i++) {		
			// Find the column with maximum norm. In case of having more than one column with almost very small diff in norm, pick the one that originally had the largest norm
			double a = max(normM); 
			norm_trace(i) = a;
			
			uvec b = find((a*ones(1, n)-normM)/a <= eps); 
			
			if(b.n_elem > 1) {
				uword idx = index_max(normM1(b)); 
				K(i) = b(idx);
			}
			else {
				K(i) = b(0);			
			}			
			
			// Pick column
			U.col(i) = A.col(K(i));

			// Orthogonalize with respect to current basis
			for (int j = 0; j < i-1; j++) {
				U.col(i) = U.col(i) - dot(U.col(j), U.col(i)) * U.col(j);
			}
			U.col(i) = U.col(i)/ norm(U.col(i), 2); 
			
			// Update column norms
			vec u = U.col(i);            
			for (int j = i-1; 0 <= j; j--) {
				u = u - dot(U.col(j), u)*U.col(j); 
			}
			rowvec r = u.t()*A;
			normM = normM - (r % r);
		}
			
		res.selected_columns = K;
		res.column_norms = norm_trace;
		
		return res;
	}

	ACTION_results run_ACTION(mat S_r, int k_min, int k_max, int thread_no, int max_it = 50, double min_delta = 0.01) {
		int feature_no = S_r.n_rows;
				
		printf("Running ACTION\n");
		
		if(k_max == -1)
			k_max = (int)S_r.n_cols;
			
		k_min = std::max(k_min, 2);
		k_max = std::min(k_max, (int)S_r.n_cols);	
					
		ACTION_results trace; 
		/*
		trace.H.resize(k_max + 1);
		trace.C.resize(k_max + 1);
		trace.selected_cols.resize(k_max + 1);
		*/

		trace.H = field<mat>(k_max + 1);
		trace.C = field<mat>(k_max + 1);
		trace.selected_cols = field<uvec>(k_max + 1);


		
		mat X_r = normalise(S_r, 1); // ATTENTION!
		 		
		int current_k = 0;		
		int total = k_min-1;
		printf("Iterating from k=%d ... %d\n", k_min, k_max);
		ParallelFor(k_min, k_max+1, thread_no, [&](size_t kk, size_t threadId) {			
			total++;
			printf("\tk = %d\n", total);
			SPA_results SPA_res = run_SPA(X_r, kk);
			trace.selected_cols[kk] = SPA_res.selected_columns;
			
			mat W = X_r.cols(trace.selected_cols[kk]);
			
			field<mat> AA_res;			
			AA_res = run_AA(X_r, W, max_it, min_delta);

			trace.C[kk] = AA_res(0);
			trace.H[kk] = AA_res(1);			
		});
		
		return trace;
	}	

	ACTION_results run_ACTION_dev(mat S_r, int k_min, int k_max, int thread_no, bool auto_stop = true, int max_it = 30, double min_delta = 0.01) {
		int feature_no = S_r.n_rows;
		
		printf("Running ACTION (developmental version)\n");
		
		if(k_max == -1)
			k_max = (int)S_r.n_cols;
			
		k_min = std::max(k_min, 2);
		k_max = std::min(k_max, (int)S_r.n_cols);	
					
		ACTION_results trace; 
		/*
		trace.H.resize(k_max + 1);
		trace.C.resize(k_max + 1);
		trace.selected_cols.resize(k_max + 1);
		*/

		trace.H = field<mat>(k_max + 1);
		trace.C = field<mat>(k_max + 1);
		trace.selected_cols = field<uvec>(k_max + 1);
		
		
		mat X_r = normalise(S_r, 1); // ATTENTION!
		 		
		int current_k = 0;				
		printf("Iterating from k=%d ... %d (auto stop = %d)\n", k_min, k_max, auto_stop);
		for(int kk = k_min; kk <= k_max; kk++) {
			printf("\tk = %d\n", kk);
			SPA_results SPA_res = run_SPA(X_r, kk);
			trace.selected_cols[kk] = SPA_res.selected_columns;
			
			mat W = X_r.cols(trace.selected_cols[kk]);
			if(kk > k_min) {
				W.cols(span(0, kk-2)) = X_r * trace.C[kk-1];
			}
			//W.print("W");
			
			
			field<mat> AA_res = run_AA(X_r, W, max_it, min_delta);
			
			if(auto_stop) {
				mat C = AA_res(0);
				bool has_trivial_arch = false;
				for(int c = 0; c < C.n_cols; c++) {
					int r = 0, nnz_counts = 0;
					while(r < C.n_rows) {
						if(C(r, c) > 0)
							nnz_counts ++;
							
						if(nnz_counts > 1)
							break;
						
						r++;
					}
					if(nnz_counts <= 1) {
						has_trivial_arch = true;
						break;
					}
				}				
				if(has_trivial_arch > 0) {
					printf("\t\tFound trivial archetypes at k = %d\n", kk);
					
					break;
				}
			}
			
			current_k = std::max(current_k, kk);

			trace.C[kk] = AA_res(0);
			trace.H[kk] = AA_res(1);			
		}
		/*
		trace.H.resize(current_k+1);
		trace.C.resize(current_k+1);
		trace.selected_cols.resize(current_k+1);
		*/
		

		return trace;
	}	
}
