#include <ACTIONet.h>

#include <thread>
#include <atomic>

#include <limits>
#define DBL_MAX std::numeric_limits<double>::max()

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

mat zscore(mat A) {	
	rowvec mu = mean(A, 0);
	rowvec sigma = stddev(A, 0);

	
	for(int j = 0; j < A.n_cols; j++) {
		A.col(j) = (A.col(j) - mu(j)) / sigma(j);
	}
	
	return A;
}

namespace ACTIONet {	

	sp_mat smoothKNN(sp_mat D, int thread_no = -1) {		
		double epsilon = 1e-6;

		int nV = D.n_rows;
		sp_mat G = D;			
		
		//#pragma omp parallel for num_threads(thread_no) 			
		//for(int i = 0; i < nV; i++) {
		ParallelFor(0, nV, thread_no, [&](size_t i, size_t threadId) {
			sp_mat v = D.col(i);
			vec vals = nonzeros(v);	
			if(vals.n_elem > 0) {
				
				double rho = min(vals);
				vec negated_shifted_vals = -(vals - rho); 
				double target = log2(vals.n_elem);
				
				// Binary search to find optimal sigma
				double sigma = 1.0;
				double lo = 0.0;
				double hi = DBL_MAX;
				
				int j;
				for(j = 0; j < 64; j ++) {
					double obj = sum(exp(negated_shifted_vals / sigma));

					if (abs(obj - target) < epsilon) {
						break;
					}

					if (target < obj) {
						hi = sigma;
						sigma = 0.5 * (lo + hi);
					}
					else {
						lo = sigma;
						if (hi == DBL_MAX) {
							sigma *= 2;
						}
						else {
							sigma = 0.5 * (lo + hi);
						}
					}				
				}
				
				double obj = sum(exp(negated_shifted_vals / sigma));			
				//printf("%d- rho = %.3f, degree = %d, log2(k) = %.2e, sigma = %.2e, residual = %.2e, iters = %d\n", i, rho, vals.n_elem, target, sigma, abs(obj - target), j);
				
				for(sp_mat::col_iterator it = G.begin_col(i); it != G.end_col(i); ++it) {
					*it = max(1e-16, exp( -max(0.0, (*it) - rho ) / sigma ));
				}			
			}
		});
		
		return(G);
	}

	field<mat> layout_ACTIONet(sp_mat &G,
		mat &S_r,
		int compactness_level = 50,
		unsigned int n_epochs = 500,
		int thread_no = 8) { 

		unsigned int nV = G.n_rows;
		
		mat initial_coordinates = S_r.rows(regspace<uvec>(0, 1));		
		
		// Convert back from similarity to distance, and then smooth them using the UMAP framework
		/*
		for(sp_mat::iterator it = G.begin(); it != G.end(); ++it) {
		  (*it) = 1.0 - (*it);
		}
		*/			
		printf("Running layout with: compactness=%d, # epochs = %d\n", compactness_level, n_epochs);
		
		G.for_each( [](sp_mat::elem_type& val) { val = 1.0 - val; } );
		
		G = smoothKNN(G, thread_no);


		if(compactness_level < 0 || compactness_level > 100)
			compactness_level = 50;
			
		double a_param = UMAP_A[compactness_level];
		double b_param = UMAP_B[compactness_level];

		// linearized list of edges (1-simplices)
		unsigned int nE = G.n_nonzero;
		vector<unsigned int> positive_head(nE);
		vector<unsigned int> positive_tail(nE);
		vector<double> epochs_per_sample(nE);		
		
		int i = 0;
		double w_max = max(max(G));
		for(sp_mat::iterator it = G.begin(); it != G.end(); ++ it) {
			epochs_per_sample[i] = w_max / (*it); // Higher the weight of the edge, the more likely it is to be sampled (inversely proportional to epochs_per_sample)
			positive_head[i] = it.row();
			positive_tail[i++] = it.col();
		}		
		
		// Initial coordinates of vertices (0-simplices)
		mat Ct = initial_coordinates;		
		vector<double> head_vec(Ct.memptr(), Ct.memptr()+Ct.n_elem);		
		vector<double> tail_vec(head_vec);		
		
		
		printf("Computing 2D layout ... "); fflush(stdout);		
		// Stores linearized coordinates [x1, y1, x2, y2, ...]
		vector<double> result;
		const apumap_gradient gradient(a_param, b_param, GAMMA);
		result = optimize_layout(
			gradient, head_vec, tail_vec, positive_head, positive_tail,
			n_epochs, nV, epochs_per_sample, LEARNING_RATE,
			NEGATIVE_SAMPLE_RATE, UMAP_SEED);	

		mat coordinates(result.data(), 2, nV);
		coordinates = zscore(trans(coordinates));		
		printf("Done\n"); fflush(stdout);

		
		/****************************
		 *  Now compute node colors *
		 ***************************/	
		Ct = trans(join_horiz(coordinates, zscore(trans(S_r.row(2)))));

		head_vec.clear(); 
		head_vec.resize(Ct.n_elem);
		std::copy(Ct.memptr(), Ct.memptr() + Ct.n_elem, head_vec.begin());
		
		tail_vec = head_vec;		
		
		printf("Compute 3D layout ... "); fflush(stdout);
		result.clear();
		result = optimize_layout(
			gradient, head_vec, tail_vec, positive_head, positive_tail,
			n_epochs, nV, epochs_per_sample, LEARNING_RATE,
			NEGATIVE_SAMPLE_RATE, UMAP_SEED);	

		mat coordinates_3D(result.data(), 3, nV);	
		coordinates_3D = zscore(trans(coordinates_3D));		
		printf("Done\n"); fflush(stdout);
	  
	  
		printf("Estimating node colors ... "); fflush(stdout);
		vec a = 128*coordinates_3D.col(0) / max(abs(coordinates_3D.col(0)));
		vec b = 128*coordinates_3D.col(1) / max(abs(coordinates_3D.col(1)));
		vec L = coordinates_3D.col(2);
		L = 25.0 + 50.0*(L - min(L)) / (max(L) - min(L));

		double r_channel, g_channel, b_channel;
		mat RGB_colors = zeros(nV, 3);
		for(int i = 0; i < nV; i++) {
			Lab2Rgb(&r_channel, &g_channel, &b_channel, L(i), a(i), b(i));			

			RGB_colors(i, 0) = min(1.0, max(0.0, r_channel));
			RGB_colors(i, 1) = min(1.0, max(0.0, g_channel));
			RGB_colors(i, 2) = min(1.0, max(0.0, b_channel));			
		}
		
		printf("done\n");
		
		field<mat> res(3);
		res(0) = coordinates;
		res(1) = coordinates_3D;
		res(2) = RGB_colors;
		
		return res;	  		
	}	
	
	// coors: 2 x N
	mat update_layout_2D(mat &coors,
		int compactness_level = 50,
		unsigned int n_epochs = 500,
		int thread_no = 8) { 
		
		unsigned int N = coors.n_cols;
		
		mat D = zeros(N, N);
		for(int i = 0; i < N; i++) {
			vec xi = coors.col(i);
			for(int j = i+1; j < N; j++) {
				vec xj = coors.col(j);
				
				vec delta = xi - xj;
				D(i, j) = D(j, i) = norm(delta, 2);
			}
		}		
		
		sp_mat G = smoothKNN(sp_mat(D), thread_no);


		if(compactness_level < 0 || compactness_level > 100)
			compactness_level = 50;
			
		double a_param = UMAP_A[compactness_level];
		double b_param = UMAP_B[compactness_level];

		// linearized list of edges (1-simplices)
		unsigned int nE = G.n_nonzero;
		vector<unsigned int> positive_head(nE);
		vector<unsigned int> positive_tail(nE);
		vector<double> epochs_per_sample(nE);		
		
		int i = 0;
		double w_max = max(max(G));
		for(sp_mat::iterator it = G.begin(); it != G.end(); ++ it) {
			epochs_per_sample[i] = w_max / (*it); // Higher the weight of the edge, the more likely it is to be sampled (inversely proportional to epochs_per_sample)
			positive_head[i] = it.row();
			positive_tail[i++] = it.col();
		}		
		
		// Initial coordinates of vertices (0-simplices)
		mat Ct = coors;
		vector<double> head_vec(Ct.memptr(), Ct.memptr()+Ct.n_elem);		
		vector<double> tail_vec(head_vec);		
		
		
		// Stores linearized coordinates [x1, y1, x2, y2, ...]
		vector<double> result;
		const apumap_gradient gradient(a_param, b_param, GAMMA);
		result = optimize_layout(
			gradient, head_vec, tail_vec, positive_head, positive_tail,
			n_epochs, N, epochs_per_sample, LEARNING_RATE,
			NEGATIVE_SAMPLE_RATE, UMAP_SEED);	

		mat coordinates(result.data(), 2, N);
				
		return coordinates;	
	}	

	
}





