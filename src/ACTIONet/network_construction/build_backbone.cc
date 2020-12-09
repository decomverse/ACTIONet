#include <ACTIONet.h>

#include <hnswlib.h>
#include <thread>
#include <atomic>

template<class Function>
inline void ParallelFor(size_t start, size_t end, size_t numThreads, Function fn) {
    if (numThreads <= 0) {
        numThreads = SYS_THREADS_DEF;
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

namespace ACTIONet {

	double threshold_vector(vec &dist, double LC) {
		vec dist_sorted = sort(dist, "ascend");
		vec beta = join_vert(LC*dist_sorted, datum::inf*ones(1));
		vec beta_sq = square(beta);

		double lambda = beta(0)+1, k = 0, Sum_beta = 0, Sum_beta_square = 0;

		for(; k < beta.n_elem-1; k++) {
			Sum_beta += beta[k];
			Sum_beta_square += (beta_sq[k]);
			lambda = (1.0 / k) * ( Sum_beta + sqrt( k  + Sum_beta*Sum_beta - k * Sum_beta_square ) ) ;

			if(lambda <= beta(k+1))
				break;
		}

		double dist_threshold = dist_sorted[(int)k];

		return(dist_threshold);
	}


	// k^{*}-Nearest Neighbors: From Global to Local (NIPS 2016)
	mat Prune_PageRank(mat &U, double density = 1.0) {

		double LC = 1.0 / density;

		mat Dist = -log(U);

		mat U_thresholded(size(U));
		for(int i = 0; i < U.n_cols; i++) {

			vec d = Dist.col(i);
			double dist_threshold = threshold_vector(d, LC);

			uvec filter_idx = find(dist_threshold <= Dist.col(i));
			vec u = U.col(i);
			u(filter_idx).zeros();

			U_thresholded.col(i) = u;
		}

		U_thresholded = normalise(U_thresholded, 1, 0);
		/*
		for(int i = 0; i < U.n_cols; i++) {
			vec u = U.col(i);

			vec u_sorted = sort(u, "descend");
			vec beta = join_cols(-LC*log(u_sorted), datum::inf*ones(1));

			vec beta_cs = cumsum(beta);
			vec beta_sq_cs = cumsum(square(beta));

			vec kk = regspace(1, beta_cs.n_elem);

			vec lambda = (1.0/kk) % (beta_cs + sqrt(kk + square(beta_cs) - kk % beta_sq_cs));
			lambda.replace(datum::nan, 0);

			uvec tmp = find(lambda <= beta);
			if(tmp.n_elem == 0) {
				printf("Can't prune column %d with density = %f\n", i+1, density);
				U_thresholded.col(i) = u;
				continue;
			}
			int k = tmp(0);
			double u_threshold = u_sorted(k);

			printf("%d- k = %d, threshold = %f\n", i+1, k+1, u_threshold);

			u(u < u_threshold).zeros();

			uvec idx = find(u >= u_threshold);
			vec x = u(idx);
			vec z = zscore(x);
			vec s = exp(z);
			vec p = s / sum(s);
			u(idx) = p;

			U_thresholded.col(i) = u;
		}
		*/

		/*

		mat U_sorted = sort(U, "descend", 0);

		mat beta = trans(-LC*log(U_sorted));
		vec beta_sum = zeros(beta.n_rows);
		vec beta_sq_sum = zeros(beta.n_rows);

		mat lambda = zeros(size(beta));
		int k = 1;
		for(; k < beta.n_cols; k++) {
			beta_sum += beta.col(k);
			beta_sq_sum += square(beta.col(k));

			lambda.col(k) = (1.0/(double)k) * ( beta_sum + sqrt(k + square(beta_sum) - k*beta_sq_sum) );
		}
		lambda.replace(datum::nan, 0);

		mat Delta = trans(lambda - beta);
		Delta.shed_row(0);


		mat U_thresholded(size(U));
		for(int v = 0; v < U.n_cols; v++) {
		//ParallelFor(0, U.n_cols, 1, [&](size_t v, size_t threadId) {
			vec delta = Delta.col(v);
			vec u = U_thresholded.col(v);

			//uvec rows = find(delta > 0, 1, "last");
			uvec tmp = find(delta < 0, 1, "first");
			printf("%d- %d\n", v, tmp.n_elem);
			if(tmp.n_elem == 0) {
				printf("Can't prune column %d with density = %f\n", v+1, density);
				U_thresholded.col(v) = u;
			} else {
				int k = tmp(0);
				double u_threshold = U_sorted(k);
				printf("\t%f\n", u_threshold);

				u(u < u_threshold).zeros();

				uvec idx = find(u >= u_threshold);
				vec x = u(idx);
				vec s = x;//exp(x); // exp(zscore(x));
				vec p = s / sum(s);
				u(idx) = p;

				U_thresholded.col(v) = u;
			}
		}
		*/

		return(U_thresholded);
	}
}
