#include "ACTIONet.h"

namespace ACTIONet {
	mat sampleUnif(int l, int m, double a, double b, int seed) {
		std::default_random_engine gen (seed);	
		std::uniform_real_distribution<double> unif(a, b);
		
		
		mat R(l, m);
		for (register int j = 0; j < m; j++) {
			for(register int i = 0; i < l; i++) {
				R(i, j) = unif(gen);
			}
		}
		return R;
	}

	void gram_schmidt(mat& A) {
		for(uword i = 0; i < A.n_cols; ++i) {
			for(uword j = 0; j < i; ++j) {
				double r = dot(A.col(i), A.col(j));
				A.col(i) -= r * A.col(j);
			}
			
			double col_norm = norm(A.col(i), 2);
			
			if(col_norm < 1E-4) {
				for(uword k = i; k < A.n_cols; ++k)
					A.col(k).zeros();

				return;
			}
			A.col(i) /= col_norm;
		}
	}
	
	field<mat> eigSVD(mat A) {
		int n = A.n_cols;
		mat B = trans(A)*A;
				
		vec d;
		mat V;
		eig_sym( d, V, B );		
		d = sqrt(d);
		
		// Compute U
		sp_mat S(n, n);
		S.diag() = 1 / d; 
		mat U = (S*trans(V))*trans(A);
		U = trans(U);
		
		field<mat> out(3);
		
		out(0) = U;
		out(1) = d;
		out(2) = V;
		
		return(out);
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


	
}
