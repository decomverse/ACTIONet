#include <ACTIONet.h>

#include <thread>
#include <atomic>


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

namespace ACTIONet {	
	arma::vec diffusion_solve_FISTA(arma::sp_mat& adj_mat, arma::vec& prob_dist, double alpha, double rho, double epsilon, int max_iter);

	mat PR_linsys(sp_mat &G, sp_mat &X, double alpha = 0.85, int thread_no = -1) {		
		X = normalise(X, 1, 0);
		
		/*
		rowvec d = sum(G, 0);
		uvec idx = find(c != 0);
		d[idx] = 1 / d[idx];
		
		sp_mat D;		
		D.diag() = d;
		
		sp_mat I = speye(size(G));
		*/
		sp_mat P = normalise(G, 1, 0);		
		sp_mat I = speye(P.n_rows, P.n_cols);
		sp_mat A = I - alpha*P;		
		//mat PR = (1-alpha)*spsolve(A, mat(X), "superlu"); 
		mat PR = (1-alpha)*spsolve(A, mat(X)); 

		
		return(PR);
	}


	mat compute_network_diffusion(sp_mat &G, sp_mat &X0, int thread_no = 4, double alpha = 0.85, int max_it = 3) {		
		thread_no = std::min(thread_no, (int)X0.n_cols);
		
		int N = G.n_rows;
		vec z = ones(N);
		vec c = vec(trans(sum(G, 0)));
		uvec idx = find(c);
		z(idx) = ones(idx.n_elem)*(1.0 - alpha);
		z = z / N;
		
		sp_mat P = alpha*normalise(G, 1, 0);
		X0 = normalise(X0, 1, 0);
		mat X = mat(X0);
				
		X0 *= N;			
		rowvec zt = trans(z);
		ParallelFor(0, X.n_cols, thread_no, [&](size_t i, size_t threadId) {			
			X.col(i) = P*X.col(i) + X0.col(i)*(zt*X.col(i));
		});
		//X = normalise(X, 1)
		
		return(X);
	}	
	
	
}





