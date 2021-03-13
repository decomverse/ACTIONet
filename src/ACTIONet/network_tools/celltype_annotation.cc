#include <ACTIONet.h>

double r8_normal_01_cdf_inverse(double p);

namespace ACTIONet {
template <class Function>
inline void ParallelFor(size_t start, size_t end, size_t thread_no,
                        Function fn) {
  if (thread_no <= 0) {
    thread_no = std::thread::hardware_concurrency();
  }

  if (thread_no == 1) {
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

    for (size_t threadId = 0; threadId < thread_no; ++threadId) {
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

mat RIN_transform(mat A, int thread_no = 4) {
  int M = A.n_rows;
  int N = A.n_cols;

  mat Zr = zeros(M, N);
  ParallelFor(0, N, thread_no, [&](size_t i, size_t threadId) {
    vec v = A.col(i);

    uvec row_perm_forward = stable_sort_index(v);
    uvec row_perm = stable_sort_index(row_perm_forward);
    vec p = (row_perm + ones(size(row_perm))) / (row_perm.n_elem + 1);

    vec v_RINT = zeros(size(p));
    for (int j = 0; j < p.n_elem; j++) {
      double norm_inv = r8_normal_01_cdf_inverse(p(j));
      v_RINT(j) = norm_inv;
    }

    Zr.col(i) = v_RINT;
  });

  return (Zr);
}

mat compute_marker_aggregate_stats(sp_mat &G, sp_mat &S, sp_mat &marker_mat,
                                   double alpha = 0.85, int max_it = 5,
                                   int thread_no = 0, bool ignore_baseline_expression = false) {
  mat stats = zeros(S.n_cols, marker_mat.n_cols);

  int n = G.n_rows;
  sp_mat o = sp_mat(ones(n, 1));
  // vec pr = compute_network_diffusion(G, o, thread_no, alpha, max_it).col(0);
  vec pr =
      compute_network_diffusion_fast(G, o, thread_no, alpha, max_it).col(0);


  for (int i = 0; i < marker_mat.n_cols; i++) {
	int marker_count = (int)sum(sum(spones(marker_mat.col(i))));
	
    int idx = 0;
    vec w = zeros(marker_count);
    vec baseline = zeros(marker_count);    
    sp_mat raw_expression(S.n_cols, marker_count);
    for (sp_mat::col_iterator it = marker_mat.begin_col(i);
         it != marker_mat.end_col(i); it++) {
      raw_expression.col(idx) = trans(S.row(it.row()));      
      w(idx) = (*it);
	  baseline(idx) = accu(raw_expression.col(idx));
      idx++;
    }
    if(!ignore_baseline_expression) {
		baseline = baseline / sum(baseline); //sqrt(sum(square(baseline)));
		w = w % baseline;
	}
    w = w / sqrt(sum(square(w)));    

    // mat imputed_expression = compute_network_diffusion(G, raw_expression,
    // thread_no, alpha, max_it);
    mat imputed_expression = compute_network_diffusion_fast(
        G, raw_expression, thread_no, alpha, max_it);
        
    
    for(int j = 0; j < imputed_expression.n_cols; j++) {
		vec ppr = imputed_expression.col(j);
		vec scores = log2(ppr / pr);
		uvec zero_idx = find(ppr == 0);
		scores(zero_idx).zeros();
		scores = scores % ppr;
/*
		// Rank-based inverse normal transformation
		uvec row_perm_forward = stable_sort_index(scores);
		uvec row_perm = stable_sort_index(row_perm_forward);
		vec p = (row_perm + ones(size(row_perm))) / (row_perm.n_elem + 1);
		vec z = zeros(size(p));
		for (int j = 0; j < p.n_elem; j++) {
		  double norm_inv = r8_normal_01_cdf_inverse(p(j));
		  z(j) = norm_inv;
		}		
		*/
		stats.col(i) += w(j) * scores;
	}
  }

  return (stats);
}

}  // namespace ACTIONet
