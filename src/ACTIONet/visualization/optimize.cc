// Adopted from the UWOT R package

#include <vector>

#include "perpendicular.h"
#include "gradient.h"
#include "optimize.h"
#include "sampler.h"


template <typename T, bool DoMove = true>
auto optimize_layout(const T &gradient, std::vector<double> &head_embedding,
                     std::vector<double> &tail_embedding,
                     const std::vector<unsigned int> &positive_head,
                     const std::vector<unsigned int> &positive_tail,
                     unsigned int n_epochs, unsigned int n_vertices,
                     const std::vector<double> &epochs_per_sample,
                     double initial_alpha, double negative_sample_rate,
                     std::size_t thread_no = 0, std::size_t grain_size = 1, unsigned int seed = 0) -> std::vector<double> {
  uwot::Sampler sampler(epochs_per_sample, negative_sample_rate);

  uwot::SgdWorker<T, DoMove> worker(
      gradient, positive_head, positive_tail, sampler, head_embedding,
      tail_embedding, head_embedding.size() / n_vertices, seed);

  const auto n_epochs_per_sample = epochs_per_sample.size();
  double alpha = initial_alpha;

	if (thread_no <= 0) {
		thread_no = std::thread::hardware_concurrency();
	}

		

  for (auto n = 0U; n < n_epochs; n++) {	  
    worker.set_alpha(alpha);
    worker.set_n(n);
    if (thread_no > 1) {
      Perpendicular::parallel_for(0, n_epochs_per_sample, worker, thread_no,
                                      grain_size);
    } else {
      worker(0, n_epochs_per_sample);
    }
    alpha = initial_alpha * (1.0 - (double(n) / double(n_epochs)));
  }
  
  return head_embedding;
}

std::vector<double>  optimize_layout_umap(
    std::vector<double> head_vec, std::vector<double> tail_vec,
    const std::vector<unsigned int> positive_head,
    const std::vector<unsigned int> positive_tail, unsigned int n_epochs,
    unsigned int n_vertices, const std::vector<double> epochs_per_sample,
    double a, double b, double gamma, double initial_alpha,
    double negative_sample_rate, bool approx_pow,
    std::size_t thread_no = 0, std::size_t grain_size = 1,
    bool move_other = true) {


  std::vector<double> result;
  if (approx_pow) {
    const uwot::apumap_gradient gradient(a, b, gamma);
    if (move_other) {
        result = optimize_layout<uwot::apumap_gradient, true>(
            gradient, head_vec, tail_vec, positive_head, positive_tail,
            n_epochs, n_vertices, epochs_per_sample, initial_alpha,
            negative_sample_rate, thread_no, grain_size);
    } else {
		result = optimize_layout<uwot::apumap_gradient, false>(
			gradient, head_vec, tail_vec, positive_head, positive_tail,
			n_epochs, n_vertices, epochs_per_sample, initial_alpha,
			negative_sample_rate, thread_no, grain_size);
    }
  } else {
    const uwot::umap_gradient gradient(a, b, gamma);
    if (move_other) {
        result = optimize_layout<uwot::umap_gradient, true>(
            gradient, head_vec, tail_vec, positive_head, positive_tail,
            n_epochs, n_vertices, epochs_per_sample, initial_alpha,
            negative_sample_rate, thread_no, grain_size);
    } else {
        result = optimize_layout<uwot::umap_gradient, false>(
            gradient, head_vec, tail_vec, positive_head, positive_tail,
            n_epochs, n_vertices, epochs_per_sample, initial_alpha,
            negative_sample_rate, thread_no, grain_size);      
    }
  }

  return(result);

}

std::vector<double> optimize_layout_tumap(
    std::vector<double> head_vec, std::vector<double> tail_vec,
    const std::vector<unsigned int> positive_head,
    const std::vector<unsigned int> positive_tail, unsigned int n_epochs,
    unsigned int n_vertices, const std::vector<double> epochs_per_sample,
    double initial_alpha, double negative_sample_rate,
    std::size_t thread_no = 0, std::size_t grain_size = 1,
    bool move_other = true) {
		
  std::vector<double> result;
  const uwot::tumap_gradient gradient;

  if (move_other) {
	  result = optimize_layout<uwot::tumap_gradient, true>(
		  gradient, head_vec, tail_vec, positive_head, positive_tail,
		  n_epochs, n_vertices, epochs_per_sample, initial_alpha,
		  negative_sample_rate, thread_no, grain_size);
  } else {
      result = optimize_layout<uwot::tumap_gradient, false>(
          gradient, head_vec, tail_vec, positive_head, positive_tail,
          n_epochs, n_vertices, epochs_per_sample, initial_alpha,
          negative_sample_rate, thread_no, grain_size);
  }

  return(result);

}

std::vector<double>  optimize_layout_largevis(
    std::vector<double> head_vec, const std::vector<unsigned int> positive_head,
    const std::vector<unsigned int> positive_tail, unsigned int n_epochs,
    unsigned int n_vertices, const std::vector<double> epochs_per_sample,
    double gamma, double initial_alpha, double negative_sample_rate,
    std::size_t thread_no = 0, std::size_t grain_size = 1) {

  const uwot::largevis_gradient gradient(gamma);

  std::vector<double> result = optimize_layout<uwot::largevis_gradient, true>(
        gradient, head_vec, head_vec, positive_head, positive_tail, n_epochs,
        n_vertices, epochs_per_sample, initial_alpha, negative_sample_rate,
        thread_no, grain_size);


  return(result);
}
