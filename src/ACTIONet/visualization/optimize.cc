// Adopted from the UWOT R package

#include <vector>

#include "perpendicular.h"
#include "gradient.h"
#include "optimize.h"
#include "sampler.h"


template <typename T, bool DoMove = true>
auto optimize_layout(const T &gradient, std::vector<float> &head_embedding,
                     std::vector<float> &tail_embedding,
                     const std::vector<unsigned int> &positive_head,
                     const std::vector<unsigned int> &positive_tail,
                     unsigned int n_epochs, unsigned int n_vertices,
                     const std::vector<float> &epochs_per_sample,
                     float initial_alpha, float negative_sample_rate,
                     std::size_t n_threads = 0, std::size_t grain_size = 1, unsigned int seed = 0) -> std::vector<float> {
  uwot::Sampler sampler(epochs_per_sample, negative_sample_rate);

  uwot::SgdWorker<T, DoMove> worker(
      gradient, positive_head, positive_tail, sampler, head_embedding,
      tail_embedding, head_embedding.size() / n_vertices, seed);

  const auto n_epochs_per_sample = epochs_per_sample.size();
  float alpha = initial_alpha;

  for (auto n = 0U; n < n_epochs; n++) {
    worker.set_alpha(alpha);
    worker.set_n(n);
    if (n_threads > 0) {
      Perpendicular::parallel_for(0, n_epochs_per_sample, worker, n_threads,
                                      grain_size);
    } else {
      worker(0, n_epochs_per_sample);
    }
    alpha = initial_alpha * (1.0 - (float(n) / float(n_epochs)));
  }
  return head_embedding;
}

std::vector<float>  optimize_layout_umap(
    std::vector<float> head_vec, std::vector<float> tail_vec,
    const std::vector<unsigned int> positive_head,
    const std::vector<unsigned int> positive_tail, unsigned int n_epochs,
    unsigned int n_vertices, const std::vector<float> epochs_per_sample,
    float a, float b, float gamma, float initial_alpha,
    float negative_sample_rate, bool approx_pow,
    std::size_t n_threads = 0, std::size_t grain_size = 1,
    bool move_other = true) {


  std::vector<float> result;
  if (approx_pow) {
    const uwot::apumap_gradient gradient(a, b, gamma);
    if (move_other) {
        result = optimize_layout<uwot::apumap_gradient, true>(
            gradient, head_vec, tail_vec, positive_head, positive_tail,
            n_epochs, n_vertices, epochs_per_sample, initial_alpha,
            negative_sample_rate, n_threads, grain_size);
    } else {
		result = optimize_layout<uwot::apumap_gradient, false>(
			gradient, head_vec, tail_vec, positive_head, positive_tail,
			n_epochs, n_vertices, epochs_per_sample, initial_alpha,
			negative_sample_rate, n_threads, grain_size);
    }
  } else {
    const uwot::umap_gradient gradient(a, b, gamma);
    if (move_other) {
        result = optimize_layout<uwot::umap_gradient, true>(
            gradient, head_vec, tail_vec, positive_head, positive_tail,
            n_epochs, n_vertices, epochs_per_sample, initial_alpha,
            negative_sample_rate, n_threads, grain_size);
    } else {
        result = optimize_layout<uwot::umap_gradient, false>(
            gradient, head_vec, tail_vec, positive_head, positive_tail,
            n_epochs, n_vertices, epochs_per_sample, initial_alpha,
            negative_sample_rate, n_threads, grain_size);      
    }
  }

  return(result);

}

std::vector<float> optimize_layout_tumap(
    std::vector<float> head_vec, std::vector<float> tail_vec,
    const std::vector<unsigned int> positive_head,
    const std::vector<unsigned int> positive_tail, unsigned int n_epochs,
    unsigned int n_vertices, const std::vector<float> epochs_per_sample,
    float initial_alpha, float negative_sample_rate,
    std::size_t n_threads = 0, std::size_t grain_size = 1,
    bool move_other = true) {
		
  std::vector<float> result;
  const uwot::tumap_gradient gradient;

  if (move_other) {
	  result = optimize_layout<uwot::tumap_gradient, true>(
		  gradient, head_vec, tail_vec, positive_head, positive_tail,
		  n_epochs, n_vertices, epochs_per_sample, initial_alpha,
		  negative_sample_rate, n_threads, grain_size);
  } else {
      result = optimize_layout<uwot::tumap_gradient, false>(
          gradient, head_vec, tail_vec, positive_head, positive_tail,
          n_epochs, n_vertices, epochs_per_sample, initial_alpha,
          negative_sample_rate, n_threads, grain_size);
  }

  return(result);

}

std::vector<float>  optimize_layout_largevis(
    std::vector<float> head_vec, const std::vector<unsigned int> positive_head,
    const std::vector<unsigned int> positive_tail, unsigned int n_epochs,
    unsigned int n_vertices, const std::vector<float> epochs_per_sample,
    float gamma, float initial_alpha, float negative_sample_rate,
    std::size_t n_threads = 0, std::size_t grain_size = 1) {

  const uwot::largevis_gradient gradient(gamma);

  std::vector<float> result = optimize_layout<uwot::largevis_gradient, true>(
        gradient, head_vec, head_vec, positive_head, positive_tail, n_epochs,
        n_vertices, epochs_per_sample, initial_alpha, negative_sample_rate,
        n_threads, grain_size);


  return(result);
}
