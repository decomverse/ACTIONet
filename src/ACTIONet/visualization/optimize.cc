// Adopted from the UWOT R package
#include <ACTIONet.h>
#include <vector>

#include "gradient.h"
#include "optimize.h"
#include "perpendicular.h"
#include "sampler.h"

template <typename T, bool DoMove = true>
auto optimize_layout(const T &gradient, std::vector<float> &head_embedding,
                     std::vector<float> &tail_embedding,
                     const std::vector<unsigned int> &positive_head,
                     const std::vector<unsigned int> &positive_tail,
                     unsigned int n_epochs, unsigned int n_vertices,
                     const std::vector<float> &epochs_per_sample,
                     double initial_alpha, double negative_sample_rate,
                     std::size_t thread_no = 0, std::size_t grain_size = 1,
                     unsigned int seed = 0) -> std::vector<float> {
  uwot::Sampler sampler(epochs_per_sample, negative_sample_rate);

  uwot::SgdWorker<T, DoMove> worker(gradient, positive_head, positive_tail,
                                    sampler, head_embedding, tail_embedding,
                                    head_embedding.size() / n_vertices, seed);

  const auto n_epochs_per_sample = epochs_per_sample.size();
  double alpha = initial_alpha;

  if (thread_no <= 0) {
    thread_no = SYS_THREADS_DEF;  // std::thread::hardware_concurrency();
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

std::vector<float> optimize_layout_umap(
    std::vector<float> head_vec, std::vector<float> tail_vec,
    const std::vector<unsigned int> positive_head,
    const std::vector<unsigned int> positive_tail, unsigned int n_epochs,
    unsigned int n_vertices, const std::vector<float> epochs_per_sample,
    double a, double b, double gamma, double initial_alpha,
    double negative_sample_rate, bool approx_pow, std::size_t thread_no = 0,
    std::size_t grain_size = 1, bool move_other = true, int seed = 0) {
  std::vector<float> result;
  if (approx_pow) {
    const uwot::apumap_gradient gradient(a, b, gamma);
    if (move_other) {
      result = optimize_layout<uwot::apumap_gradient, true>(
          gradient, head_vec, tail_vec, positive_head, positive_tail, n_epochs,
          n_vertices, epochs_per_sample, initial_alpha, negative_sample_rate,
          thread_no, grain_size, seed);
    } else {
      result = optimize_layout<uwot::apumap_gradient, false>(
          gradient, head_vec, tail_vec, positive_head, positive_tail, n_epochs,
          n_vertices, epochs_per_sample, initial_alpha, negative_sample_rate,
          thread_no, grain_size, seed);
    }
  } else {
    const uwot::umap_gradient gradient(a, b, gamma);
    if (move_other) {
      result = optimize_layout<uwot::umap_gradient, true>(
          gradient, head_vec, tail_vec, positive_head, positive_tail, n_epochs,
          n_vertices, epochs_per_sample, initial_alpha, negative_sample_rate,
          thread_no, grain_size, seed);
    } else {
      result = optimize_layout<uwot::umap_gradient, false>(
          gradient, head_vec, tail_vec, positive_head, positive_tail, n_epochs,
          n_vertices, epochs_per_sample, initial_alpha, negative_sample_rate,
          thread_no, grain_size, seed);
    }
  }

  return (result);
}

std::vector<float> optimize_layout_tumap(
    std::vector<float> head_vec, std::vector<float> tail_vec,
    const std::vector<unsigned int> positive_head,
    const std::vector<unsigned int> positive_tail, unsigned int n_epochs,
    unsigned int n_vertices, const std::vector<float> epochs_per_sample,
    double initial_alpha, double negative_sample_rate,
    std::size_t thread_no = 0, std::size_t grain_size = 1,
    bool move_other = true, int seed = 0) {
  std::vector<float> result;
  const uwot::tumap_gradient gradient;

  if (move_other) {
    result = optimize_layout<uwot::tumap_gradient, true>(
        gradient, head_vec, tail_vec, positive_head, positive_tail, n_epochs,
        n_vertices, epochs_per_sample, initial_alpha, negative_sample_rate,
        thread_no, grain_size, seed);
  } else {
    result = optimize_layout<uwot::tumap_gradient, false>(
        gradient, head_vec, tail_vec, positive_head, positive_tail, n_epochs,
        n_vertices, epochs_per_sample, initial_alpha, negative_sample_rate,
        thread_no, grain_size, seed);
  }

  return (result);
}

std::vector<float> optimize_layout_largevis(
    std::vector<float> head_vec, const std::vector<unsigned int> positive_head,
    const std::vector<unsigned int> positive_tail, unsigned int n_epochs,
    unsigned int n_vertices, const std::vector<float> epochs_per_sample,
    double gamma, double initial_alpha, double negative_sample_rate,
    std::size_t thread_no = 0, std::size_t grain_size = 1, int seed = 0) {
  const uwot::largevis_gradient gradient(gamma);

  std::vector<float> result = optimize_layout<uwot::largevis_gradient, true>(
      gradient, head_vec, head_vec, positive_head, positive_tail, n_epochs,
      n_vertices, epochs_per_sample, initial_alpha, negative_sample_rate,
      thread_no, grain_size, seed);

  return (result);
}
 // BSD 2-Clause License
//
// Copyright 2020 James Melville
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// OF SUCH DAMAGE.

#ifndef UWOT_OPTIMIZE_H
#define UWOT_OPTIMIZE_H

#include <limits>
#include <utility>

#include <pcg_random.hpp>

#include "sampler.h"
#include "tauprng.h"

namespace uwot {

inline auto clamp(float v, float lo, float hi) -> float {
  float t = v < lo ? lo : v;
  return t > hi ? hi : t;
}

// Gradient: the type of gradient used in the optimization
// DoMoveVertex: true if both ends of a positive edge should be updated
template <typename Gradient, bool DoMoveVertex>
struct SgdWorker {
  int n;  // epoch counter
  float alpha;
  const Gradient gradient;
  const std::vector<unsigned int> positive_head;
  const std::vector<unsigned int> positive_tail;
  uwot::Sampler sampler;
  std::vector<float> &head_embedding;
  std::vector<float> &tail_embedding;
  std::size_t ndim;
  std::size_t head_nvert;
  std::size_t tail_nvert;
  float dist_eps;

  thread_local pcg32 rng;

  SgdWorker(const Gradient &gradient, std::vector<unsigned int> positive_head,
            std::vector<unsigned int> positive_tail, uwot::Sampler &sampler,
            std::vector<float> &head_embedding,
            std::vector<float> &tail_embedding, std::size_t ndim,
            uint64_t seed = std::mt19937_64::default_seed)
      :

        n(0),
        alpha(0.0),
        gradient(gradient),
        positive_head(positive_head),
        positive_tail(positive_tail),

        sampler(sampler),

        head_embedding(head_embedding),
        tail_embedding(tail_embedding),
        ndim(ndim),
        head_nvert(head_embedding.size() / ndim),
        tail_nvert(tail_embedding.size() / ndim),
        dist_eps(std::numeric_limits<float>::epsilon()) {
    rng.seed(seed);
  }

  void operator()(std::size_t begin, std::size_t end) {
    std::vector<float> dys(ndim);
	
    std::uniform_int_distribution<int> uniform_dist(0, tail_nvert - 1);

	long long ss = 0, tt = 0;
	double g1 = 0, g2 = 0;
	
    for (auto i = begin; i < end; i++) {
      if (!sampler.is_sample_edge(i, n)) {
        continue;
      }
		tt += i;
		
      std::size_t dj = ndim * positive_head[i];
      std::size_t dk = ndim * positive_tail[i];

      float dist_squared = 0.0;
      for (std::size_t d = 0; d < ndim; d++) {
        float diff = head_embedding[dj + d] - tail_embedding[dk + d];
        dys[d] = diff;
        dist_squared += diff * diff;
      }
      dist_squared = (std::max)(dist_eps, dist_squared);

      float grad_coeff = gradient.grad_attr(dist_squared);
		g1 += grad_coeff;
      for (std::size_t d = 0; d < ndim; d++) {
        float grad_d = alpha * clamp(grad_coeff * dys[d], Gradient::clamp_lo,
                                     Gradient::clamp_hi);
        head_embedding[dj + d] += grad_d;
		if(DoMoveVertex)
			tail_embedding[dk + d] -= grad_d;

        //move_other_vertex<DoMoveVertex>(tail_embedding, grad_d, d, dk);
      }

      std::size_t n_neg_samples = sampler.get_num_neg_samples(i, n);
      for (std::size_t p = 0; p < n_neg_samples; p++) {
        int r = uniform_dist(rng);
        ss += r;
        std::size_t dkn = r * ndim;
        if (dj == dkn) {
          continue;
        }
        float dist_squared = 0.0;
        for (std::size_t d = 0; d < ndim; d++) {
          float diff = head_embedding[dj + d] - tail_embedding[dkn + d];
          dys[d] = diff;
          dist_squared += diff * diff;
        }
        dist_squared = (std::max)(dist_eps, dist_squared);

        float grad_coeff = gradient.grad_rep(dist_squared);
		g2 += grad_coeff;

        for (std::size_t d = 0; d < ndim; d++) {
          float grad_d = alpha * clamp(grad_coeff * dys[d], Gradient::clamp_lo,
                                       Gradient::clamp_hi);

          head_embedding[dj + d] += grad_d;
        }
      }
      sampler.next_sample(i, n_neg_samples);
    }
    
    printf("ss = %ld, tt = %ld, g1 = %e, g2 = %e\n", ss, tt, g1, g2);
  }

  void set_n(int n) { this->n = n; }

  void set_alpha(float alpha) { this->alpha = alpha; }

  void reseed(int new_seed) { rng.seed(new_seed); }
};
}  // namespace uwot

#endif  // UWOT_OPTIMIZE_H
      
