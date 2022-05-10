#include <vector>
#include <ACTIONet.h>

#include "uwot/coords.h"
#include "uwot/epoch.h"
#include "uwot/gradient.h"
#include "uwot/optimize.h"
#include "uwot/sampler.h"

#include "rng.h"
#include "rparallel.h"


// Template class specialization to handle different rng/batch combinations
template <bool DoBatch = true> struct BatchRngFactory {
  using PcgFactoryType = batch_pcg_factory;
  using TauFactoryType = batch_tau_factory;
};
template <> struct BatchRngFactory<false> {
  using PcgFactoryType = pcg_factory;
  using TauFactoryType = tau_factory;
};

namespace ACTIONet {
struct UmapFactory {
  bool move_other;
  bool pcg_rand;
  std::vector<float> &head_embedding;
  std::vector<float> &tail_embedding;
  const std::vector<unsigned int> &positive_head;
  const std::vector<unsigned int> &positive_tail;
  const std::vector<unsigned int> &positive_ptr;
  unsigned int n_epochs;
  unsigned int n_head_vertices;
  unsigned int n_tail_vertices;
  const std::vector<float> &epochs_per_sample;
  float initial_alpha;
  float negative_sample_rate;
  bool batch;
  std::size_t n_threads;
  std::size_t grain_size;
  bool verbose;
  float adam_alpha, adam_beta1, adam_beta2, adam_eps, sgd_alpha;
  string opt_name = "adam";
  std::mt19937_64 engine;

  UmapFactory(bool move_other, bool pcg_rand,
              std::vector<float> &head_embedding,
              std::vector<float> &tail_embedding,
              const std::vector<unsigned int> &positive_head,
              const std::vector<unsigned int> &positive_tail,
              const std::vector<unsigned int> &positive_ptr,
              unsigned int n_epochs, unsigned int n_head_vertices,
              unsigned int n_tail_vertices,
              const std::vector<float> &epochs_per_sample, float initial_alpha,
              float negative_sample_rate, bool batch,
              std::size_t n_threads, std::size_t grain_size,
              bool verbose,
            float adam_alpha, float adam_beta1, float adam_beta2, float adam_eps,
            float sgd_alpha, string opt_name, std::mt19937_64 &engine)
      : move_other(move_other), pcg_rand(pcg_rand),
        head_embedding(head_embedding), tail_embedding(tail_embedding),
        positive_head(positive_head), positive_tail(positive_tail),
        positive_ptr(positive_ptr), n_epochs(n_epochs),
        n_head_vertices(n_head_vertices), n_tail_vertices(n_tail_vertices),
        epochs_per_sample(epochs_per_sample), initial_alpha(initial_alpha),
        negative_sample_rate(negative_sample_rate),
        batch(batch), n_threads(n_threads), grain_size(grain_size),
        verbose(verbose),
            adam_alpha(adam_alpha), adam_beta1(adam_beta1), adam_beta2(adam_beta2), adam_eps(adam_eps),
            sgd_alpha(sgd_alpha), opt_name(opt_name), engine(engine)    
         {}

  template <typename Gradient> void create(const Gradient &gradient) {
    if (move_other) {
      create_impl<true>(gradient, pcg_rand, batch);
    } else {
      create_impl<false>(gradient, pcg_rand, batch);
    }
  }

  template <bool DoMove, typename Gradient>
  void create_impl(const Gradient &gradient, bool pcg_rand, bool batch) {
    if (batch) {
      create_impl<BatchRngFactory<true>, DoMove>(gradient, pcg_rand, batch);
    } else {
      create_impl<BatchRngFactory<false>, DoMove>(gradient, pcg_rand, batch);
    }
  }

  template <typename BatchRngFactory, bool DoMove, typename Gradient>
  void create_impl(const Gradient &gradient, bool pcg_rand, bool batch) {
    if (pcg_rand) {
      create_impl<typename BatchRngFactory::PcgFactoryType, DoMove>(gradient,
                                                                    batch);
    } else {
      create_impl<typename BatchRngFactory::TauFactoryType, DoMove>(gradient,
                                                                    batch);
    }
  }

  auto create_adam(float adam_alpha, float adam_beta1, float adam_beta2, float adam_eps) -> uwot::Adam {
    return uwot::Adam(adam_alpha, adam_beta1, adam_beta2, adam_eps, head_embedding.size());
  }

  auto create_sgd(float sgd_alpha = 1.0) -> uwot::Sgd {
    return uwot::Sgd(sgd_alpha);
  }

  template <typename RandFactory, bool DoMove, typename Gradient>
  void create_impl(const Gradient &gradient, bool batch) {
    
    if (batch) {
      if (opt_name == "adam") {
        auto opt = create_adam(adam_alpha, adam_beta1, adam_beta2, adam_eps);
        create_impl_batch_opt<decltype(opt), RandFactory, DoMove, Gradient>(
            gradient, opt, batch);
      } else if (opt_name == "sgd") {
        auto opt = create_sgd(sgd_alpha);
        create_impl_batch_opt<decltype(opt), RandFactory, DoMove, Gradient>(
            gradient, opt, batch);
      } 
    } else {
      const std::size_t ndim = head_embedding.size() / n_head_vertices;
      uwot::Sampler sampler(epochs_per_sample, negative_sample_rate);
      uwot::InPlaceUpdate<DoMove> update(head_embedding, tail_embedding,
                                         initial_alpha);
      uwot::EdgeWorker<Gradient, decltype(update), RandFactory> worker(
          gradient, update, positive_head, positive_tail, sampler, ndim,
          n_tail_vertices, n_threads, engine);
      create_impl(worker, gradient);
    }
  }

  template <typename Opt, typename RandFactory, bool DoMove, typename Gradient>
  void create_impl_batch_opt(const Gradient &gradient, Opt &opt, bool batch) {
    uwot::Sampler sampler(epochs_per_sample, negative_sample_rate);
    const std::size_t ndim = head_embedding.size() / n_head_vertices;
    uwot::BatchUpdate<DoMove, decltype(opt)> update(
        head_embedding, tail_embedding, opt);

    uwot::NodeWorker<Gradient, decltype(update), RandFactory> worker(
        gradient, update, positive_head, positive_tail, positive_ptr, sampler,
        ndim, n_tail_vertices, engine);
    create_impl(worker, gradient);
  }

  template <typename Worker, typename Gradient>
  void create_impl(Worker &worker, const Gradient &gradient) {

    if (n_threads > 0) {
      RParallel parallel(n_threads, grain_size);
      create_impl(worker, gradient, parallel);
    } else {
      RSerial serial;
      create_impl(worker, gradient, serial);
    }
  }

  template <typename Worker, typename Gradient, typename Parallel>
  void create_impl(Worker &worker, const Gradient &gradient, Parallel &parallel) {
    uwot::optimize_layout(worker, n_epochs, parallel);
  }
};

void create_umap(UmapFactory &umap_factory, float a, float b, float gamma, bool approx_pow) {
  if (approx_pow) {
    const uwot::apumap_gradient gradient(a, b, gamma);
   umap_factory.create(gradient);
  } else {
    const uwot::umap_gradient gradient(a, b, gamma);
    umap_factory.create(gradient);
  }
}

void create_tumap(UmapFactory &umap_factory) {
  const uwot::tumap_gradient gradient;
  umap_factory.create(gradient);
}

void create_pacmap(UmapFactory &umap_factory, float a, float b) {
  const uwot::pacmap_gradient gradient(a, b);
  umap_factory.create(gradient);
}

void create_largevis(UmapFactory &umap_factory, float gamma) {
  const uwot::largevis_gradient gradient(gamma);
  umap_factory.create(gradient);
}

auto vec_to_coords(std::vector<float> &head_vec,
                 std::vector<float> &tail_vec) -> uwot::Coords {
  if (tail_vec.size() == 0) {
    return uwot::Coords(head_vec);
  } else {
    return uwot::Coords(head_vec, tail_vec);
  }
}

auto vec_to_coords(std::vector<float> &head_vec) -> uwot::Coords {
  return uwot::Coords(head_vec);
}

std::vector<float> optimize_layout_interface(
    std::vector<float> &head_embedding, std::vector<float> &tail_embedding,
    const std::vector<unsigned int> positive_head,
    const std::vector<unsigned int> positive_tail,
    const std::vector<unsigned int> positive_ptr, unsigned int n_epochs,
    unsigned int n_head_vertices, unsigned int n_tail_vertices,
    const std::vector<float> epochs_per_sample, const std::string &method,
    float initial_alpha, float a, float b, float gamma, bool approx_pow, float negative_sample_rate,
    bool pcg_rand = true, bool batch = true, std::size_t n_threads = 0,
    std::size_t grain_size = 1, bool move_other = true, bool verbose = false,
    float adam_alpha = 1.0, float adam_beta1 = 0.5, float adam_beta2 = 0.9, float adam_eps = 1e-7,
    float sgd_alpha = 1.0, string opt_name = "adam", int seed = 0) {

    std::mt19937_64 engine(seed);
  auto coords = vec_to_coords(head_embedding, tail_embedding);

  const std::size_t ndim = head_embedding.size() / n_head_vertices;

  UmapFactory umap_factory(move_other, pcg_rand, coords.get_head_embedding(),
                          coords.get_tail_embedding(), positive_head,
                          positive_tail, positive_ptr, n_epochs,
                          n_head_vertices, n_tail_vertices, epochs_per_sample,
                          initial_alpha, negative_sample_rate, batch,
                          n_threads, grain_size, verbose,
                          adam_alpha, adam_beta1, adam_beta2, adam_eps, sgd_alpha, opt_name, engine
                          );

  //printf("epochs = %d, head_embedding = %d, positive_head = %d, positive_ptr = %d\n", n_epochs, head_embedding.size(), positive_head.size(), positive_ptr.size());

  if (method == "umap") {
    create_umap(umap_factory, a, b, gamma, approx_pow);
  } else if (method == "tumap") {
    create_tumap(umap_factory);
  } else if (method == "largevis") {
    create_largevis(umap_factory, gamma);
  } else if (method == "pacmap") {
    create_pacmap(umap_factory, a, b);
  } else {
    fprintf(stderr, "Unknown layout method %s", method.c_str());
  }

  return(coords.get_head_embedding());
}
}