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

#ifndef UWOT_EPOCH_H
#define UWOT_EPOCH_H

#include "sampler.h"
#include "update.h"

namespace uwot
{
  template <typename Worker, typename Parallel>
  void optimize_layout(Worker &worker, unsigned int n_epochs,
                       Parallel &parallel)
  {
    for (auto n = 0U; n < n_epochs; n++)
    {
      run_epoch(worker, n, n_epochs, parallel);
    }
  }

  template <typename Worker, typename Parallel>
  void run_epoch(Worker &worker, std::size_t epoch, std::size_t n_epochs,
                 Parallel &parallel)
  {
    worker.epoch_begin(epoch, n_epochs);
    parallel.pfor(worker.n_items, worker);
    worker.epoch_end(epoch, n_epochs, parallel);
  }

  // Gradient: the type of gradient used in the optimization
  // Update: type of update to the embedding coordinates
  template <typename Gradient, typename Update, typename RngFactory>
  struct EdgeWorker
  {
    const Gradient gradient;
    Update &update;
    const std::vector<unsigned int> &positive_head;
    const std::vector<unsigned int> &positive_tail;
    uwot::Sampler sampler;
    std::size_t ndim;
    std::size_t n_tail_vertices;
    std::size_t n_items;
    std::size_t n_threads;
    RngFactory rng_factory;
    std::mt19937_64 engine;

    std::vector<float> old_head_embedding;
    std::vector<float> old_tail_embedding;

    EdgeWorker(const Gradient &gradient, Update &update,
               const std::vector<unsigned int> &positive_head,
               const std::vector<unsigned int> &positive_tail,
               uwot::Sampler &sampler, std::size_t ndim,
               std::size_t n_tail_vertices, std::size_t n_threads, std::mt19937_64 &engine)
        : gradient(gradient), update(update), positive_head(positive_head),
          positive_tail(positive_tail), sampler(sampler), ndim(ndim),
          n_tail_vertices(n_tail_vertices), n_items(positive_head.size()),
          n_threads(std::max(n_threads, std::size_t{1})),
          rng_factory(this->n_threads), engine(engine) {}

    void epoch_begin(std::size_t epoch, std::size_t n_epochs)
    {
      rng_factory.reseed(engine);
      sampler.epoch = epoch;
      update.epoch_begin(epoch, n_epochs);
    }

    template <typename Parallel>
    void epoch_end(std::size_t epoch, std::size_t n_epochs, Parallel &parallel)
    {
      update.epoch_end(epoch, n_epochs, parallel);
    }

    void operator()(std::size_t begin, std::size_t end, std::size_t thread_id)
    {
      // Each window gets its own PRNG state, to prevent locking inside the loop.
      auto prng = rng_factory.create(thread_id);

      // displacement between two points, cost of reallocating inside the loop
      // is noticeable, also cheaper to calculate it once in the d2 calc
      std::vector<float> disp(ndim);
      for (auto edge = begin; edge < end; edge++)
      {
        process_edge(update, gradient, sampler, prng, positive_head,
                     positive_tail, ndim, n_tail_vertices, edge, thread_id, disp, old_head_embedding, old_tail_embedding);
      }
    }
  };

  template <typename Gradient, typename Update, typename RngFactory>
  struct NodeWorker
  {
    const Gradient gradient;
    Update &update;

    std::vector<float> old_head_embedding;
    std::vector<float> old_tail_embedding;

    const std::vector<unsigned int> &positive_head;
    const std::vector<unsigned int> &positive_tail;
    const std::vector<unsigned int> &positive_ptr;
    uwot::Sampler sampler;
    std::size_t ndim;
    std::size_t n_tail_vertices;
    std::size_t n_items;
    RngFactory rng_factory;
    std::mt19937_64 engine;

    NodeWorker(const Gradient &gradient, Update &update,
               const std::vector<unsigned int> &positive_head,
               const std::vector<unsigned int> &positive_tail,
               const std::vector<unsigned int> &positive_ptr,
               uwot::Sampler &sampler, std::size_t ndim,
               std::size_t n_tail_vertices, std::mt19937_64 &engine)
        : gradient(gradient), update(update), positive_head(positive_head),
          positive_tail(positive_tail), positive_ptr(positive_ptr),
          sampler(sampler), ndim(ndim), n_tail_vertices(n_tail_vertices),
          n_items(positive_ptr.size() - 1), rng_factory(n_tail_vertices), engine(engine) {}

    void epoch_begin(std::size_t epoch, std::size_t n_epochs)
    {
      rng_factory.reseed(engine);
      sampler.epoch = epoch;
      update.epoch_begin(epoch, n_epochs);

      old_head_embedding.assign(update.head_embedding.begin(), update.head_embedding.end());
      old_tail_embedding.assign(update.tail_embedding.begin(), update.tail_embedding.end());
    }

    template <typename Parallel>
    void epoch_end(std::size_t epoch, std::size_t n_epochs, Parallel &parallel)
    {
      update.epoch_end(epoch, n_epochs, parallel);

      for (int i = 0; i < update.head_embedding.size(); i++)
      {
        update.head_embedding[i] = (update.head_embedding[i] + old_head_embedding[i]) / 2.0;
        update.tail_embedding[i] = (update.tail_embedding[i] + old_tail_embedding[i]) / 2.0;
      }
    }

    void operator()(std::size_t begin, std::size_t end, std::size_t thread_id)
    {
      std::vector<float> disp(ndim);
      for (auto p = begin; p < end; p++)
      {
        for (auto edge = positive_ptr[p]; edge < positive_ptr[p + 1]; edge++)
        {
          if (edge < 10)
          {
            // rng_factory.print(positive_head[edge], positive_tail[edge]);
            // auto prng = rng_factory.create2(positive_head[edge], positive_tail[edge]);
            // printf("%d => %d\n", edge, prng(n_tail_vertices));
          }

          auto prng = rng_factory.create2(positive_head[edge], positive_tail[edge]);
          process_edge(update, gradient, sampler, prng, positive_head,
                       positive_tail, ndim, n_tail_vertices, edge, thread_id,
                       disp, old_head_embedding, old_tail_embedding);
        }
      }
    }
  };
} // namespace uwot

#endif // UWOT_EPOCH_H