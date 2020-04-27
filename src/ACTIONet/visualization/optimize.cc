#include <limits>
#include <random>

#include <gradient.h>
#include <sampler.h>
#include <tauprng.h>
#include <math.h>

// Function to decide whether to move both vertices in an edge
// Default empty version does nothing: used in umap_transform when
// some of the vertices should be held fixed
template <bool DoMoveVertex = false>
void move_other_vertex(
    std::vector<double>& embedding,
    const double grad_d,
    const std::size_t i,
    const std::size_t nrj) {
}

// Specialization to move the vertex: used in umap when both
// vertices in an edge should be moved
template <>
void move_other_vertex<true>(
    std::vector<double>& embedding,
    const double grad_d,
    const std::size_t i,
    const std::size_t nrj) {
  embedding[nrj + i] -= grad_d;
}

const double clamp(const double v, const double lo, const double hi) {
  const double t = v < lo ? lo : v;
  return t > hi ? hi : t;
}

// Gradient: the type of gradient used in the optimization
// DoMoveVertex: true if both ends of a positive edge should be updated
template <typename Gradient,
          bool DoMoveVertex = true>
struct SgdWorker{
  int n; // epoch counter
  double alpha;
  const Gradient gradient;
  const std::vector<unsigned int> positive_head;
  const std::vector<unsigned int> positive_tail;
  Sampler sampler;
  std::vector<double>& head_embedding;
  std::vector<double>& tail_embedding;
  const std::size_t ndim;
  const std::size_t head_nvert;
  const std::size_t tail_nvert;
  std::mt19937 rng;
  std::uniform_int_distribution<long> gen;
  const double dist_eps;
  
  SgdWorker(
    const Gradient& gradient,
    const std::vector<unsigned int>& positive_head,
    const std::vector<unsigned int>& positive_tail,
    Sampler& sampler,
    std::vector<double>& head_embedding,
    std::vector<double>& tail_embedding,
    const std::size_t ndim,
    unsigned int seed) :
    
    n(0), alpha(0.0), gradient(gradient),
    positive_head(positive_head), positive_tail(positive_tail),
    
    sampler(sampler),
    
    head_embedding(head_embedding),
    tail_embedding(tail_embedding),
    ndim(ndim), 
    head_nvert(head_embedding.size() / ndim), 
    tail_nvert(tail_embedding.size() / ndim),
    rng(seed), gen(-2147483647, 2147483646),
    dist_eps(std::numeric_limits<double>::epsilon())
  { }
  
  void  operator()(std::size_t begin, std::size_t end) {
    // Each window gets its own fast PRNG state, so no locking needed inside the loop.
    // Want separate seeds though, so seed the fast PRNG with three random numbers
    // taken from the mt19937 generator, which is shared across windows, so is locked.
    // Probably this is a bit of a waste of time:
    // Could use the mt19937 seed, the begin and the epoch number as seeds?
    // Doesn't waste much time, though.
    long s1, s2, s3;
    {
      s1 = gen(rng);
      s2 = gen(rng); // technically this needs to always be > 7
      s3 = gen(rng); // should be > 15
    }
    tau_prng prng(s1, s2, s3);
    
    std::vector<double> dys(ndim);
    for (std::size_t i = begin; i < end; i++) {
      if (!sampler.is_sample_edge(i, n)) {
        continue;
      }
      
      const std::size_t dj = ndim * positive_head[i];
      const std::size_t dk = ndim * positive_tail[i];

      double dist_squared = 0.0;
      for (std::size_t d = 0; d < ndim; d++) {
        const double diff = head_embedding[dj + d] - tail_embedding[dk + d];
        dys[d] = diff;
        dist_squared += diff * diff;
      }
      dist_squared = std::max(dist_eps, dist_squared);
      
      const double grad_coeff = gradient.grad_attr(dist_squared);
      for (std::size_t d = 0; d < ndim; d++) {
        const double grad_d = alpha * clamp(grad_coeff * dys[d], 
                                            Gradient::clamp_lo,
                                            Gradient::clamp_hi);
        head_embedding[dj + d] += grad_d;
		
        move_other_vertex<DoMoveVertex>(tail_embedding, grad_d, d, dk);
      }

      const std::size_t n_neg_samples = sampler.get_num_neg_samples(i, n);
      for (std::size_t p = 0; p < n_neg_samples; p++) {
        const std::size_t dkn = (prng() % tail_nvert) * ndim;
        if (dj == dkn) {
          continue;
        }

        double dist_squared = 0.0;
        for (std::size_t d = 0; d < ndim; d++) {
          const double diff = head_embedding[dj + d] - tail_embedding[dkn + d];
          dys[d] = diff;
          dist_squared += diff * diff;
        }
        dist_squared = std::max(dist_eps, dist_squared);
        
        const double grad_coeff = gradient.grad_rep(dist_squared);
        for (std::size_t d = 0; d < ndim; d++) {
          const double grad_d = alpha * clamp(grad_coeff * dys[d], 
                                              Gradient::clamp_lo,
                                              Gradient::clamp_hi);
          head_embedding[dj + d] += grad_d;
        }
      }
      
      sampler.next_sample(i, n_neg_samples);
    }
  }
  
  void set_n(int n) {
    this->n = n;
  }
  
  void set_alpha(double alpha) {
    this->alpha = alpha;
  }
};


std::vector<double> optimize_layout(
    const apumap_gradient& gradient,
    std::vector<double>& head_embedding,
    std::vector<double>& tail_embedding,
    const std::vector<unsigned int>& positive_head,
    const std::vector<unsigned int>& positive_tail,
    unsigned int n_epochs, 
    unsigned int n_vertices,
    const std::vector<double>& epochs_per_sample,
    double initial_alpha,
    double negative_sample_rate,
    unsigned int seed) 
{
	/*
	gradient.print();
	
	for(int e = 0; e < 10; e++) {
		printf("e%d: %d %d %.2f\n", e, positive_head[e], positive_tail[e], epochs_per_sample[e]);
	}

	
	for(int v = 0; v < 10; v++) {
		printf("v%d: %.2f %.2f\n", v+1, head_embedding[2*v], head_embedding[2*v+1]);
	}
	*/
	
  Sampler sampler(epochs_per_sample, negative_sample_rate);
  
  SgdWorker<apumap_gradient, true> worker(gradient, 
                              positive_head, positive_tail,
                              sampler,
                              head_embedding, tail_embedding,
                              head_embedding.size() / n_vertices,
                              seed);
  
  
  const auto n_epochs_per_sample = epochs_per_sample.size();
  double alpha = initial_alpha;
  
  for (auto n = 0U; n < n_epochs; n++) {
    worker.set_alpha(alpha);
    worker.set_n(n);
    
	worker(0, n_epochs_per_sample);
    alpha = initial_alpha * (1.0 - (double(n) / double(n_epochs)));
    
  }

  
  return head_embedding;
}
