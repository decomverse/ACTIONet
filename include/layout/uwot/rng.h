#ifndef UWOT_RNG_H
#define UWOT_RNG_H

#include <limits>

// linked from dqrng
#include "convert_seed.h"
#include "pcg_random.hpp"

#include "tauprng.h"

#include <math/StatsLib/stats.hpp>

// NOT THREAD SAFE
// based on code in the dqsample package
static uint64_t random64(std::mt19937_64 &engine)
{
  return static_cast<uint64_t>(
      stats::runif(0, (std::numeric_limits<uint64_t>::max)(), engine));
}

// NOT THREAD SAFE
static uint32_t random32(std::mt19937_64 &engine)
{
  return static_cast<uint32_t>(
      stats::runif(0, (std::numeric_limits<uint32_t>::max)(), engine));
}

struct batch_tau_factory
{
  std::size_t n_rngs;
  std::vector<uint64_t> seeds;
  static const constexpr std::size_t seeds_per_rng = 3;

  batch_tau_factory() : n_rngs(1), seeds(seeds_per_rng * n_rngs) {}
  batch_tau_factory(std::size_t n_rngs)
      : n_rngs(n_rngs), seeds(seeds_per_rng * n_rngs) {}

  void reseed(std::mt19937_64 engine)
  {
    for (std::size_t i = 0; i < seeds.size(); i++)
    {
      seeds[i] = random64(engine);
    }
  }

  uwot::tau_prng create(std::size_t n)
  {
    const std::size_t idx = n * seeds_per_rng;
    return uwot::tau_prng(seeds[idx], seeds[idx + 1], seeds[idx + 2]);
  }

  uwot::tau_prng create2(std::size_t n1, std::size_t n2)
  {
    return uwot::tau_prng(seeds[n1 * seeds_per_rng], seeds[n1 * seeds_per_rng + 1], seeds[n2 * seeds_per_rng]);
  }

  void print(std::size_t n1, std::size_t n2)
  {
  }
};

struct pcg_prng
{
  pcg32 gen;

  pcg_prng(uint64_t seed) { gen.seed(seed); }

  // return a value in (0, n]
  inline std::size_t operator()(std::size_t n)
  {
    std::size_t result = gen(n);
    return result;
  }
};

struct batch_pcg_factory
{
  std::size_t n_rngs;

  std::vector<uint32_t> seeds;
  static const constexpr std::size_t seeds_per_rng = 2;

  batch_pcg_factory() : n_rngs(1), seeds(seeds_per_rng * n_rngs) {}
  batch_pcg_factory(std::size_t n_rngs)
      : n_rngs(n_rngs), seeds(seeds_per_rng * n_rngs) {}

  void reseed(std::mt19937_64 engine)
  {
    for (std::size_t i = 0; i < seeds.size(); i++)
    {
      seeds[i] = random32(engine);
    }
  }

  pcg_prng create(std::size_t n)
  {
    uint32_t pcg_seeds[2] = {seeds[n * seeds_per_rng],
                             seeds[n * seeds_per_rng + 1]};
    // printf("%d- <%d, %d>\n", n, pcg_seeds[0], pcg_seeds[1]);

    return pcg_prng(dqrng::convert_seed<uint64_t>(pcg_seeds, 2));
  }

  pcg_prng create2(std::size_t n1, std::size_t n2)
  {
    uint32_t pcg_seeds[2] = {seeds[n1 * seeds_per_rng],
                             seeds[n2 * seeds_per_rng]};
    // printf("%d- <%d, %d>\n", n, pcg_seeds[0], pcg_seeds[1]);

    return pcg_prng(dqrng::convert_seed<uint64_t>(pcg_seeds, 2));
  }

  void print(std::size_t n1, std::size_t n2)
  {
    uint32_t pcg_seeds[2] = {seeds[n1 * seeds_per_rng],
                             seeds[n2 * seeds_per_rng]};

    printf("<%d, %d> -> <%d, %d> -> %d\n", n1, n2, pcg_seeds[0], pcg_seeds[1], dqrng::convert_seed<uint64_t>(pcg_seeds, 2));
  }
};

// For backwards compatibility in non-batch mode
struct tau_factory
{
  uint64_t seed1;
  uint64_t seed2;

  tau_factory(std::size_t) : seed1(0), seed2(0)
  {
  }

  void reseed(std::mt19937_64 engine)
  {
    seed1 = random64(engine);
    seed2 = random64(engine);
  }

  uwot::tau_prng create(std::size_t seed)
  {
    return uwot::tau_prng(seed1, seed2, uint64_t{seed});
  }

  uwot::tau_prng create2(std::size_t seed2, std::size_t seed3)
  {
    return uwot::tau_prng(seed1, uint64_t{seed2}, uint64_t{seed3});
  }

  void print(std::size_t n1, std::size_t n2)
  {
  }
};

struct pcg_factory
{
  uint32_t seed1;

  pcg_factory(std::size_t) : seed1(0) {}

  void reseed(std::mt19937_64 engine) { seed1 = random32(engine); }

  pcg_prng create(std::size_t seed)
  {
    uint32_t seeds[2] = {seed1, static_cast<uint32_t>(seed)};
    return pcg_prng(dqrng::convert_seed<uint64_t>(seeds, 2));
  }

  pcg_prng create2(std::size_t seed1, std::size_t seed2)
  {
    uint32_t seeds[2] = {static_cast<uint32_t>(seed1), static_cast<uint32_t>(seed2)};
    return pcg_prng(dqrng::convert_seed<uint64_t>(seeds, 2));
  }

  void print(std::size_t n1, std::size_t n2)
  {
  }
};

#endif // UWOT_RNG_H
