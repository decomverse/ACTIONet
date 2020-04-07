#include <sampler.h>

Sampler::Sampler(const std::vector<double>& epochs_per_sample, 
                 const double negative_sample_rate) :
  
  epochs_per_sample(epochs_per_sample),
  epoch_of_next_sample(epochs_per_sample),
  epochs_per_negative_sample(epochs_per_sample.size()),
  epoch_of_next_negative_sample(epochs_per_sample.size())
{
  const std::size_t esz = epochs_per_sample.size();
  const double nsr = 1.0 / negative_sample_rate;
  for (std::size_t i = 0; i < esz; i++) {
    epochs_per_negative_sample[i] = epochs_per_sample[i] * nsr;
    epoch_of_next_negative_sample[i] = epochs_per_negative_sample[i];
  }
}

bool Sampler::is_sample_edge(const std::size_t i, const std::size_t n) const {
  return epoch_of_next_sample[i] <= n;
}

const std::size_t Sampler::get_num_neg_samples(
    const std::size_t i, 
    const std::size_t n) const 
{
  return static_cast<std::size_t>(
    (n - epoch_of_next_negative_sample[i]) / epochs_per_negative_sample[i]);
}

void Sampler::next_sample(const std::size_t i, 
                          const std::size_t num_neg_samples) {
  epoch_of_next_sample[i] += epochs_per_sample[i];
  epoch_of_next_negative_sample[i] += 
    num_neg_samples * epochs_per_negative_sample[i];
}
