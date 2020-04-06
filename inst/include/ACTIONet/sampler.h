#ifndef UWOT_SAMPLER_H
#define UWOT_SAMPLER_H

#include <vector>

// Weighted edge sampler
class Sampler {
public:
  Sampler(const std::vector<double>& epochs_per_sample, 
          const double negative_sample_rate);
  
  bool is_sample_edge(const std::size_t i, const std::size_t n) const;
  const std::size_t get_num_neg_samples(const std::size_t i, const std::size_t n) const;
  void next_sample(const std::size_t i, const std::size_t num_neg_samples);
  
private:
  std::vector<double> epochs_per_sample;
  std::vector<double> epoch_of_next_sample;
  std::vector<double> epochs_per_negative_sample;
  std::vector<double> epoch_of_next_negative_sample;
};

#endif // UWOT_SAMPLER_H
