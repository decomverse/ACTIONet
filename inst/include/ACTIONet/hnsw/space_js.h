#pragma once
#include <cmath>
#include "fastlog.h"
#include "hnswlib.h"

using std::numeric_limits;

namespace hnswlib {

static double JSD_metric(const void *pVect1_p, const void *pVect2_p,
                         const void *params) {
  double *log_vec = (double *)params;
  size_t N = (size_t)log_vec[0];

  double *pVect1 = (double *)pVect1_p;
  double *pVect2 = (double *)pVect2_p;
  double half = 0.5;

  double sum1 = 0, sum2 = 0;
  double res1 = 1.0, res2 = 1.0;
  for (size_t i = 0; i < N; i++) {
    double p = pVect1[i];
    double q = pVect2[i];

    res1 -= p;
    res2 -= q;

    double m = (p + q) * half;

    /*
    int p_idx = (int)floor(p *1000000.0);
    int q_idx = (int)floor(q *1000000.0);
    int m_idx = (int)floor(m *1000000.0);

    double lg_p = log_vec[p_idx];
    double lg_q = log_vec[q_idx];
    double lg_m = log_vec[m_idx];
    */

    double lg_p = fasterlog2(p);
    double lg_q = fasterlog2(q);
    double lg_m = fasterlog2(m);

    sum1 += (p * lg_p) + (q * lg_q);
    sum2 += m * lg_m;
  }

  double JS = std::max(half * sum1 - sum2, 0.0);

  JS += half * (std::max(res1, 0.0) +
                std::max(res2, 0.0));  // If there is residual, consider that as
                                       // as being in mismatched dimension

  return (double)sqrt(JS);
}

class JSDSpace : public SpaceInterface<double> {
  DISTFUNC<double> fstdistfunc_;
  size_t data_size_;
  double params[2];

 public:
  JSDSpace(size_t dim) {
    fstdistfunc_ = JSD_metric;
    data_size_ = dim * sizeof(double);
    params[0] = dim;
  }

  size_t get_data_size() { return data_size_; }

  DISTFUNC<double> get_dist_func() { return fstdistfunc_; }

  void *get_dist_func_param() { return (void *)params; }

  ~JSDSpace() {}
};

}  // namespace hnswlib
