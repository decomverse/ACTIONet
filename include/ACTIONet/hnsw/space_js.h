#pragma once
#include <cmath>
#include "fastlog.h"
#include "hnswlib.h"

#define LOGLEN 1000000

using std::numeric_limits;

namespace hnswlib {

static float JSD_metric(const void *pVect1_p, const void *pVect2_p,
                         const void *params) {
  float *log_vec = (float *)params;
  size_t N = (size_t) log_vec[LOGLEN+1];
  

  float *pVect1 = (float *)pVect1_p;
  float *pVect2 = (float *)pVect2_p;

  float sum1 = 0, sum2 = 0;
  float res1 = 1.0, res2 = 1.0;
  for (size_t i = 0; i < N; i++) {
    float p = pVect1[i];
    float q = pVect2[i];

    res1 -= p;
    res2 -= q;

    float m = 0.5 * (p + q);

	
    int p_idx = (int)floor(p *LOGLEN);
    int q_idx = (int)floor(q *LOGLEN);
    int m_idx = (int)floor(m *LOGLEN);

    float lg_p = log_vec[p_idx];
    float lg_q = log_vec[q_idx];
    float lg_m = log_vec[m_idx];

/*
    float lg_p = fasterlog2(p);
    float lg_q = fasterlog2(q);
    float lg_m = fasterlog2(m);
*/    
    float plgp = (p == 0? 0:p*lg_p);
    float qlgq = (q == 0? 0:q*lg_q);
    float mlgm = (m == 0? 0:m*lg_m);

	//printf("%d- %f %f %f -> %f %f %f\n", i, p, q, m, plgp, qlgq, mlgm);
    
    sum1 += (plgp) + (qlgq);
    sum2 += mlgm;
  }

  float JS = 0.5*sum1 - sum2;  
  JS = 0 < JS? JS:0;  
  JS += 0.5 * ( (0 < res1? res1:0) + (0 < res2? res2:0) );
  
  return (float)sqrt(JS);
}

class JSDSpace : public SpaceInterface<float> {
  DISTFUNC<float> fstdistfunc_;
  size_t data_size_;
  float params[LOGLEN+2];

 public:
  JSDSpace(size_t dim) {
    fstdistfunc_ = JSD_metric;
    data_size_ = dim * sizeof(float);
        
	params[0] = 0;	
	for( int i = 1; i <= LOGLEN; i++) {
		params[i] = (float)fasterlog2((float)i / LOGLEN);
	}
	params[LOGLEN+1] = dim;
  }

  size_t get_data_size() { return data_size_; }

  DISTFUNC<float> get_dist_func() { return fstdistfunc_; }

  void *get_dist_func_param() { return (void *)params; }

  ~JSDSpace() {}
};

}  // namespace hnswlib
