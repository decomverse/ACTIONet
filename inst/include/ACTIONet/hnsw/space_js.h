#pragma once
#include "hnswlib.h"
#include "fastlog.h"
#include <cmath>

using std::numeric_limits;

namespace hnswlib {

    static double JSD_metric(const void *pVect1_p, const void *pVect2_p, const void *params) {
		double *log_vec = (double *)params;
        size_t N = (size_t) log_vec[1000001];
        
		double * pVect1 = (double *)pVect1_p;
		double * pVect2 = (double *)pVect2_p;
		double half = 0.5;
		
		double sum1 = 0, sum2 = 0;
		for (size_t i = 0; i < N; i++) {
			double p = pVect1[i];
			double q = pVect2[i];
			double m = (p + q)*half;

			int p_idx = (int)floor(p *1000000.0);
			int q_idx = (int)floor(q *1000000.0);
			int m_idx = (int)floor(m *1000000.0);

			double lg_p = log_vec[p_idx];
			double lg_q = log_vec[q_idx];
			double lg_m = log_vec[m_idx];

			sum1 += (p * lg_p) + (q * lg_q);
			sum2 +=  m * lg_m;		  
		}
        
        double JS = std::max(half*sum1 - sum2, 0.0);
        
		return (double)sqrt(JS);
	}
	
	
    class JSDSpace : public SpaceInterface<double> {

        DISTFUNC<double> fstdistfunc_;
        size_t data_size_;
        double params[1000002];
        
    public:
        JSDSpace(size_t dim) {
            fstdistfunc_ = JSD_metric;
            data_size_ = dim * sizeof(double);
            
            for(register int i = 0; i <= 1000000; i++) {
				params[i] = (double)log2((double)i / 1000000.0);
			}
			params[1000001] = dim;
			params[0] = 0;			
        }

        size_t get_data_size() {
            return data_size_;
        }

        DISTFUNC<double> get_dist_func() {
            return fstdistfunc_;
        }

        void *get_dist_func_param() {
            return (void *)params;
        }

        ~JSDSpace() {}
    };

}
