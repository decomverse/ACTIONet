#include <ACTIONet.h>
//#include <mini_cblas.h>
#include <cassert>

// Re-implemented from: Fast and Robust Archetypal Analysis for Representation Learning
namespace ACTIONet {
	double pstart = 1.0, pmax = 1e14, pincfactor = 1.5, pinctol = 1e-7, pstoptol = 1e-7;	
	float cg_tol = 1e-6; int cg_maxIter = 1000;
	
	vec SimplexProjection(vec &y, double r = 1) {
		int n = y.n_elem;		

		double s = 0 , lambda = 0;
		vec yproj = sort(y, "descend");		
		for(int i = 0; i < n; i++) {
			s += yproj[i];
			lambda = (s - r) / (i+1);
			if ( (i < (n-1)) && ( (lambda < yproj[i]) && (yproj[i + 1] <= lambda) ) ) {
				break;
			}
		}
		
		for(int i = 0; i < n; i++) {
			yproj[i] = std::max(y[i] - lambda, 0.0);
		}
		
		return(yproj);
	}
	
	mat run_simplex_regression_proxdist(mat &X, mat &Y, int pmaxiter = 100, int pincmaxiter = 200) {
		// Precompute eigen-decomposition of X'X
		vec s;
		mat V;
		eig_sym( s, V, trans(X)*X );		
		mat Vt = trans(V);
		
		// Precompute X'y vectors
		mat XtY = trans(X)*Y; 

		
		// pre-allocate working arrays	
		vec beta(X.n_cols), beta_proj(X.n_cols), beta_outer(X.n_cols), beta_inner(X.n_cols);				
		mat Beta = zeros(X.n_cols, Y.n_cols);		
		
		beta = Beta.col(0);
		vec beta_proj0 = SimplexProjection(beta);
		for(int c = 0; c < Y.n_cols; c++) {
			beta_proj = beta_proj0;
			vec Xty = XtY.col(c); // p x 1	
		
			double p = pstart;		
			for (int iter = 1; iter <= pmaxiter; iter++) {
				beta_outer = beta; // record solution at previous p
				for (int iter_inner = 1; iter_inner <= pincmaxiter; iter_inner++) {
					beta_inner = beta;
					
					// solve β = (X'X + ρ * I) \ (X'y + ρ * βproj)
					beta = V * ( (Vt * (Xty + p * beta_proj)) / (s + p) );
					
					// project β
					beta_proj = SimplexProjection(beta);				
					
					// decide to increase p or not
					if ( norm(beta - beta_inner) < (pinctol * (norm(beta_inner) + 1)) ) {
						break;
					}
				}
				
				// decide to stop or not
				double d1 = norm(beta - beta_outer);
				double d2 = norm(beta - beta_proj);
				if ( (d1 < pstoptol * (norm(beta_outer) + 1)) && 
					(d2 < pstoptol * (norm(beta_proj) + 1)) ) {
					break;
				} else {
					p = std::min(pincfactor * p, pmax);
				}
			}
			
			Beta.col(c) = SimplexProjection(beta);		
		}
		
		return(Beta);
	}
	
	
}

