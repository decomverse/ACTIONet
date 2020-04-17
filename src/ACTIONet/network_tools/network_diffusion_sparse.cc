#include <ACTIONet.h>

namespace ACTIONet {	
	// Page 556
    sp_mat get_Q(sp_mat& A, vec& ds_inv, double alpha) {
       
        sp_mat A_norm = A;
        double c = (1-alpha) / 2.0;        
		for(sp_mat::iterator it = A_norm.begin(); it != A_norm.end(); ++it) {
			(*it) *= c*(ds_inv[it.row()]*ds_inv[it.col()]);
		}
        sp_mat I = speye(size(A));
		sp_mat Q = ((1 + alpha)/2.0)*I - A_norm;
		
		return(Q);
    }

	// Used in computing gradients
    sp_mat get_shift(sp_mat& S, vec& ds_inv, double alpha) {       		
        sp_mat shift = -alpha*S;
		for(sp_mat::iterator it = shift.begin(); it != shift.end(); ++it) {
			(*it) *= (ds_inv[it.row()]);
		}
		
		return( shift );
    }

	// Page 559
    sp_mat get_gradients(sp_mat& Q, sp_mat& X, sp_mat& shift) {
        return( (Q * X) + shift );
    }

	// Page 564
	sp_mat proximal_step(sp_mat& X, sp_mat& gradients, vec& kappa) {
        sp_mat Delta = -(X - gradients); 
        
        int counts = 0;
		for(sp_mat::iterator it = Delta.begin(); it != Delta.end(); ++it) {
			double delta = - ((*it) + kappa[it.row()]);			
			if( delta > 0 ) {
				counts ++;
				(*it) = delta;			
			}
		}
		
		Delta.transform( [](double val) { return (val < 0? 0:val); } );		
       		
		return(Delta);
    }

    sp_mat compute_sparse_network_diffusion(sp_mat& A, sp_mat& S, double alpha = 0.85, double rho = 1e-4, double epsilon = 0.001, int max_iter = 20) {
        
        alpha = 1.0 - alpha; // Make larger alpha mean deeper diffusion
        
        (A.diag()).zeros(); // remove self-loops
        
        vec d = vec(sum(A, 1));
        vec ds = sqrt(d);        
        vec ds_inv = ones(size(ds));
        uvec idx = find(ds > 0);        
		ds_inv(idx) = 1.0 / ds(idx);
        				
        sp_mat Q = get_Q(A, ds_inv, alpha);
        sp_mat shift = get_shift(S, ds_inv, alpha); // -alpha*D^(-1/2)*s
        
		sp_mat X(size(S)), next_X(size(S)); // X0 = 0
        sp_mat gradients = shift; // G0 = -alpha*D^(-1/2)*s
        
		vec kappa = alpha * rho * ds;
		
		printf("Running network diffusion for %d vectors\n", X.n_cols);
		for(int it = 1; it <= max_iter; it++) {
            next_X = proximal_step(X, gradients, kappa);            			
            gradients = get_gradients(Q, next_X, shift);            
            
            double delta_magnitude = norm(next_X - X, "fro");
            printf("\t%d- %.2e\n", it, delta_magnitude);
            
            if (delta_magnitude < epsilon) 
				break;
				
            X = next_X;
        }
        
        // Re-scale
		for(sp_mat::iterator it = X.begin(); it != X.end(); ++it) {
			(*it) *= ds[it.row()];
		}        
		
		return(X);		
    }    
}
