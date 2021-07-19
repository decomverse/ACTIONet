#include <ACTIONet.h>

namespace ACTIONet {
// Page 556
sp_mat get_Q(sp_mat& A, vec& ds_inv, double alpha) {
  sp_mat A_norm = A;
  double c = (1 - alpha) / 2.0;
  for (sp_mat::iterator it = A_norm.begin(); it != A_norm.end(); ++it) {
    (*it) *= c * (ds_inv[it.row()] * ds_inv[it.col()]);
  }
  sp_mat I = speye(size(A));
  sp_mat Q = ((1.0 + alpha) / 2.0) * I - A_norm;

  return (Q);
}

// Used in computing gradients
sp_mat get_shift(sp_mat& U, vec& ds_inv, double alpha) {
  sp_mat shift = -alpha * U;
  for (sp_mat::iterator it = shift.begin(); it != shift.end(); ++it) {
    (*it) *= (ds_inv[it.row()]);
  }

  return (shift);
}

sp_mat proximal_step(sp_mat& Y, sp_mat& gradients, vec& kappa) {
  sp_mat Delta = Y - gradients;
  sp_mat Z = abs(Delta);
  for (sp_mat::iterator it = Z.begin(); it != Z.end(); ++it) {
    (*it) -= (kappa[it.row()]);
  }
  Z.transform([](double val) { return (val < 0 ? 0 : val); });

  return (arma::sign(Delta) % Z);
}

sp_mat compute_sparse_network_diffusion(sp_mat& A, sp_mat& U,
                                        double alpha = 0.85, double rho = 1e-4,
                                        double epsilon = 0.001,
                                        int max_iter = 20) {
  alpha = 1.0 - alpha;  // Make larger alpha mean deeper diffusion
  double beta = (1 - sqrt(alpha)) / (1 + sqrt(alpha));

  (A.diag()).zeros();  // remove self-loops

  vec d = vec(sum(A, 1));
  vec ds = sqrt(d);
  vec ds_inv = ones(size(ds));
  uvec idx = find(ds > 0);
  ds_inv(idx) = 1.0 / ds(idx);

  sp_mat Q = get_Q(
      A, ds_inv,
      alpha);  // .5 * ((1 + alpha)*I - (1 - alpha)*(D^(-1/2) * A * D^(-1/2)));
  sp_mat Shift = get_shift(U, ds_inv, alpha);  // [-alpha*D^(-1/2)*s_j]

  vec kappa = alpha * rho * ds;
  sp_mat Y(size(U)), YY(size(U));
  printf("Running network diffusion for %d vectors (%d iterations)\n", U.n_cols,
         max_iter);
  for (int it = 1; it <= max_iter; it++) {
    sp_mat gradients = (Q * Y) + Shift;
    sp_mat next_YY = proximal_step(Y, gradients, kappa);

    Y = next_YY + beta * (next_YY - YY);

    double delta_magnitude = norm(next_YY - YY, "fro");
    printf("\t%d- %.2e\n", it, delta_magnitude);
    if (delta_magnitude < epsilon) break;

    YY = next_YY;
  }

  for (sp_mat::iterator it = YY.begin(); it != YY.end(); ++it) {
    (*it) *= ds[it.row()];
  }
  return YY;
}
}  // namespace ACTIONet
