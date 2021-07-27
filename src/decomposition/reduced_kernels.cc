#include "ACTIONet.h"

namespace ACTIONet
{
  field<mat> perturbedSVD(field<mat> SVD_results, mat &A, mat &B)
  {
    // printf("\tComputing perturbed SVD ... "); fflush(stdout);
    stdout_printf("\tComputing perturbed SVD ... "); // fflush(stdout);
    int n = A.n_rows;

    mat U = SVD_results(0);
    vec s = SVD_results(1);
    mat V = SVD_results(2);

    int dim = U.n_cols;

    vec s_prime;
    mat U_prime, V_prime;

    mat M = U.t() * A;
    mat A_ortho_proj = A - U * M;
    mat P = A_ortho_proj; // = orth(A_ortho_proj);
    gram_schmidt(P);
    mat R_P = P.t() * A_ortho_proj;

    mat N = V.t() * B;
    mat B_ortho_proj = B - V * N;
    mat Q = B_ortho_proj; // orth(B_ortho_proj);
    gram_schmidt(Q);
    mat R_Q = Q.t() * B_ortho_proj;

    mat K1 = zeros(s.n_elem + A.n_cols, s.n_elem + A.n_cols);
    for (int i = 0; i < s.n_elem; i++)
    {
      K1(i, i) = s(i);
    }

    mat K2 = join_vert(M, R_P) * trans(join_vert(N, R_Q));

    mat K = K1 + K2;

    svd(U_prime, s_prime, V_prime, K);

    mat U_updated = join_horiz(U, P) * U_prime;
    mat V_updated = join_horiz(V, Q) * V_prime;
    stdout_printf("done\n"); // fflush(stdout);

    field<mat> output(5);
    output(0) = U_updated.cols(0, dim - 1);
    output(1) = s_prime(span(0, dim - 1));
    output(2) = V_updated.cols(0, dim - 1);

    if ((SVD_results.n_elem == 5) && (SVD_results(3).n_elem != 0))
    {
      output(3) = join_rows(SVD_results(3), A);
      output(4) = join_rows(SVD_results(4), B);
    }
    else
    {
      output(3) = A;
      output(4) = B;
    }

    return output;
  }

  field<mat> PCA2SVD(sp_mat &S, field<mat> PCA_results)
  {
    int n = S.n_rows;

    stdout_printf("PCA => SVD (sparse)\n"); // fflush(stdout);
    mat U = PCA_results(0);
    vec s = PCA_results(1);
    mat V = PCA_results(2);

    int dim = U.n_cols;

    mat A = ones(S.n_rows, 1);
    mat B = mat(trans(mean(S, 0)));

    field<mat> perturbed_SVD = perturbedSVD(PCA_results, A, B);

    return perturbed_SVD;
  }

  field<mat> PCA2SVD(mat &S, field<mat> PCA_results)
  {
    int n = S.n_rows;

    stdout_printf("PCA => SVD (dense)\n"); // fflush(stdout);
    mat U = PCA_results(0);
    vec s = PCA_results(1);
    mat V = PCA_results(2);

    int dim = U.n_cols;

    mat A = ones(S.n_rows, 1);
    mat B = mat(trans(mean(S, 0)));

    field<mat> perturbed_SVD = perturbedSVD(PCA_results, A, B);

    return perturbed_SVD;
  }

  field<mat> SVD2PCA(sp_mat &S, field<mat> SVD_results)
  {
    int n = S.n_rows;

    stdout_printf("SVD => PCA (sparse)\n"); // fflush(stdout);
    mat U = SVD_results(0);
    vec s = SVD_results(1);
    mat V = SVD_results(2);

    int dim = U.n_cols;

    mat A = ones(S.n_rows, 1);
    mat B = -mat(trans(mean(S, 0)));

    field<mat> perturbed_SVD = perturbedSVD(SVD_results, A, B);

    return perturbed_SVD;
  }

  field<mat> SVD2PCA(mat &S, field<mat> SVD_results)
  {
    int n = S.n_rows;

    stdout_printf("SVD => PCA (dense)\n"); // fflush(stdout);
    mat U = SVD_results(0);
    vec s = SVD_results(1);
    mat V = SVD_results(2);

    int dim = U.n_cols;

    mat A = ones(S.n_rows, 1);
    mat B = -mat(trans(mean(S, 0)));

    field<mat> perturbed_SVD = perturbedSVD(SVD_results, A, B);

    return perturbed_SVD;
  }

  field<mat> SVD2ACTIONred(sp_mat &S, field<mat> SVD_results)
  {
    stdout_printf("SVD => ACTIONred (sparse)\n"); // fflush(stdout);
    int n = S.n_rows;
    int dim = SVD_results(0).n_cols;

    // Update 1: Orthogonalize columns w.r.t. background (mean)
    vec mu = vec(mean(S, 1));
    vec v = mu / norm(mu, 2);
    vec a1 = v;
    vec b1 = -trans(S) * v;

    // Update 2: Center columns of orthogonalized matrix before performing SVD
    vec c = vec(trans(mean(S, 0)));
    double a1_mean = mean(a1);
    vec a2 = ones(S.n_rows);
    vec b2 = -(a1_mean * b1 + c);

    mat A = join_rows(a1, a2);
    mat B = join_rows(b1, b2);

    field<mat> perturbed_SVD = perturbedSVD(SVD_results, A, B);

    return perturbed_SVD;
  }

  field<mat> SVD2ACTIONred(mat &S, field<mat> SVD_results)
  {
    stdout_printf("SVD => ACTIONred (dense)\n"); // fflush(stdout);
    int n = S.n_rows;
    int dim = SVD_results(0).n_cols;

    // Update 1: Orthogonalize columns w.r.t. background (mean)
    vec mu = vec(mean(S, 1));
    vec v = mu / norm(mu, 2);
    vec a1 = v;
    vec b1 = -trans(S) * v;

    // Update 2: Center columns of orthogonalized matrix before performing SVD
    vec c = vec(trans(mean(S, 0)));
    double a1_mean = mean(a1);
    vec a2 = ones(S.n_rows);
    vec b2 = -(a1_mean * b1 + c);

    mat A = join_rows(a1, a2);
    mat B = join_rows(b1, b2);

    field<mat> perturbed_SVD = perturbedSVD(SVD_results, A, B);

    return perturbed_SVD;
  }

  field<mat> PCA2ACTIONred(sp_mat &S, field<mat> PCA_results)
  {
    stdout_printf("Reverting column-centering ... "); // fflush(stdout);
    field<mat> SVD_results = PCA2SVD(S, PCA_results);
    stdout_printf("done\n");

    field<mat> output = SVD2ACTIONred(S, SVD_results);
    return output;
  }

  field<mat> PCA2ACTIONred(mat &S, field<mat> PCA_results)
  {
    stdout_printf("Reverting column-centering ... "); // fflush(stdout);
    field<mat> SVD_results = PCA2SVD(S, PCA_results);
    stdout_printf("done\n");

    field<mat> output = SVD2ACTIONred(S, SVD_results);
    return output;
  }

  field<mat> reduce_kernel(sp_mat &S, int dim, int iter = 5, int seed = 0,
                           int SVD_algorithm = HALKO_ALG,
                           bool prenormalize = false,
                           int verbose = 1)
  {
    int n = S.n_rows;
    if (prenormalize)
      S = normalise(S, 2);

    stdout_printf(
        "Computing reduced ACTION kernel (sparse):\n"); // fflush(stdout);

    stdout_printf("\tPerforming SVD on original matrix: "); // fflush(stdout);
    vec s;
    mat U, V;
    field<mat> SVD_results;
    switch (SVD_algorithm)
    {
    case FULL_SVD:
      svd_econ(U, s, V, mat(S));
      SVD_results(0) = U;
      SVD_results(1) = s;
      SVD_results(2) = V;
      break;
    case IRLB_ALG:
      SVD_results = IRLB_SVD(S, dim, iter, seed, verbose);
      break;
    case HALKO_ALG:
      SVD_results = HalkoSVD(S, dim, iter, seed, verbose);
      break;
    case FENG_ALG:
      SVD_results = FengSVD(S, dim, iter, seed, verbose);
      break;
    default:
      // fprintf(stderr, "Unknown SVD algorithm chosen (%d). Switching to
      // Halko.\n", SVD_algorithm);
      stderr_printf("Unknown SVD algorithm chosen (%d). Switching to Halko.\n",
                    SVD_algorithm);
      SVD_results = HalkoSVD(S, dim, iter, seed, verbose);
      break;
    }

    // Update 1: Orthogonalize columns w.r.t. background (mean)
    vec mu = vec(mean(S, 1));
    vec v = mu / norm(mu, 2);
    vec a1 = v;
    vec b1 = -trans(S) * v;

    // Update 2: Center columns of orthogonalized matrix before performing SVD
    vec c = vec(trans(mean(S, 0)));
    double a1_mean = mean(a1);
    vec a2 = ones(S.n_rows);
    vec b2 = -(a1_mean * b1 + c);

    mat A = join_rows(a1, a2);
    mat B = join_rows(b1, b2);

    field<mat> perturbed_SVD = perturbedSVD(SVD_results, A, B);

    return perturbed_SVD;
  }

  field<mat> reduce_kernel(mat &S, int dim, int iter = 5, int seed = 0,
                           int SVD_algorithm = HALKO_ALG,
                           bool prenormalize = false,
                           int verbose = 1)
  {
    int n = S.n_rows;
    if (prenormalize)
      S = normalise(S, 2);

    stdout_printf(
        "Computing reduced ACTION kernel (dense):\n"); // fflush(stdout);

    stdout_printf("\tPerforming SVD on original matrix: "); // fflush(stdout);
    vec s;
    mat U, V;
    field<mat> SVD_results;
    switch (SVD_algorithm)
    {
    case FULL_SVD:
      svd_econ(U, s, V, S);
      SVD_results(0) = U;
      SVD_results(1) = s;
      SVD_results(2) = V;
      break;
    case IRLB_ALG:
      SVD_results = IRLB_SVD(S, dim, iter, seed, verbose);
      break;
    case HALKO_ALG:
      SVD_results = HalkoSVD(S, dim, iter, seed, verbose);
      break;
    case FENG_ALG:
      SVD_results = FengSVD(S, dim, iter, seed, verbose);
      break;
    default:
      // fprintf(stderr, "Unknown SVD algorithm chosen (%d). Switching to
      // Halko.\n", SVD_algorithm);
      stderr_printf("Unknown SVD algorithm chosen (%d). Switching to Halko.\n",
                    SVD_algorithm);
      SVD_results = HalkoSVD(S, dim, iter, seed, verbose);
      break;
    }

    // Update 1: Orthogonalize columns w.r.t. background (mean)
    vec mu = vec(mean(S, 1));
    vec v = mu / norm(mu, 2);
    vec a1 = v;
    vec b1 = -trans(S) * v;

    // Update 2: Center columns of orthogonalized matrix before performing SVD
    vec c = vec(trans(mean(S, 0)));
    double a1_mean = mean(a1);
    vec a2 = ones(S.n_rows);
    vec b2 = -(a1_mean * b1 + c);

    mat A = join_rows(a1, a2);
    mat B = join_rows(b1, b2);

    field<mat> perturbed_SVD = perturbedSVD(SVD_results, A, B);

    return perturbed_SVD;
  }

  field<mat> ACTIONred2SVD(field<mat> SVD_results)
  {
    stdout_printf("ACTION kernel => SVD\n"); // fflush(stdout);

    mat A = -1 * SVD_results(3); // Reverting
    mat B = SVD_results(4);

    field<mat> perturbed_SVD = perturbedSVD(SVD_results, A, B);

    return perturbed_SVD;
  }

  field<mat> deflate_reduction(field<mat> SVD_results, mat &A, mat &B)
  {
    stdout_printf("\tDeflating reduction ... "); // fflush(stdout);

    vec mu_A = vec(trans(mean(A, 0)));
    vec mu = B * mu_A;

    A = join_rows(ones(A.n_rows), A);
    B = join_rows(-mu, B);
    stdout_printf("done\n");
    FLUSH;

    field<mat> perturbed_SVD = perturbedSVD(SVD_results, A, B);
    return (perturbed_SVD);
  }

  field<mat> orthogonalize_batch_effect(sp_mat &S, field<mat> SVD_results,
                                        mat &design)
  {
    stdout_printf("Orthogonalizing batch effect (sparse):\n"); // fflush(stdout);

    mat Z = mat(S * design);
    gram_schmidt(Z);

    mat A = Z;
    mat B = -mat(trans(trans(Z) * S));

    field<mat> perturbed_SVD = deflate_reduction(SVD_results, A, B);
    FLUSH;
    return (perturbed_SVD);
  }

  field<mat> orthogonalize_batch_effect(mat &S, field<mat> SVD_results,
                                        mat &design)
  {
    stdout_printf("Orthogonalizing batch effect: (dense):\n"); // fflush(stdout);

    mat Z = mat(S * design);
    gram_schmidt(Z);

    mat A = Z;
    mat B = -mat(trans(trans(Z) * S));

    field<mat> perturbed_SVD = deflate_reduction(SVD_results, A, B);
    FLUSH;
    return (perturbed_SVD);
  }

  field<mat> orthogonalize_basal(sp_mat &S, field<mat> SVD_results,
                                 mat &basal_state)
  {
    stdout_printf("Orthogonalizing basal (sparse):\n"); // fflush(stdout);

    mat Z = basal_state;
    gram_schmidt(Z);

    mat A = Z;
    mat B = -mat(trans(trans(Z) * S));

    field<mat> perturbed_SVD = deflate_reduction(SVD_results, A, B);
    FLUSH;
    return (perturbed_SVD);
  }

  field<mat> orthogonalize_basal(mat &S, field<mat> SVD_results,
                                 mat &basal_state)
  {
    stdout_printf("Orthogonalizing basal (dense):\n"); // fflush(stdout);

    mat Z = basal_state;
    gram_schmidt(Z);

    mat A = Z;
    mat B = -mat(trans(trans(Z) * S));

    printf("PSVD now");

    field<mat> perturbed_SVD = deflate_reduction(SVD_results, A, B);
    FLUSH;
    return (perturbed_SVD);
  }
} // namespace ACTIONet
