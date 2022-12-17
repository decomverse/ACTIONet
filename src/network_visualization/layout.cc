#include <ACTIONet.h>

#include <atomic>
#include <cfloat>
#include <thread>

#include "coords.h"
#include "epoch.h"
#include "gradient.h"
#include "optimize.h"
#include "rng.h"
#include "rparallel.h"
#include "sampler.h"

const double UMAP_A[101] = {
    1.93280839781719, 1.89560586588002, 1.85873666431227, 1.82221007490834,
    1.78603612060048, 1.75022496320214, 1.71478579945151, 1.67972997626197,
    1.64506544270902, 1.610800661285, 1.57694346052399, 1.54350101780511,
    1.51047986323257, 1.47788588612333, 1.44572435168023, 1.41399925414561,
    1.38271638006498, 1.35187804260518, 1.3214872860387, 1.29154663185922,
    1.26205810311418, 1.23302325071067, 1.20444317424075, 1.17631854866857,
    1.14864964274379, 1.12143634262879, 1.09467817152021, 1.0683743100033,
    1.04252361298475, 1.01712481754341, 0.992175611624647, 0.967674513244996,
    0.943619207179927, 0.920007077834315, 0.896835219021839, 0.874100443595699,
    0.851800999392949, 0.829931994792615, 0.808490430178554, 0.787472613514984,
    0.766873638278737, 0.746690990400437, 0.726919886947928, 0.707556026044195,
    0.688594985599233, 0.670032232635194, 0.651864066568649, 0.634084192553475,
    0.616688494561969, 0.599672088669339, 0.583030020204371, 0.5667572718654,
    0.550848768322639, 0.535299383967892, 0.520103947257001, 0.505257246260431,
    0.490754031684977, 0.476589022213249, 0.46275690208242, 0.449252325341552,
    0.436069912245555, 0.423205974605747, 0.4106531652521, 0.39840668039948,
    0.386461380891047, 0.374811984314975, 0.363453224264704, 0.352379851902848,
    0.341586644916259, 0.331068403184832, 0.320819956874279, 0.31083616902857,
    0.301110995958752, 0.291641183389757, 0.282420831386121, 0.273444955588216,
    0.264708614833586, 0.256206914916444, 0.247935008593902, 0.239888099677924,
    0.232061441819675, 0.224450342118235, 0.217050162160312, 0.209856317524031,
    0.202864281204524, 0.196069583611474, 0.189467814398248, 0.183054621446351,
    0.176825713015038, 0.17077685928726, 0.164903890637922, 0.159202699934773,
    0.153669242222215, 0.148299535941784, 0.143089661250278, 0.138035764053223,
    0.133134049958711, 0.12838079222654, 0.123772324007265, 0.119305671122251,
    0.114976081494676};
const double UMAP_B[101] = {
    0.790494973419029, 0.80063784415826, 0.810876441425738, 0.821199202674006,
    0.831595366275022, 0.84205539236769, 0.852571713401325, 0.863135518043442,
    0.873741680140683, 0.884384956993888, 0.895060878257082, 0.905765637284042,
    0.916495998501859, 0.927249214280422, 0.938022954467018, 0.948815759038301,
    0.95962499558526, 0.970449732070657, 0.981288783823989, 0.992141168965973,
    1.00300608092206, 1.01388286515112, 1.02477099750548, 1.03567006898871,
    1.04657977025277, 1.05749987674998, 1.06843023939592, 1.07937077470387,
    1.09032145585694, 1.10128169075827, 1.11225322117536, 1.12323470900213,
    1.13422639755358, 1.14522861434516, 1.15624176559097, 1.16726633179917,
    1.17830241385901, 1.18934945144456, 1.20040819996369, 1.21147891097075,
    1.22256381651844, 1.23366041866219, 1.24477022428392, 1.2558936051142,
    1.26703094885274, 1.27818265467871, 1.28934756395537, 1.30052872175886,
    1.31172539107843, 1.32293800168803, 1.3341669930459, 1.34541281413396,
    1.35667592718974, 1.36795680610473, 1.37925594017143, 1.39057383474783,
    1.40191101858967, 1.41326804557094, 1.42464550789942, 1.436044048272,
    1.44746436980037, 1.45890393087319, 1.47036701291879, 1.48185337703821,
    1.49336326709497, 1.50489726618312, 1.51645596605121, 1.52803997486173,
    1.53964990048402, 1.55128637349183, 1.56295003156298, 1.57464152150044,
    1.58636409305622, 1.59811350189048, 1.60989278253114, 1.62170263415549,
    1.63354377154668, 1.64541692037945, 1.65732282325244, 1.66926223230814,
    1.68123591907029, 1.69324466615879, 1.70528927262371, 1.71737055545595,
    1.72948934595558, 1.74164649289645, 1.75384285823827, 1.76607932576738,
    1.77835679827623, 1.79067619009556, 1.80303844043406, 1.81544450541945,
    1.82789536263139, 1.84039200538657, 1.85293545544251, 1.86552674229068,
    1.87816693701183, 1.89085711093115, 1.90359837758981, 1.91638829237987,
    1.92923479503841};

#define NEGATIVE_SAMPLE_RATE 3.0
#define UMAP_SEED 0
#define GAMMA 1.0
#define ADAM_ALPHA 1.0 /*same as learning_rate*/
#define ADAM_BETA1 0.5 /*only adam: between 0 and 1*/
#define ADAM_BETA2 0.9 /*only adam: between 0 and 1*/
#define ADAM_EPS 1e-7  /*only adam: between 1e-8 and 1e-3*/

sp_mat smoothKNN(sp_mat &D, int max_iter = 64, double epsilon = 1e-6,
                 double bandwidth = 1.0, double local_connectivity = 1.0,
                 double min_k_dist_scale = 1e-3, double min_sim = 1e-8,
                 int thread_no = 0)
{
  ;

  int nV = D.n_cols;
  sp_mat G = D;

  //#pragma omp parallel for num_threads(thread_no)
  for (int i = 0; i < nV; i++)
  {
    //  ParallelFor(0, nV, thread_no, [&](size_t i, size_t threadId) {
    sp_mat v = D.col(i);
    vec vals = nonzeros(v);
    if (vals.n_elem > local_connectivity)
    {
      double rho = min(vals);
      vec negated_shifted_vals = -(vals - rho);
      uvec deflate = find(vals <= rho);
      negated_shifted_vals(deflate).zeros();

      double target = log2(vals.n_elem + 1);

      // Binary search to find optimal sigma
      double sigma = 1.0;
      double lo = 0.0;
      double hi = DBL_MAX;

      int j;
      for (j = 0; j < max_iter; j++)
      {
        double obj = sum(exp(negated_shifted_vals / sigma));

        if (abs(obj - target) < epsilon)
        {
          break;
        }

        if (target < obj)
        {
          hi = sigma;
          sigma = 0.5 * (lo + hi);
        }
        else
        {
          lo = sigma;
          if (hi == DBL_MAX)
          {
            sigma *= 2;
          }
          else
          {
            sigma = 0.5 * (lo + hi);
          }
        }
      }

      // double obj = sum(exp(negated_shifted_vals / sigma));
      double mean_dist = mean(mean_dist);
      sigma = std::max(min_k_dist_scale * mean_dist, sigma);

      for (sp_mat::col_iterator it = G.begin_col(i); it != G.end_col(i); ++it)
      {
        *it = max(min_sim, exp(-max(0.0, (*it) - rho) / (sigma * bandwidth)));
      }
    }
    else
    {
      for (sp_mat::col_iterator it = G.begin_col(i); it != G.end_col(i); ++it)
      {
        *it = 1.0;
      }
    }
  }

  return (G);
}

// Template class specialization to handle different rng/batch combinations
template <bool DoBatch = true>
struct BatchRngFactory
{
  using PcgFactoryType = batch_pcg_factory;
  using TauFactoryType = batch_tau_factory;
};
template <>
struct BatchRngFactory<false>
{
  using PcgFactoryType = pcg_factory;
  using TauFactoryType = tau_factory;
};

struct UmapFactory
{
  bool move_other;
  bool pcg_rand;
  std::vector<float> &head_embedding;
  std::vector<float> &tail_embedding;
  const std::vector<unsigned int> &positive_head;
  const std::vector<unsigned int> &positive_tail;
  const std::vector<unsigned int> &positive_ptr;
  unsigned int n_epochs;
  unsigned int n_head_vertices;
  unsigned int n_tail_vertices;
  const std::vector<float> &epochs_per_sample;
  float initial_alpha;
  float negative_sample_rate;
  bool batch;
  std::size_t n_threads;
  std::size_t grain_size;
  string opt_name;
  double alpha, beta1, beta2, eps;
  std::mt19937_64 engine;

  UmapFactory(bool move_other, bool pcg_rand,
              std::vector<float> &head_embedding,
              std::vector<float> &tail_embedding,
              const std::vector<unsigned int> &positive_head,
              const std::vector<unsigned int> &positive_tail,
              const std::vector<unsigned int> &positive_ptr,
              unsigned int n_epochs, unsigned int n_head_vertices,
              unsigned int n_tail_vertices,
              const std::vector<float> &epochs_per_sample, float initial_alpha,
              float negative_sample_rate, bool batch, std::size_t n_threads,
              std::size_t grain_size, string opt_name, double alpha,
              double beta1, double beta2, double eps, std::mt19937_64 &engine)
      : move_other(move_other),
        pcg_rand(pcg_rand),
        head_embedding(head_embedding),
        tail_embedding(tail_embedding),
        positive_head(positive_head),
        positive_tail(positive_tail),
        positive_ptr(positive_ptr),
        n_epochs(n_epochs),
        n_head_vertices(n_head_vertices),
        n_tail_vertices(n_tail_vertices),
        epochs_per_sample(epochs_per_sample),
        initial_alpha(initial_alpha),
        negative_sample_rate(negative_sample_rate),
        batch(batch),
        n_threads(n_threads),
        grain_size(grain_size),
        alpha(alpha),
        beta1(beta1),
        beta2(beta2),
        eps(eps),
        engine(engine),
        opt_name(opt_name) {}

  template <typename Gradient>
  void create(const Gradient &gradient)
  {
    if (move_other)
    {
      create_impl<true>(gradient, pcg_rand, batch);
    }
    else
    {
      create_impl<false>(gradient, pcg_rand, batch);
    }
  }

  template <bool DoMove, typename Gradient>
  void create_impl(const Gradient &gradient, bool pcg_rand, bool batch)
  {
    if (batch)
    {
      create_impl<BatchRngFactory<true>, DoMove>(gradient, pcg_rand, batch);
    }
    else
    {
      create_impl<BatchRngFactory<false>, DoMove>(gradient, pcg_rand, batch);
    }
  }

  template <typename BatchRngFactory, bool DoMove, typename Gradient>
  void create_impl(const Gradient &gradient, bool pcg_rand, bool batch)
  {
    if (pcg_rand)
    {
      create_impl<typename BatchRngFactory::PcgFactoryType, DoMove>(gradient,
                                                                    batch);
    }
    else
    {
      create_impl<typename BatchRngFactory::TauFactoryType, DoMove>(gradient,
                                                                    batch);
    }
  }

  uwot::Adam create_adam()
  {
    return uwot::Adam(alpha, beta1, beta2, eps, head_embedding.size());
  }

  uwot::Sgd create_sgd() { return uwot::Sgd(alpha); }

  template <typename RandFactory, bool DoMove, typename Gradient>
  void create_impl(const Gradient &gradient, bool batch)
  {
    if (batch)
    {
      if (opt_name == "adam")
      {
        auto opt = create_adam();
        create_impl_batch_opt<decltype(opt), RandFactory, DoMove, Gradient>(
            gradient, opt, batch);
      }
      else if (opt_name == "sgd")
      {
        auto opt = create_sgd();
        create_impl_batch_opt<decltype(opt), RandFactory, DoMove, Gradient>(
            gradient, opt, batch);
      }
      else
      {
        stderr_printf("Unknown optimization method: %s\n", opt_name.c_str());
        FLUSH;
        return;
      }
    }
    else
    {
      const std::size_t ndim = head_embedding.size() / n_head_vertices;
      uwot::Sampler sampler(epochs_per_sample, negative_sample_rate);
      uwot::InPlaceUpdate<DoMove> update(head_embedding, tail_embedding,
                                         initial_alpha);
      uwot::EdgeWorker<Gradient, decltype(update), RandFactory> worker(
          gradient, update, positive_head, positive_tail, sampler, ndim,
          n_tail_vertices, n_threads, engine);
      create_impl(worker, gradient);
    }
  }

  template <typename Opt, typename RandFactory, bool DoMove, typename Gradient>
  void create_impl_batch_opt(const Gradient &gradient, Opt &opt, bool batch)
  {
    uwot::Sampler sampler(epochs_per_sample, negative_sample_rate);
    const std::size_t ndim = head_embedding.size() / n_head_vertices;
    uwot::BatchUpdate<DoMove, decltype(opt)> update(head_embedding,
                                                    tail_embedding, opt);
    uwot::NodeWorker<Gradient, decltype(update), RandFactory> worker(
        gradient, update, positive_head, positive_tail, positive_ptr, sampler,
        ndim, n_tail_vertices, engine);
    create_impl(worker, gradient);
  }

  template <typename Worker, typename Gradient>
  void create_impl(Worker &worker, const Gradient &gradient)
  {
    if (n_threads > 0)
    {
      RParallel parallel(n_threads, grain_size);
      create_impl(worker, gradient, parallel);
    }
    else
    {
      RSerial serial;
      create_impl(worker, gradient, serial);
    }
  }

  template <typename Worker, typename Gradient, typename Parallel>
  void create_impl(Worker &worker, const Gradient &gradient,
                   Parallel &parallel)
  {
    uwot::optimize_layout(worker, n_epochs, parallel);
  }
};

void create_umap(UmapFactory &umap_factory, double a, double b, double gamma,
                 bool approx_pow)
{
  if (approx_pow)
  {
    const uwot::apumap_gradient gradient(a, b, gamma);
    umap_factory.create(gradient);
  }
  else
  {
    const uwot::umap_gradient gradient(a, b, gamma);
    umap_factory.create(gradient);
  }
}

void create_tumap(UmapFactory &umap_factory)
{
  const uwot::tumap_gradient gradient;
  umap_factory.create(gradient);
}

void create_pacmap(UmapFactory &umap_factory, double a, double b)
{
  const uwot::pacmap_gradient gradient(a, b);
  umap_factory.create(gradient);
}

void create_largevis(UmapFactory &umap_factory, double gamma)
{
  const uwot::largevis_gradient gradient(gamma);
  umap_factory.create(gradient);
}

template <typename Float>
std::pair<Float, Float> find_ab(Float spread, Float min_dist, Float grid = 300,
                                Float limit = 0.5, int iter = 50,
                                Float tol = 1e-6)
{
  Float x_half = std::log(limit) * -spread + min_dist;
  Float d_half = limit / -spread;

  // Compute the x and y coordinates of the expected distance curve.
  std::vector<Float> grid_x(grid), grid_y(grid), log_x(grid);
  const Float delta = spread * 3 / grid;
  for (int g = 0; g < grid; ++g)
  {
    grid_x[g] =
        (g + 1) * delta; // +1 to avoid meaningless least squares result at x =
                         // 0, where both curves have y = 1 (and also the
                         // derivative w.r.t. b is not defined).
    log_x[g] = std::log(grid_x[g]);
    grid_y[g] =
        (grid_x[g] <= min_dist ? 1
                               : std::exp(-(grid_x[g] - min_dist) / spread));
  }

  // Starting estimates.
  Float b = -d_half * x_half / (1 / limit - 1) / (2 * limit * limit);
  Float a = (1 / limit - 1) / std::pow(x_half, 2 * b);

  std::vector<Float> observed_y(grid), xpow(grid);
  auto compute_ss = [&](Float A, Float B) -> Float
  {
    for (int g = 0; g < grid; ++g)
    {
      xpow[g] = std::pow(grid_x[g], 2 * B);
      observed_y[g] = 1 / (1 + A * xpow[g]);
    }

    Float ss = 0;
    for (int g = 0; g < grid; ++g)
    {
      ss += (grid_y[g] - observed_y[g]) * (grid_y[g] - observed_y[g]);
    }

    return ss;
  };
  Float ss = compute_ss(a, b);

  for (int it = 0; it < iter; ++it)
  {
    // Computing the first and second derivatives of the sum of squared
    // differences.
    Float da = 0, db = 0, daa = 0, dab = 0, dbb = 0;
    for (int g = 0; g < grid; ++g)
    {
      const Float &x = grid_x[g];
      const Float &gy = grid_y[g];
      const Float &oy = observed_y[g];

      const Float &x2b = xpow[g];
      const Float logx2 = log_x[g] * 2;
      const Float delta = oy - gy;

      // -(2 * (x^(2 * b)/(1 + a * x^(2 * b))^2 * (1/(1 + a * x^(2 * b)) - y)))
      da += -2 * x2b * oy * oy * delta;

      // -(2 * (a * (x^(2 * b) * (log(x) * 2))/(1 + a * x^(2 * b))^2 * (1/(1 + a
      // * x^(2 * b)) - y)))
      db += -2 * a * x2b * logx2 * oy * oy * delta;

      // 2 * (
      //     x^(2 * b)/(1 + a * x^(2 * b))^2 * (x^(2 * b)/(1 + a * x^(2 * b))^2)
      //     + x^(2 * b) * (2 * (x^(2 * b) * (1 + a * x^(2 * b))))/((1 + a *
      //     x^(2 * b))^2)^2 * (1/(1 + a * x^(2 * b)) - y)
      // )
      daa += 2 * (x2b * oy * oy * x2b * oy * oy +
                  x2b * 2 * x2b * oy * oy * oy * delta);

      //-(2 *
      //    (
      //        (
      //            (x^(2 * b) * (log(x) * 2))/(1 + a * x^(2 * b))^2
      //            - a * (x^(2 * b) * (log(x) * 2)) * (2 * (x^(2 * b) * (1 + a
      //            * x^(2 * b))))/((1 + a * x^(2 * b))^2)^2
      //        )
      //        * (1/(1 + a * x^(2 * b)) - y)
      //        - a * (x^(2 * b) * (log(x) * 2))/(1 + a * x^(2 * b))^2 * (x^(2 *
      //        b)/(1 + a * x^(2 * b))^2)
      //    )
      //)
      dab +=
          -2 *
          ((x2b * logx2 * oy * oy - a * x2b * logx2 * 2 * x2b * oy * oy * oy) *
               delta -
           a * x2b * logx2 * oy * oy * x2b * oy * oy);

      // -(2 *
      //     (
      //         (
      //             a * (x^(2 * b) * (log(x) * 2) * (log(x) * 2))/(1 + a * x^(2
      //             * b))^2
      //             - a * (x^(2 * b) * (log(x) * 2)) * (2 * (a * (x^(2 * b) *
      //             (log(x) * 2)) * (1 + a * x^(2 * b))))/((1 + a * x^(2 *
      //             b))^2)^2
      //         )
      //         * (1/(1 + a * x^(2 * b)) - y)
      //         - a * (x^(2 * b) * (log(x) * 2))/(1 + a * x^(2 * b))^2 * (a *
      //         (x^(2 * b) * (log(x) * 2))/(1 + a * x^(2 * b))^2)
      //     )
      // )
      dbb += -2 * (((a * x2b * logx2 * logx2 * oy * oy) -
                    (a * x2b * logx2 * 2 * a * x2b * logx2 * oy * oy * oy)) *
                       delta -
                   a * x2b * logx2 * oy * oy * a * x2b * logx2 * oy * oy);
    }

    // Applying the Newton iterations with damping.
    Float determinant = daa * dbb - dab * dab;
    const Float delta_a = (da * dbb - dab * db) / determinant;
    const Float delta_b = (-da * dab + daa * db) / determinant;

    Float ss_next = 0;
    Float factor = 1;
    for (int inner = 0; inner < 10; ++inner, factor /= 2)
    {
      ss_next = compute_ss(a - factor * delta_a, b - factor * delta_b);
      if (ss_next < ss)
      {
        break;
      }
    }

    if (ss && 1 - ss_next / ss > tol)
    {
      a -= factor * delta_a;
      b -= factor * delta_b;
      ss = ss_next;
    }
    else
    {
      break;
    }
  }

  return std::make_pair(a, b);
}

namespace ACTIONet
{

  field<mat> layoutNetwork_xmap(sp_mat &G, mat &initial_position,
                                bool presmooth_network, const std::string &method,
                                double min_dist, double spread, double gamma,
                                unsigned int n_epochs, int thread_no, int seed,
                                double learning_rate, int sim2dist)
  {
    if (thread_no <= 0)
    {
      thread_no = SYS_THREADS_DEF;
    }

    auto found = find_ab(spread, min_dist);
    double a = found.first;
    double b = found.second;

    // a = UMAP_A[50];
    // b = UMAP_B[50];

    stdout_printf(
        "Laying-out input network: method = %s, a = %.3f, b = %.3f (epochs = %d, "
        "threads=%d)\n",
        method.c_str(), a, b, n_epochs, thread_no);

    bool move_other = true;
    std::size_t grain_size = 1;
    bool pcg_rand = true;
    bool approx_pow = true;
    bool batch = true;
    string opt_name = "adam";
    double alpha = ADAM_ALPHA, beta1 = ADAM_BETA1, beta2 = ADAM_BETA2,
           eps = ADAM_EPS, negative_sample_rate = NEGATIVE_SAMPLE_RATE;

    field<mat> res(3);
    std::mt19937_64 engine(seed);

    mat init_coors;
    if (initial_position.n_rows != G.n_rows)
    {
      stderr_printf(
          "Number of rows in the initial_position should match with the number "
          "of vertices in G\n");
      FLUSH;
      return (res);
    }

    // Encode positive edges of the graph
    sp_mat H = G;
    if (presmooth_network == true)
    {
      stdout_printf("%\tSmoothing similarities (sim2dist = %d) ... ", sim2dist);
      if (sim2dist == 1)
      {
        H.for_each([](sp_mat::elem_type &val)
                   { val = 1 - val; });
      }
      else if (sim2dist == 2)
      {
        H.for_each([](sp_mat::elem_type &val)
                   { val = (1 - val) * (1 - val); });
      }
      else if (sim2dist == 3)
      {
        H.for_each([](sp_mat::elem_type &val)
                   { val = -log(val); });
      }
      else
      {
        H.for_each([](sp_mat::elem_type &val)
                   { val = 1 / val; });
      }

      int max_iter = 64;
      double epsilon = 1e-6, bandwidth = 1.0, local_connectivity = 1.0,
             min_k_dist_scale = 1e-3, min_sim = 1e-8;

      H = smoothKNN(H, max_iter, epsilon, bandwidth, local_connectivity,
                    min_k_dist_scale, min_sim, thread_no);
      stdout_printf("done\n");
    }

    double w_max = max(max(H));
    H.clean(w_max / n_epochs);

    sp_mat Ht = trans(H);
    Ht.sync();

    unsigned int nV = H.n_rows;
    unsigned int nE = H.n_nonzero;
    unsigned int nD = init_coors.n_rows;

    vector<unsigned int> positive_head(nE);
    vector<unsigned int> positive_tail(nE);
    vector<float> epochs_per_sample(nE);

    std::vector<unsigned int> positive_ptr(Ht.n_cols + 1);

    int i = 0;
    if (batch == false)
    {
      for (sp_mat::iterator it = H.begin(); it != H.end(); ++it)
      {
        epochs_per_sample[i] = w_max / (*it);
        positive_head[i] = it.row();
        positive_tail[i] = it.col();
        i++;
      }
    }
    else
    {
      for (sp_mat::iterator it = Ht.begin(); it != Ht.end(); ++it)
      {
        epochs_per_sample[i] = w_max / (*it);
        positive_tail[i] = it.row();
        positive_head[i] = it.col();
        i++;
      }
      for (int k = 0; k < Ht.n_cols + 1; k++)
      {
        positive_ptr[k] = Ht.col_ptrs[k];
      }
    }

    mat coors2D = initial_position.cols(0, 1);
    init_coors = trans(zscore(coors2D));

    // Initial coordinates of vertices (0-simplices)
    vector<float> head_embedding(init_coors.n_elem);
    fmat sub_coor = conv_to<fmat>::from(init_coors);
    memcpy(head_embedding.data(), sub_coor.memptr(),
           sizeof(float) * head_embedding.size());
    // vector<float> tail_embedding(head_embedding);
    uwot::Coords coords = uwot::Coords(head_embedding);

    UmapFactory umap_factory(
        move_other, pcg_rand, coords.get_head_embedding(),
        coords.get_tail_embedding(), positive_head, positive_tail, positive_ptr,
        n_epochs, nV, nV, epochs_per_sample, learning_rate, negative_sample_rate,
        batch, thread_no, grain_size, opt_name, alpha, beta1, beta2, eps, engine);

    stdout_printf("Computing 2D layout ... ");
    FLUSH;
    if (method == "umap")
    {
      create_umap(umap_factory, a, b, gamma, approx_pow);
    }
    else if (method == "tumap")
    {
      create_tumap(umap_factory);
    }
    else if (method == "largevis")
    {
      create_largevis(umap_factory, gamma);
    }
    else if (method == "pacmap")
    {
      create_pacmap(umap_factory, a, b);
    }
    else
    {
      stderr_printf("Unknown method: %s\n", method.c_str());
      FLUSH;
      return (res);
    }
    fmat coordinates_float(coords.get_head_embedding().data(), 2, nV);
    mat coordinates_2D = trans(conv_to<mat>::from(coordinates_float));
    stdout_printf("done\n");
    FLUSH;

    mat coordinates_3D = zeros(nV, 3);
    mat RGB_colors = zeros(nV, 3);
    if (initial_position.n_cols > 2)
    {
      mat coors3D = join_rows(coordinates_2D, initial_position.col(2));
      init_coors = trans(zscore(coors3D));
      head_embedding.clear();
      head_embedding.resize(nV * 3);
      sub_coor = conv_to<fmat>::from(init_coors);
      memcpy(head_embedding.data(), sub_coor.memptr(),
             sizeof(float) * head_embedding.size());
      coords = uwot::Coords(head_embedding);

      UmapFactory umap_factory_3D(
          move_other, pcg_rand, coords.get_head_embedding(),
          coords.get_tail_embedding(), positive_head, positive_tail, positive_ptr,
          n_epochs / 2, nV, nV, epochs_per_sample, learning_rate,
          negative_sample_rate, batch, thread_no, grain_size, opt_name, alpha,
          beta1, beta2, eps, engine);

      stdout_printf("Computing 3D layout ... ");
      FLUSH;
      if (method == "umap")
      {
        create_umap(umap_factory_3D, a, b, gamma, approx_pow);
      }
      else if (method == "tumap")
      {
        create_tumap(umap_factory_3D);
      }
      else if (method == "largevis")
      {
        create_largevis(umap_factory_3D, gamma);
      }
      else if (method == "pacmap")
      {
        create_pacmap(umap_factory_3D, a, b);
      }
      else
      {
        stderr_printf("Unknown method: %s\n", method.c_str());
        FLUSH;
        return (res);
      }
      fmat coordinates_float(coords.get_head_embedding().data(), 3, nV);
      coordinates_3D = trans(conv_to<mat>::from(coordinates_float));
      stdout_printf("done\n");
      FLUSH;

      stdout_printf("Computing de novo node colors ... ");
      FLUSH;
      mat U;
      vec s;
      mat V;
      svd_econ(U, s, V, coordinates_3D, "left", "std");

      mat Z = normalize_scores(U.cols(0, 2), 1, thread_no);

      vec a = 75 * Z.col(0);
      vec b = 75 * Z.col(1);

      vec L = Z.col(2);
      L = 25.0 + 70.0 * (L - min(L)) / (max(L) - min(L));

      double r_channel, g_channel, b_channel;
      for (int i = 0; i < nV; i++)
      {
        Lab2Rgb(&r_channel, &g_channel, &b_channel, L(i), a(i), b(i));

        RGB_colors(i, 0) = min(1.0, max(0.0, r_channel));
        RGB_colors(i, 1) = min(1.0, max(0.0, g_channel));
        RGB_colors(i, 2) = min(1.0, max(0.0, b_channel));
      }
      stdout_printf("done\n");
      FLUSH;
    }

    res(0) = coordinates_2D;
    res(1) = coordinates_3D;
    res(2) = RGB_colors;

    return (res);
  }

  mat transform_layout(sp_mat &G, mat &reference_layout, bool presmooth_network,
                       const std::string &method, double min_dist, double spread,
                       double gamma, unsigned int n_epochs, int thread_no,
                       int seed, double learning_rate, int sim2dist)
  {
    mat coordinates;
    if (thread_no <= 0)
    {
      thread_no = SYS_THREADS_DEF;
    }

    auto found = find_ab(spread, min_dist);
    double a = found.first;
    double b = found.second;

    // a = UMAP_A[50];
    // b = UMAP_B[50];

    if (reference_layout.n_rows != G.n_cols)
    {
      stderr_printf(
          "Number of rows in the reference_layout should match with the number "
          "of columns in G\n");
      FLUSH;
      return (coordinates);
    }

    int Nq = G.n_rows, Nr = G.n_cols, D = reference_layout.n_cols;
    stdout_printf(
        "Transforming graph G with %d vertices, using a reference of %d "
        "vertices, in a %dD dimensions (%d threads)\n",
        Nq, Nr, D, thread_no);
    stdout_printf("\tmethod = %s, a = %.3f, b = %.3f (epochs = %d, threads=%d)\n",
                  method.c_str(), a, b, n_epochs, thread_no);

    sp_mat W = normalize_adj(G, 1);
    mat query_layout = spmat_mat_product_parallel(W, reference_layout, thread_no);

    bool move_other = false;
    std::size_t grain_size = 1;
    bool pcg_rand = true;
    bool approx_pow = true;
    bool batch = true;
    string opt_name = "adam";
    double alpha = ADAM_ALPHA, beta1 = ADAM_BETA1, beta2 = ADAM_BETA2,
           eps = ADAM_EPS, negative_sample_rate = NEGATIVE_SAMPLE_RATE;

    field<mat> res(3);
    std::mt19937_64 engine(seed);

    // Encode positive edges of the graph
    sp_mat H = G;

    double w_max = max(max(H));
    H.clean(w_max / n_epochs);

    sp_mat Ht = trans(H);
    Ht.sync();

    unsigned int nE = H.n_nonzero;
    vector<unsigned int> positive_head(nE);
    vector<unsigned int> positive_tail(nE);
    vector<float> epochs_per_sample(nE);

    std::vector<unsigned int> positive_ptr(Ht.n_cols + 1);

    int i = 0;
    if (batch == false)
    {
      for (sp_mat::iterator it = H.begin(); it != H.end(); ++it)
      {
        epochs_per_sample[i] = w_max / (*it);
        positive_head[i] = it.row();
        positive_tail[i] = it.col();
        i++;
      }
    }
    else
    {
      for (sp_mat::iterator it = Ht.begin(); it != Ht.end(); ++it)
      {
        epochs_per_sample[i] = w_max / (*it);
        positive_tail[i] = it.row();
        positive_head[i] = it.col();
        i++;
      }
      for (int k = 0; k < Ht.n_cols + 1; k++)
      {
        positive_ptr[k] = Ht.col_ptrs[k];
      }
    }

    query_layout = trans(query_layout);
    reference_layout = trans(reference_layout);

    // Initial coordinates of vertices (0-simplices)
    vector<float> head_embedding(query_layout.n_elem);
    fmat sub_coor = conv_to<fmat>::from(query_layout);
    memcpy(head_embedding.data(), sub_coor.memptr(),
           sizeof(float) * head_embedding.size());
    vector<float> tail_embedding(reference_layout.n_elem);
    sub_coor = conv_to<fmat>::from(reference_layout);
    memcpy(tail_embedding.data(), sub_coor.memptr(),
           sizeof(float) * tail_embedding.size());
    uwot::Coords coords = uwot::Coords(head_embedding, tail_embedding);

    UmapFactory umap_factory(
        move_other, pcg_rand, coords.get_head_embedding(),
        coords.get_tail_embedding(), positive_head, positive_tail, positive_ptr,
        n_epochs, Nq, Nr, epochs_per_sample, learning_rate, negative_sample_rate,
        batch, thread_no, grain_size, opt_name, alpha, beta1, beta2, eps, engine);

    stdout_printf("Transforming layout ... ");
    FLUSH;
    if (method == "umap")
    {
      create_umap(umap_factory, a, b, gamma, approx_pow);
    }
    else if (method == "tumap")
    {
      create_tumap(umap_factory);
    }
    else if (method == "largevis")
    {
      create_largevis(umap_factory, gamma);
    }
    else if (method == "pacmap")
    {
      create_pacmap(umap_factory, a, b);
    }
    else
    {
      stderr_printf("Unknown method: %s\n", method.c_str());
      FLUSH;
      return (coordinates);
    }
    fmat coordinates_float(coords.get_head_embedding().data(), 2, Nq);
    coordinates = trans(conv_to<mat>::from(coordinates_float));
    stdout_printf("done\n");
    FLUSH;

    return (coordinates);
  }

} // namespace ACTIONet
