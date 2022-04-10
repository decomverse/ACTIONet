#include <ACTIONet.h>

#include <atomic>
#include <cfloat>
#include <thread>

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

#define NEGATIVE_SAMPLE_RATE 5.0
#define LEARNING_RATE 1.0
#define UMAP_SEED 0
#define GAMMA 1.0

std::vector<float> optimize_layout_umap(
    std::vector<float> head_vec, std::vector<float> tail_vec,
    const std::vector<unsigned int> positive_head,
    const std::vector<unsigned int> positive_tail, unsigned int n_epochs,
    unsigned int n_vertices, const std::vector<float> epochs_per_sample,
    double a, double b, double gamma, double initial_alpha,
    double negative_sample_rate, bool approx_pow, std::size_t n_threads,
    std::size_t grain_size, bool move_other, int seed);

std::vector<float> optimize_layout_tumap(
    std::vector<float> head_vec, std::vector<float> tail_vec,
    const std::vector<unsigned int> positive_head,
    const std::vector<unsigned int> positive_tail, unsigned int n_epochs,
    unsigned int n_vertices, const std::vector<float> epochs_per_sample,
    double initial_alpha, double negative_sample_rate,
    std::size_t thread_no = 0, std::size_t grain_size = 1,
    bool move_other = true, int seed = 0);

std::vector<float> optimize_layout_largevis(
    std::vector<float> head_vec, const std::vector<unsigned int> positive_head,
    const std::vector<unsigned int> positive_tail, unsigned int n_epochs,
    unsigned int n_vertices, const std::vector<float> epochs_per_sample,
    double gamma, double initial_alpha, double negative_sample_rate,
    std::size_t thread_no = 0, std::size_t grain_size = 1, int seed = 0);

namespace ACTIONet
{

  sp_mat smoothKNN(sp_mat D, int thread_no = -1)
  {
    double epsilon = 1e-6;

    int nV = D.n_cols;
    sp_mat G = D;

    //#pragma omp parallel for num_threads(thread_no)
    for (int i = 0; i < nV; i++)
    {
      //  ParallelFor(0, nV, thread_no, [&](size_t i, size_t threadId) {
      sp_mat v = D.col(i);
      vec vals = nonzeros(v);
      if (vals.n_elem > 0)
      {
        double rho = min(vals);
        vec negated_shifted_vals = -(vals - rho);
        double target = log2(vals.n_elem);

        // Binary search to find optimal sigma
        double sigma = 1.0;
        double lo = 0.0;
        double hi = DBL_MAX;

        int j;
        for (j = 0; j < 64; j++)
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

        double obj = sum(exp(negated_shifted_vals / sigma));

        for (sp_mat::col_iterator it = G.begin_col(i); it != G.end_col(i); ++it)
        {
          *it = max(1e-16, exp(-max(0.0, (*it) - rho) / sigma));
        }
      }
    }

    return (G);
  }

  field<mat> layoutNetwork_umap(sp_mat &G, mat initial_position, int compactness_level,
                                unsigned int n_epochs,
                                int layout_alg, int thread_no,
                                int seed)
  {
    if (thread_no <= 0)
    {
      thread_no = SYS_THREADS_DEF;
    }

    stdout_printf("Computing layout (%d threads):\n", thread_no); FLUSH;

    field<mat> res(3);

    mat init_coors = tzscoret(initial_position.rows(0, 2));

    sp_mat H = G;
    H.for_each([](sp_mat::elem_type &val)
               { val = 1.0 / val; });
    H = smoothKNN(H, thread_no);

    unsigned int nV = H.n_rows;

    if (layout_alg == UMAP_LAYOUT)
    {
      compactness_level = (compactness_level < 0) ? 0 : compactness_level;
      compactness_level = (compactness_level > 100) ? 100 : compactness_level;
      stdout_printf("\tParameters: compactness = %d, layout_epochs = %d\n",
                    compactness_level, n_epochs); FLUSH;
    }
    else
    {
      stdout_printf("\tParameters: layout_epochs = %d\n", n_epochs); FLUSH;
    }

    double a_param = UMAP_A[compactness_level];
    double b_param = UMAP_B[compactness_level];

    // linearized list of edges (1-simplices)
    unsigned int nE = H.n_nonzero;
    vector<unsigned int> positive_head(nE);
    vector<unsigned int> positive_tail(nE);
    vector<float> epochs_per_sample(nE);

    int i = 0;
    double w_max = max(max(H));
    for (sp_mat::iterator it = H.begin(); it != H.end(); ++it)
    {
      epochs_per_sample[i] = w_max / (*it);
      positive_head[i] = it.row();
      positive_tail[i] = it.col();
      i++;
    }

    // Initial coordinates of vertices (0-simplices)
    vector<float> head_vec(init_coors.n_cols * 2);
    fmat sub_coor = conv_to<fmat>::from(init_coors.rows(0, 1));
    float *ptr = sub_coor.memptr();
    memcpy(head_vec.data(), ptr, sizeof(float) * head_vec.size());
    vector<float> tail_vec(head_vec);

    stdout_printf("\tComputing 2D layout ... ");

    // Stores linearized coordinates [x1, y1, x2, y2, ...]
    vector<float> result;

    switch (layout_alg)
    {
    case TUMAP_LAYOUT:
      result = optimize_layout_tumap(
          head_vec, tail_vec, positive_head, positive_tail, n_epochs, nV,
          epochs_per_sample, LEARNING_RATE, NEGATIVE_SAMPLE_RATE, thread_no, 1,
          true, seed);
      break;
    case UMAP_LAYOUT:
      result = optimize_layout_umap(
          head_vec, tail_vec, positive_head, positive_tail, n_epochs, nV,
          epochs_per_sample, a_param, b_param, GAMMA, LEARNING_RATE,
          NEGATIVE_SAMPLE_RATE, false, thread_no, 1, true, seed);
      break;
    default:
      stderr_printf(
          "Unknown layout algorithm chosen (%d). Switching to TUMAP.\n",
          layout_alg);
      result = optimize_layout_tumap(
          head_vec, tail_vec, positive_head, positive_tail, n_epochs, nV,
          epochs_per_sample, LEARNING_RATE, NEGATIVE_SAMPLE_RATE, thread_no, 1,
          true, seed);
      break;
    }

    fmat coordinates_float(result.data(), 2, nV);
    mat coordinates = conv_to<mat>::from(coordinates_float);
    coordinates = trans(coordinates);

    stdout_printf("done\n"); FLUSH;

    /****************************
   *  Compute 3D Embedding	*
   ***************************/
    head_vec.clear();
    head_vec.resize(init_coors.n_cols * 3);
    sub_coor = conv_to<fmat>::from(tzscoret(join_vert(trans(coordinates), init_coors.row(2))));
    ptr = sub_coor.memptr();
    memcpy(head_vec.data(), ptr, sizeof(float) * head_vec.size());
    tail_vec = head_vec;

    stdout_printf("\tComputing 3D layout ... ");

    result.clear();

    switch (layout_alg)
    {
    case TUMAP_LAYOUT:
      result = optimize_layout_tumap(
          head_vec, tail_vec, positive_head, positive_tail, n_epochs, nV,
          epochs_per_sample, LEARNING_RATE, NEGATIVE_SAMPLE_RATE, thread_no, 1,
          true, seed);
      break;
    case UMAP_LAYOUT:
      result = optimize_layout_umap(
          head_vec, tail_vec, positive_head, positive_tail, n_epochs, nV,
          epochs_per_sample, a_param, b_param, GAMMA, LEARNING_RATE,
          NEGATIVE_SAMPLE_RATE, false, thread_no, 1, true, seed);
      break;
    default:
      result = optimize_layout_tumap(
          head_vec, tail_vec, positive_head, positive_tail, n_epochs, nV,
          epochs_per_sample, LEARNING_RATE, NEGATIVE_SAMPLE_RATE, thread_no, 1,
          true, seed);
      break;
    }

    fmat coordinates_3D_float(result.data(), 3, nV);
    mat coordinates_3D = conv_to<mat>::from(coordinates_3D_float);
    coordinates_3D = trans(coordinates_3D);

    stdout_printf("done\n"); FLUSH;

    /****************************
   *  Now compute node colors *
   ***************************/
    stdout_printf("\tComputing node colors ... ");

    mat U;
    vec s;
    mat V;
    svd_econ(U, s, V, coordinates_3D, "left", "std");

    mat Z = normalize_scores(U, 2, thread_no);

    vec a = 75 * Z.col(0);
    vec b = 75 * Z.col(1);

    vec L = Z.col(2);
    L = 25.0 + 70.0 * (L - min(L)) / (max(L) - min(L));

    double r_channel, g_channel, b_channel;
    mat RGB_colors = zeros(nV, 3);
    for (int i = 0; i < nV; i++)
    {
      Lab2Rgb(&r_channel, &g_channel, &b_channel, L(i), a(i), b(i));

      RGB_colors(i, 0) = min(1.0, max(0.0, r_channel));
      RGB_colors(i, 1) = min(1.0, max(0.0, g_channel));
      RGB_colors(i, 2) = min(1.0, max(0.0, b_channel));
    }

    stdout_printf("done\n"); FLUSH;

    res(0) = coordinates;
    res(1) = coordinates_3D;
    res(2) = RGB_colors;
    return res;
  }

  field<mat> layoutNetwork_forced_atlas(sp_mat &G, mat initial_position, int compactness_level,
                                        unsigned int n_epochs, int thread_no,
                                        int seed)
  {
    if (thread_no <= 0)
    {
      thread_no = SYS_THREADS_DEF;
    }

    stdout_printf("Computing layout (%d threads):\n", thread_no);
    field<mat> res(3);

    int D = initial_position.n_rows;
    mat init_coors = initial_position.rows(0, min(2, D));

    /****************************
   *  Compute 3D Embedding	*
   ***************************/
    if (2 < D)
    {
      // Empty
    }
    else
    {
      // Empty
    }
    /****************************
   *  Now compute node colors *
   ***************************/

    return res;
  }

  //G: MxM, inter_graph: MxN, reference_coordinates: DxN
  field<mat> transform_layout(sp_mat &G, sp_mat &inter_graph, mat reference_coordinates, int compactness_level,
                              unsigned int n_epochs,
                              int layout_alg, int thread_no,
                              int seed)
  {
    if (thread_no <= 0)
    {
      thread_no = SYS_THREADS_DEF;
    }

    int M = G.n_rows, N = inter_graph.n_cols, D = reference_coordinates.n_rows;
    stdout_printf("Transforming graph G with %d vertices, using a reference of %d vertices, in a %dD dimensions (%d threads)\n", M, N, D, thread_no);

    field<mat> res(3);

    mat interpolated_coordinates = reference_coordinates * trans(inter_graph);

    sp_mat H = G;
    H.for_each([](sp_mat::elem_type &val)
               { val = 1.0 / val; });
    H = smoothKNN(H, thread_no);

    unsigned int nV = H.n_rows;

    if (layout_alg == UMAP_LAYOUT)
    {
      compactness_level = (compactness_level < 0) ? 0 : compactness_level;
      compactness_level = (compactness_level > 100) ? 100 : compactness_level;
      stdout_printf("\tParameters: compactness = %d, layout_epochs = %d\n",
                    compactness_level, n_epochs); FLUSH;
    }
    else
    {
      stdout_printf("\tParameters: layout_epochs = %d\n", n_epochs); FLUSH;
    }

    double a_param = UMAP_A[compactness_level];
    double b_param = UMAP_B[compactness_level];

    // linearized list of edges (1-simplices)
    unsigned int nE = H.n_nonzero;
    vector<unsigned int> positive_head(nE);
    vector<unsigned int> positive_tail(nE);
    vector<float> epochs_per_sample(nE);

    int i = 0;
    double w_max = max(max(H));
    for (sp_mat::iterator it = H.begin(); it != H.end(); ++it)
    {
      epochs_per_sample[i] = w_max / (*it);
      positive_head[i] = it.row();
      positive_tail[i] = it.col();
      i++;
    }

    vector<float> head_vec(interpolated_coordinates.n_cols * 2);
    fmat sub_coor = conv_to<fmat>::from(interpolated_coordinates.rows(0, 1));
    float *ptr = sub_coor.memptr();
    memcpy(head_vec.data(), ptr, sizeof(float) * head_vec.size());

    vector<float> tail_vec(reference_coordinates.n_cols * 2);
    sub_coor = conv_to<fmat>::from(reference_coordinates.rows(0, 1));
    ptr = sub_coor.memptr();
    memcpy(tail_vec.data(), ptr, sizeof(float) * tail_vec.size());

    stdout_printf("\tComputing 2D layout ... ");
    // Stores linearized coordinates [x1, y1, x2, y2, ...]
    vector<float> result;

    switch (layout_alg)
    {
    case TUMAP_LAYOUT:
      result = optimize_layout_tumap(
          head_vec, tail_vec, positive_head, positive_tail, n_epochs, nV,
          epochs_per_sample, LEARNING_RATE, NEGATIVE_SAMPLE_RATE, thread_no, 1,
          false, seed);
      break;
    case UMAP_LAYOUT:
      result = optimize_layout_umap(
          head_vec, tail_vec, positive_head, positive_tail, n_epochs, nV,
          epochs_per_sample, a_param, b_param, GAMMA, LEARNING_RATE,
          NEGATIVE_SAMPLE_RATE, false, thread_no, 1, false, seed);
      break;
    default:
      stderr_printf(
          "Unknown layout algorithm chosen (%d). Switching to TUMAP.\n",
          layout_alg);
      result = optimize_layout_tumap(
          head_vec, tail_vec, positive_head, positive_tail, n_epochs, nV,
          epochs_per_sample, LEARNING_RATE, NEGATIVE_SAMPLE_RATE, thread_no, 1,
          true, seed);
      break;
    }

    fmat coordinates_float(result.data(), 2, nV);
    mat coordinates = conv_to<mat>::from(coordinates_float);
    coordinates = trans(coordinates);
    res(0) = coordinates;

    stdout_printf("done\n"); FLUSH;

    if (D == 3)
    {
      /****************************
     *  Compute 3D Embedding	*
     ***************************/
      stdout_printf("\tComputing 3D layout ... ");
      result.clear();

      head_vec.clear();
      head_vec.resize(interpolated_coordinates.n_cols * 3);
      sub_coor = conv_to<fmat>::from(interpolated_coordinates.rows(0, 2));
      ptr = sub_coor.memptr();
      memcpy(head_vec.data(), ptr, sizeof(float) * head_vec.size());

      tail_vec.clear();
      tail_vec.resize(reference_coordinates.n_cols * 3);
      sub_coor = conv_to<fmat>::from(reference_coordinates.rows(0, 2));
      ptr = sub_coor.memptr();
      memcpy(tail_vec.data(), ptr, sizeof(float) * tail_vec.size());

      switch (layout_alg)
      {
      case TUMAP_LAYOUT:
        result = optimize_layout_tumap(
            head_vec, tail_vec, positive_head, positive_tail, n_epochs, nV,
            epochs_per_sample, LEARNING_RATE, NEGATIVE_SAMPLE_RATE, thread_no, 1,
            true, seed);
        break;
      case UMAP_LAYOUT:
        result = optimize_layout_umap(
            head_vec, tail_vec, positive_head, positive_tail, n_epochs, nV,
            epochs_per_sample, a_param, b_param, GAMMA, LEARNING_RATE,
            NEGATIVE_SAMPLE_RATE, false, thread_no, 1, false, seed);
        break;
      default:
        result = optimize_layout_tumap(
            head_vec, tail_vec, positive_head, positive_tail, n_epochs, nV,
            epochs_per_sample, LEARNING_RATE, NEGATIVE_SAMPLE_RATE, thread_no, 1,
            false, seed);
        break;
      }

      fmat coordinates_3D_float(result.data(), 3, nV);
      mat coordinates_3D = conv_to<mat>::from(coordinates_3D_float);
      coordinates_3D = trans(coordinates_3D);

      stdout_printf("done\n"); FLUSH;

      /****************************
     *  Now compute node colors *
     ***************************/
      stdout_printf("\tComputing node colors ... ");

      mat U;
      vec s;
      mat V;
      svd_econ(U, s, V, coordinates_3D);
      mat Z = zscore(U);

      vec a = 75 * Z.col(0);
      vec b = 75 * Z.col(1);

      vec L = Z.col(2);
      L = 25.0 + 70.0 * (L - min(L)) / (max(L) - min(L));

      double r_channel, g_channel, b_channel;
      mat RGB_colors = zeros(nV, 3);
      for (int i = 0; i < nV; i++)
      {
        Lab2Rgb(&r_channel, &g_channel, &b_channel, L(i), a(i), b(i));

        RGB_colors(i, 0) = min(1.0, max(0.0, r_channel));
        RGB_colors(i, 1) = min(1.0, max(0.0, g_channel));
        RGB_colors(i, 2) = min(1.0, max(0.0, b_channel));
      }

      stdout_printf("done\n");
      FLUSH;
      res(1) = coordinates_3D;
      res(2) = RGB_colors; //RGB_colors;
    }
    return res;
  }

  field<mat> layoutNetwork(sp_mat &G, mat initial_position, string algorithm, int compactness_level,
                           unsigned int n_epochs, int thread_no,
                           int seed)
  {
    field<mat> res;
    if (algorithm == "tumap")
    {
      res = layoutNetwork_umap(G, initial_position, compactness_level, n_epochs, TUMAP_LAYOUT, thread_no, seed);
    }
    else if (algorithm == "umap")
    {
      res = layoutNetwork_umap(G, initial_position, compactness_level, n_epochs, UMAP_LAYOUT, thread_no, seed);
    }
    else if (algorithm == "forced_atlas")
    {
      res = layoutNetwork_forced_atlas(G, initial_position, compactness_level, n_epochs, thread_no, seed);
    }
    else
    { // Default to TUMAP
      res = layoutNetwork_umap(G, initial_position, compactness_level, n_epochs, TUMAP_LAYOUT, thread_no, seed);
    }

    return (res);
  }

} // namespace ACTIONet