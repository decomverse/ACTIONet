#ifndef ACTIONETCORE_H
#define ACTIONETCORE_H

#include <algorithm>
#include <cmath>
#include <ctime>
#include <errno.h>
#include <getopt.h>
#include <limits.h>
#include <list>
#include <map>
#include <math.h>
#include <set>
#include <stack>
#include <stdio.h>
#include <string.h>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <unordered_map>
#include <vector>

#include <colorspace.h>
#include <gradient.h>
#include <sampler.h>
#include <tauprng.h>

#include <Rcpp.h>
#include <RcppArmadillo.h>

using namespace arma;
using namespace std;

#define TUMAP_LAYOUT 0
#define UMAP_LAYOUT 1
#define GRAPHVIS_LAYOUT 2

#define ACTIONet_AND 1
#define ACTIONet_OR 2

struct Trace {
  vec parents;
  vec dist;
  vec visited_vertices;
};

Trace dijkstra(sp_mat Dist, int src);

mat zscore(mat A);

const double UMAP_A[101] = {
    1.93280839781719,  1.89560586588002,  1.85873666431227,  1.82221007490834,
    1.78603612060048,  1.75022496320214,  1.71478579945151,  1.67972997626197,
    1.64506544270902,  1.610800661285,    1.57694346052399,  1.54350101780511,
    1.51047986323257,  1.47788588612333,  1.44572435168023,  1.41399925414561,
    1.38271638006498,  1.35187804260518,  1.3214872860387,   1.29154663185922,
    1.26205810311418,  1.23302325071067,  1.20444317424075,  1.17631854866857,
    1.14864964274379,  1.12143634262879,  1.09467817152021,  1.0683743100033,
    1.04252361298475,  1.01712481754341,  0.992175611624647, 0.967674513244996,
    0.943619207179927, 0.920007077834315, 0.896835219021839, 0.874100443595699,
    0.851800999392949, 0.829931994792615, 0.808490430178554, 0.787472613514984,
    0.766873638278737, 0.746690990400437, 0.726919886947928, 0.707556026044195,
    0.688594985599233, 0.670032232635194, 0.651864066568649, 0.634084192553475,
    0.616688494561969, 0.599672088669339, 0.583030020204371, 0.5667572718654,
    0.550848768322639, 0.535299383967892, 0.520103947257001, 0.505257246260431,
    0.490754031684977, 0.476589022213249, 0.46275690208242,  0.449252325341552,
    0.436069912245555, 0.423205974605747, 0.4106531652521,   0.39840668039948,
    0.386461380891047, 0.374811984314975, 0.363453224264704, 0.352379851902848,
    0.341586644916259, 0.331068403184832, 0.320819956874279, 0.31083616902857,
    0.301110995958752, 0.291641183389757, 0.282420831386121, 0.273444955588216,
    0.264708614833586, 0.256206914916444, 0.247935008593902, 0.239888099677924,
    0.232061441819675, 0.224450342118235, 0.217050162160312, 0.209856317524031,
    0.202864281204524, 0.196069583611474, 0.189467814398248, 0.183054621446351,
    0.176825713015038, 0.17077685928726,  0.164903890637922, 0.159202699934773,
    0.153669242222215, 0.148299535941784, 0.143089661250278, 0.138035764053223,
    0.133134049958711, 0.12838079222654,  0.123772324007265, 0.119305671122251,
    0.114976081494676};
const double UMAP_B[101] = {
    0.790494973419029, 0.80063784415826,  0.810876441425738, 0.821199202674006,
    0.831595366275022, 0.84205539236769,  0.852571713401325, 0.863135518043442,
    0.873741680140683, 0.884384956993888, 0.895060878257082, 0.905765637284042,
    0.916495998501859, 0.927249214280422, 0.938022954467018, 0.948815759038301,
    0.95962499558526,  0.970449732070657, 0.981288783823989, 0.992141168965973,
    1.00300608092206,  1.01388286515112,  1.02477099750548,  1.03567006898871,
    1.04657977025277,  1.05749987674998,  1.06843023939592,  1.07937077470387,
    1.09032145585694,  1.10128169075827,  1.11225322117536,  1.12323470900213,
    1.13422639755358,  1.14522861434516,  1.15624176559097,  1.16726633179917,
    1.17830241385901,  1.18934945144456,  1.20040819996369,  1.21147891097075,
    1.22256381651844,  1.23366041866219,  1.24477022428392,  1.2558936051142,
    1.26703094885274,  1.27818265467871,  1.28934756395537,  1.30052872175886,
    1.31172539107843,  1.32293800168803,  1.3341669930459,   1.34541281413396,
    1.35667592718974,  1.36795680610473,  1.37925594017143,  1.39057383474783,
    1.40191101858967,  1.41326804557094,  1.42464550789942,  1.436044048272,
    1.44746436980037,  1.45890393087319,  1.47036701291879,  1.48185337703821,
    1.49336326709497,  1.50489726618312,  1.51645596605121,  1.52803997486173,
    1.53964990048402,  1.55128637349183,  1.56295003156298,  1.57464152150044,
    1.58636409305622,  1.59811350189048,  1.60989278253114,  1.62170263415549,
    1.63354377154668,  1.64541692037945,  1.65732282325244,  1.66926223230814,
    1.68123591907029,  1.69324466615879,  1.70528927262371,  1.71737055545595,
    1.72948934595558,  1.74164649289645,  1.75384285823827,  1.76607932576738,
    1.77835679827623,  1.79067619009556,  1.80303844043406,  1.81544450541945,
    1.82789536263139,  1.84039200538657,  1.85293545544251,  1.86552674229068,
    1.87816693701183,  1.89085711093115,  1.90359837758981,  1.91638829237987,
    1.92923479503841};

#define NEGATIVE_SAMPLE_RATE 5.0
#define LEARNING_RATE 1.0
#define UMAP_SEED 0
#define GAMMA 1.0

vector<double> optimize_layout(
    const apumap_gradient &gradient, vector<double> &head_embedding,
    vector<double> &tail_embedding, const vector<unsigned int> &positive_head,
    const vector<unsigned int> &positive_tail, unsigned int n_epochs,
    unsigned int n_vertices, const vector<double> &epochs_per_sample,
    double initial_alpha, double negative_sample_rate, unsigned int seed);

namespace ACTIONetcore {

struct multilevel_archetypal_decomposition {
  vec landmark_cells;  // Closest cell to each archetype. It's -1, if archetype
                       // is close to more than one vertex (i.e., it's hub or
                       // nonspecific)
  uvec selected_archs; // If hub removal requested, this will hold the indices
                       // of retained archetypes
  mat C_stacked; // Stacking of C matrices, after potentially removing the hub
                 // archetypes
  mat H_stacked; // Stacking of H matrices, after potentially removing the hub
                 // archetypes
  mat archetype_profile; // gene x # archetype (after potentially removing hub
                         // archetypes). Computed as S*C_stacked
  mat backbone; // archetype x archetype graph (after potentially removing hub
                // archetypes).
};

mat computeFullSim(mat &H, int thread_no);
field<sp_mat> buildAdaptiveACTIONet(mat &H_stacked, double LC, double M,
                                    double ef_construction, double ef,
                                    int thread_no, int sym_method);

field<mat> layoutACTIONet(sp_mat &G, mat &S_r, int compactness_level,
                          unsigned int n_epochs, int thread_no);
mat update_layout_2D(mat &coors, int compactness_level, unsigned int n_epochs,
                     int thread_no);

multilevel_archetypal_decomposition reconstructArchetypes(sp_mat S,
                                                          vector<mat> C_trace,
                                                          vector<mat> H_trace,
                                                          double z_score);

mat assessFeatureSets(sp_mat &S, field<uvec> feature_sets, int rand_perm);

mat phenotypeEnrichment(mat &H_stacked, mat &phenotype_associations,
                        int rand_perm_no);

// Basic network algorithms
mat PR_iter(sp_mat &G, sp_mat &X0, double alpha, int max_it, int thread_no);
vec sweepcut(sp_mat &A, vec s);
mat MWM(mat &G);

} // namespace ACTIONetcore

#endif
