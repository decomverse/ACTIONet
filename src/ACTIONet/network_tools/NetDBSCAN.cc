#include "ACTIONet.h"
#include <queue>

#define UNDEFINED -1
#define NOISE 0

namespace ACTIONet {

// -1: Undefined, 0: Noise, 1...: Cluster IDs
vec NetDBSCAN(sp_mat &G, int minPts, double eps = 0.5,
              double alpha_val = 0.85) {
  int n = G.n_rows;
  G.diag().zeros(); // Remove self-loops

  // Process vertices in a particular order to give higher priority to the
  // inner, the most well-connected vertices
  printf("Ordering vertices ... ");
  vec cn = conv_to<vec>::from(compute_core_number(G));
  sp_mat X0(G.n_rows, 1);
  X0.col(0) = cn;
  mat pr = compute_network_diffusion(G, X0, 1, alpha_val, 3);
  uvec perm = sort_index(pr.col(0), "descend");
  printf("done\n");

  // Main body of the Net-DBSCAN
  printf("Clustering nodes\n");
  unsigned int C = 0;
  vec clusters = UNDEFINED * ones(n);
  vec visited = zeros(n);

  vector<int> N;
  N.reserve(n); // Allocate memory for the maximum size
  for (int i = 0; i < n; i++) {

    int v = perm[i]; // Visit vertices in order of connectivity
    if (visited[v])
      continue;     // Skip previously visited vertices
    visited[v] = 1; // Mark vertex as visited

    N.clear(); // Flush previous content
    int counts = 0;
    for (sp_mat::col_iterator col_it = G.begin_col(v); col_it != G.end_col(v);
         ++col_it) {
      counts++;
      int u = col_it.row();
      if ((eps < (*col_it)) &&
          ((clusters[u] == UNDEFINED) || (clusters[u] == NOISE)))
        N.push_back(u);
    }

    if (N.size() < minPts) { // Mark it as "noise"
      clusters[v] = NOISE;
      continue;
    }

    clusters[v] = ++C;

    // Mark all neighbors as visited and add them to the seed "queue"
    vector<int> S;
    while (!N.empty()) {
      int u = N.back();
      if ((clusters[u] == UNDEFINED) || (clusters[u] == NOISE)) {
        visited[u] = 1;
        S.push_back(u);
      }
      N.pop_back();
    }
    while (!S.empty()) {
      int u = S.back();
      S.pop_back();

      visited[u] = 1;

      if (clusters[u] == NOISE) {
        clusters[u] = C;
        continue;
      }
      clusters[u] = C;

      N.clear(); // Flush previous content
      for (sp_mat::col_iterator col_it = G.begin_col(u); col_it != G.end_col(u);
           ++col_it) {
        int w = col_it.row();
        if ((eps < (*col_it)) &&
            ((clusters[w] == UNDEFINED) || (clusters[w] == NOISE)))
          N.push_back(w);
      }
      if (N.size() >= minPts) {
        while (!N.empty()) {
          int w = N.back();
          if ((clusters[w] == UNDEFINED) || (clusters[w] == NOISE)) {
            visited[w] = 1;
            S.push_back(w);
          }
          N.pop_back();
        }
      }
    }
  }

  return (clusters);
}
} // namespace ACTIONet
