/************************************************
** This is a local search solver for Minimum Vertex Cover.
************************************************/

/************************************************
** Date:  2011.7.1
** NuMVC (Two Stage Exchange and Weighting with Forgetting)
** Author: Shaowei Cai, shaowei_cai@126.com
**       School of EECS, Peking University
**       Beijing, China
**
** Date:  2011.10.28
** Modify: Shaowei Cai
** use dynamic memory for v_adj[][] and v_edge[][], tidy codes.
**
** Date: 2017.11.17
** Modify: Felix Moessbauer
** Object oriented design, C++11 support
** Use heap for construction of initial solution
************************************************/

#ifndef NuMVC_INCLUDED
#define NuMVC_INCLUDED

#include <algorithm>
#include <chrono>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <utility>
#include <vector>

#include "internal/indexed_heap.hpp"

namespace libmvc {

class NuMVC {
 public:
  using Edge = std::pair<int, int>;

 private:
  using timepoint_t = std::chrono::time_point<std::chrono::system_clock>;
  using duration_ms = std::chrono::milliseconds;

  using compare_t = std::function<bool(const int &, const int &)>;
  using heap_t = internal::Indexed_Heap<int, compare_t>;

 private:
  static constexpr int try_step = 100000;

  bool verbose = false;

  timepoint_t start, finish;

  /*parameters of algorithm*/
  duration_ms cutoff_time;  // time limit
  long long step;
  int optimal_size;  // terminate the algorithm before step limit if it finds a
                     // vertex cover of optimal_size

  /*parameters of the instance*/
  int v_num;  //|V|: 0...v-1
  int e_num;  //|E|: 0...e-1

  /*structures about edge*/
  std::vector<Edge> edge;
  std::vector<int> edge_weight;
  int default_edge_weight = 1;

  /*structures about vertex*/
  std::vector<int> dscore;  // dscore of v
  std::vector<long long> time_stamp;
  int best_cov_v;

  // from vertex to it's edges and neighbors
  std::vector<int> v_beg_idx;  // v_edges and v_adj is flattened 2-d array,
                               // hence store indices
  std::vector<int> v_edges;    // edges related to v, v_edges[i][k] means vertex
                               // v_i's k_th edge
  std::vector<int>
      v_adj;  // v_adj[v_i][k] = v_j (actually, that is v_i's k_th neighbor)

  /* structures about solution */
  // current candidate solution
  int c_size;                    // cardinality of C
  std::vector<char> v_in_c;      // a flag indicates whether a vertex is in C
  std::vector<int> remove_cand;  // remove candidates, an array consists of only
                                 // vertices in C, not including tabu_remove
  std::vector<int> index_in_remove_cand;
  int remove_cand_size;

  // best solution found
  int best_c_size;
  std::vector<char>
      best_v_in_c;  // a flag indicates whether a vertex is in best solution

  duration_ms best_comp_time;
  long best_step;

  // uncovered edge stack
  std::vector<int> uncov_stack;  // store the uncov edge number
  int uncov_stack_fill_pointer;
  std::vector<int>
      index_in_uncov_stack;  // which position is an edge in the uncov_stack

  // CC and taboo, pack densly to avoid cache misses
  std::vector<bool> conf_change;
  int tabu_remove = 0;

  // smooth
  static constexpr float p_scale = 0.3;  // w=w*p
  int delta_total_weight = 0;
  int ave_weight = 1;

  // random
  std::mt19937 mt_rand;

 public:
  /**
   * Construct solver instance by importing a graph in DIMACS format
   */
  template <typename Is, typename Duration>
  NuMVC(
      /// Input Stream with graph in DIMACS format
      Is &str,
      /// Size of optimal vertex cover (set to 0 if not known)
      int optimal_size,
      /// Stop calculation after this duration (chrono duration)
      Duration cutoff_time,
      /// Print messages during calculation
      bool verbose = false,
      /// seed for random number generator
      unsigned int rnd_seed =
          std::chrono::high_resolution_clock::now().time_since_epoch().count())
      : verbose(verbose),
        cutoff_time(std::chrono::duration_cast<duration_ms>(cutoff_time)),
        optimal_size(optimal_size) {
    set_random_seed(rnd_seed);
    build_instance(str);
    init_sol();
  }

  template <typename Duration>
  NuMVC(
      /// Graph in edge-list format
      const std::vector<std::pair<int, int>> &edges,
      /// Number of vertices in graph
      const int &num_vertices,
      /// Size of optimal vertex cover (set to 0 if not known)
      int optimal_size,
      /// Stop calculation after this duration (chrono duration)
      Duration cutoff_time,
      /// Print messages during calculation
      bool verbose = false,
      /// seed for random number generator
      unsigned int rnd_seed =
          std::chrono::high_resolution_clock::now().time_since_epoch().count())
      : verbose(verbose),
        cutoff_time(std::chrono::duration_cast<duration_ms>(cutoff_time)),
        optimal_size(optimal_size) {
    set_random_seed(rnd_seed);
    build_instance(num_vertices, edges);
    init_sol();
  }

 private:
  void init_internal(int num_vertices, int num_edges) {
    edge.resize(num_edges);
    edge_weight.resize(num_edges, default_edge_weight);
    dscore.resize(num_vertices);
    time_stamp.resize(num_vertices);

    v_in_c.resize(num_vertices);
    remove_cand.resize(num_vertices);
    index_in_remove_cand.resize(num_vertices);

    best_v_in_c.resize(num_vertices);

    // uncovered edge stack
    uncov_stack.resize(num_edges);  // store the uncov edge number
    index_in_uncov_stack.resize(
        num_edges);  // which position is an edge in the uncov_stack

    // CC and taboo
    conf_change.resize(num_vertices, 1);
  }

  inline int degree(int v) const {
    return v_beg_idx[v+1] - v_beg_idx[v];
  }

  inline void update_best_sol() {
    best_v_in_c = v_in_c;
    best_c_size = c_size;
    finish = std::chrono::system_clock::now();
    best_comp_time = std::chrono::duration_cast<duration_ms>(finish - start);
    best_step = step;
  }

  /** choose a vertex u in C with the highest dscore,
   * breaking ties in favor of the oldest one;
   */
  inline void update_best_cov_v() {
    best_cov_v  = remove_cand[0];
    auto dscore_best = dscore[best_cov_v];

    for (int i = 0; i < remove_cand_size - 1; ++i) {
      int v = remove_cand[i];
      if (dscore[v] < dscore_best) continue;
      if (dscore[v] > dscore_best){
        if (v != tabu_remove){
          best_cov_v  = v;
          dscore_best = dscore[v];
        }
      } else if (time_stamp[v] < time_stamp[best_cov_v]){
        if (v != tabu_remove){
          best_cov_v = v;
          dscore_best = dscore[v];
        }
      }
    }
  }

  template <typename Is>
  void build_instance(Is &str) {
    char line[1024];
    char tempstr1[10];
    char tempstr2[10];
    int e;

    char tmp;
    int v1, v2;

    /*** build problem data structures of the instance ***/
    str.getline(line, 1024);
    while (line[0] != 'p') str.getline(line, 1024);
    sscanf(line, "%s %s %d %d", tempstr1, tempstr2, &v_num, &e_num);
    init_internal(v_num, e_num);

    std::vector<int> v_degree(v_num);
    for (e = 0; e < e_num; ++e) {
      str >> tmp >> v1 >> v2;
      // start counting at zero
      --v1;
      --v2;
      ++(v_degree[v1]);
      ++(v_degree[v2]);

      edge[e] = {v1, v2};
    }
    update_instance_internal(v_degree);
  }

  void build_instance(const int &num_vertices,
                      const std::vector<std::pair<int, int>> &edges) {
    v_num = num_vertices;
    e_num = edges.size();

    init_internal(v_num, e_num);

    std::vector<int> v_degree(v_num);
    for (unsigned int e = 0; e < edges.size(); ++e) {
      ++(v_degree[edges[e].first]);
      ++(v_degree[edges[e].second]);
    }
    edge = edges;
    update_instance_internal(v_degree);
  }

  /**
   * builds internal data structures for processing this instance.
   * has to be called after loading a new instance using
   * build_instance
   */
  void update_instance_internal(const std::vector<int> & v_degree) {
    // indices are the partial sums
    // first offset is 0, last is partial including last element
    v_beg_idx.reserve(e_num);
    v_beg_idx.push_back(0);
    std::partial_sum(v_degree.begin(), v_degree.end(),
                     std::back_inserter(v_beg_idx));

    v_edges.resize(v_beg_idx.back());
    v_adj.resize(v_beg_idx.back());

    int v1, v2;
    std::vector<int> v_degree_tmp(v_num);
    for (int e = 0; e < e_num; ++e) {
      v1 = edge[e].first;
      v2 = edge[e].second;

      v_edges[v_beg_idx[v1] + v_degree_tmp[v1]] = e;
      v_edges[v_beg_idx[v2] + v_degree_tmp[v2]] = e;

      v_adj[v_beg_idx[v1] + v_degree_tmp[v1]] = v2;
      v_adj[v_beg_idx[v2] + v_degree_tmp[v2]] = v1;

      ++(v_degree_tmp[v1]);
      ++(v_degree_tmp[v2]);
    }
  }

  inline void reset_remove_cand() {
    int j = 0;
    for (int v = 0; v < v_num; ++v) {
      if (v_in_c[v]) {  // && v!=tabu_remove)
        remove_cand[j] = v;
        index_in_remove_cand[v] = j;
        ++j;
      } else
        index_in_remove_cand[v] = 0;
    }
    remove_cand_size = j;
  }

  // kick out the worst vertex in current cover
  void update_target_size() {
    --c_size;

    int max_improvement = std::numeric_limits<int>::min();
    int max_vertex = 0;  // vertex with the highest improvement in C

    for (int v = 0; v < v_num; ++v) {
      if (!v_in_c[v]) continue;
      if (dscore[v] > max_improvement) {
        max_improvement = dscore[v];
        max_vertex = v;
      }
    }
    remove(max_vertex);

    reset_remove_cand();
  }

  inline void uncover(int e) {
    index_in_uncov_stack[e] = uncov_stack_fill_pointer;
    uncov_stack[uncov_stack_fill_pointer++] = e;
  }

  inline void cover(int e) {
    // since the edge is satisfied, its position can be reused to store the
    // last_uncov_edge
    int last_uncov_edge = uncov_stack[--uncov_stack_fill_pointer];
    int index = index_in_uncov_stack[e];
    uncov_stack[index] = last_uncov_edge;
    index_in_uncov_stack[last_uncov_edge] = index;
  }

  void init_sol() {
    start = std::chrono::system_clock::now();

    auto dscore_cmp = [&dscore = dscore](const int &a, const int &b) {
      return (dscore[a] < dscore[b]);
    };
    heap_t v_heap(dscore_cmp);
    v_heap.resize(v_num);

    /*** build solution data structures of the instance ***/
    // init vertex cover
    // conf_change = 1 (already set in resize call)
    for (int e = 0; e < e_num; ++e) {
      const auto weight = edge_weight[e];
      dscore[edge[e].first] += weight;
      dscore[edge[e].second] += weight;
    }

      // init uncovered edge stack and cover_vertrex_count_of_edge array

#if 0
        uncov_stack_fill_pointer = 0;
        for (int e=0; e<e_num; ++e)
            uncover(e);
#else
    // tweak to uncover all
    std::iota(index_in_uncov_stack.begin(), index_in_uncov_stack.end(), 0);
    memcpy(uncov_stack.data(), index_in_uncov_stack.data(), e_num);
    uncov_stack_fill_pointer = e_num;
#endif

    for (int v = 0; v < v_num; ++v) {
      v_heap.push(v);
    }

    int i = 0;
    while (uncov_stack_fill_pointer > 0) {
      int best_v = v_heap.top();
      v_heap.pop();
      if (dscore[best_v] > 0) {
        add_init(best_v, v_heap);
        ++i;
      }
    }

    c_size = i;
    update_best_sol();

    finish = std::chrono::system_clock::now();
    auto init_sol_time =
        std::chrono::duration_cast<duration_ms>(finish - start);
    if (verbose) {
      std::cout << "Initial solution size: " << c_size << std::endl;
      std::cout << "Initial solution time: " << init_sol_time.count() << "ms"
                << std::endl;
    }

    reset_remove_cand();
    update_best_cov_v();
  }

  // add a vertex to current cover
  void add(int v) {
    v_in_c[v] = true;
    dscore[v] *= -1;

    int edge_count = degree(v);
    int idx_v = v_beg_idx[v];

    for (int i = 0; i < edge_count; ++i) {
      int e = v_edges[idx_v+i];  // v's i'th edge
      int n = v_adj[idx_v+i];    // v's i'th neighbor

      if (!v_in_c[n]) {  // this adj isn't in cover set
        dscore[n] -= edge_weight[e];
        conf_change[n] = 1;

        cover(e);
      } else {
        dscore[n] += edge_weight[e];
      }
    }
  }

  void add_init(int v, heap_t &v_heap) {
    v_in_c[v] = true;
    dscore[v] *= -1;

    int edge_count = degree(v);
    int idx_v = v_beg_idx[v];

    for (int i = 0; i < edge_count; ++i) {
      int e = v_edges[idx_v + i];  // v's i'th edge
      int n = v_adj[idx_v + i];    // v's i'th neighbor

      if (!v_in_c[n]) {  // this adj isn't in cover set
        bool inheap = v_heap.count(n) > 0;
        if (inheap) {
          v_heap.erase(v_heap[n]);
        }
        dscore[n] -= edge_weight[e];
        conf_change[n] = 1;
        if (inheap) {
          v_heap.push(n);
        }
        cover(e);
      } else {
        dscore[n] += edge_weight[e];
      }
    }
  }

  void remove(int v) {
    v_in_c[v] = false;
    dscore[v] *= -1;
    conf_change[v] = 0;

    int edge_count = degree(v);
    int idx_v = v_beg_idx[v];

    for (int i = 0; i < edge_count; ++i) {
      int e = v_edges[idx_v + i];
      int n = v_adj[idx_v + i];

      if (!v_in_c[n]) {  // this adj isn't in cover set
        dscore[n] += edge_weight[e];
        conf_change[n] = 1;

        uncover(e);
      } else {
        dscore[n] -= edge_weight[e];
      }
    }
  }

  void forget_edge_weights() {
    int new_total_weight = 0;

    std::fill(dscore.begin(), dscore.end(), 0);

    // scale_ave=ave_weight*q_scale;
    for (int e = 0; e < e_num; ++e) {
      edge_weight[e] = edge_weight[e] * p_scale;

      new_total_weight += edge_weight[e];

      // update dscore
      if (!(v_in_c[edge[e].first] || v_in_c[edge[e].second])) {
        dscore[edge[e].first] += edge_weight[e];
        dscore[edge[e].second] += edge_weight[e];
      } else if (v_in_c[edge[e].first] != v_in_c[edge[e].second]) {
        if (v_in_c[edge[e].first])
          dscore[edge[e].first] -= edge_weight[e];
        else
          dscore[edge[e].second] -= edge_weight[e];
      }
    }
    ave_weight = new_total_weight / e_num;
  }

  void update_edge_weight() {
    for (int i = 0; i < uncov_stack_fill_pointer; ++i) {
      int e = uncov_stack[i];

      edge_weight[e] += 1;
      dscore[edge[e].first] += 1;
      dscore[edge[e].second] += 1;
    }

    delta_total_weight += uncov_stack_fill_pointer;
    if (delta_total_weight >= e_num) {
      ave_weight += 1;
      delta_total_weight -= e_num;
    }

    // smooth weights
    int threshold = v_num / 2;
    if (ave_weight >= threshold) {
      forget_edge_weights();
    }
  }

 public:
  /**
   * calculate minimum vertex cover
   */
  void cover_LS() { cover_LS(nullptr); }

  /**
   * calculate minimum vertex cover and call callback after
   * each iteration. If callback returns true, stop calculation.
   */
  void cover_LS(
      const std::function<bool(const NuMVC &, bool)> &callback_on_update) {

    // initialize random pool
    step = 1;
    start = std::chrono::system_clock::now();

    while (true) {
      /* ### if there is no uncovered edge ### */
      if (uncov_stack_fill_pointer <= 0) {
        update_best_sol();  // C* := C
        if (callback_on_update != nullptr && callback_on_update(*this, true)) {
          return;
        }

        if (c_size <= optimal_size) return;

        update_target_size();  // remove a vertex with the highest dscore from
                               // C;
        continue;
      }

      /* monitor status */
      if (step % try_step == 0) {
        finish = std::chrono::system_clock::now();
        auto elapsed_ms =
            std::chrono::duration_cast<duration_ms>(finish - start);
        // call monitor to support external cancelling
        if (callback_on_update != nullptr && callback_on_update(*this, false)) {
          return;
        }
        if (elapsed_ms >= cutoff_time) {
          if (verbose) {
            std::cout << "Time limit reached:" << elapsed_ms.count() << "ms"
                      << std::endl;
          }
          return;
        }
      }

      /* choose a vertex u in C with the highest dscore,
        breaking ties in favor of the oldest one; */
      update_best_cov_v();

      /* C := C\{u}, confChange(u) := 0 and confChange(z) := 1 for each z in
       * N(u); */
      remove(best_cov_v);

      /* choose an uncovered edge e randomly; */

      auto next_rnd = mt_rand();
      int e = uncov_stack[next_rnd % uncov_stack_fill_pointer];

      /* choose a vertex v in e such that confChange(v) = 1 with higher dscore,
        breaking ties in favor of the older one; */
      int v1 = edge[e].first;
      int v2 = edge[e].second;
      int best_add_v;

      if (conf_change[v1] == 0)
        best_add_v = v2;
      else if (conf_change[v2] == 0)
        best_add_v = v1;
      else {
        if (dscore[v1] > dscore[v2] ||
            (dscore[v1] == dscore[v2] && time_stamp[v1] < time_stamp[v2]))
          best_add_v = v1;
        else
          best_add_v = v2;
      }

      /* C := C plus {v}, confChange(z) := 1 for each z in N(v); */
      add(best_add_v);

      // update remove candidate index. speed optimization.
      int index = index_in_remove_cand[best_cov_v];
      index_in_remove_cand[best_cov_v] = 0;

      remove_cand[index] = best_add_v;
      index_in_remove_cand[best_add_v] = index;

      // update timestamp used to determine the age of each virtex
      time_stamp[best_add_v] = step;
      time_stamp[best_cov_v] = step;

      tabu_remove = best_add_v;

      /* w(e) := w(e) + 1 for each uncovered edge e; */
      /* if w >= y then w(e) := [p*w(e)] for each edge e; */
      update_edge_weight();

      ++step;
    }
    return;
  }

  /**
   * Check if the calculated solution is a valid vertex cover
   */
  bool check_solution() const {
    int e;
    for (e = 0; e < e_num; ++e) {
      if (!(best_v_in_c[edge[e].first] || best_v_in_c[edge[e].second])) {
        std::cout << "uncovered edge " << e << std::endl;
        return false;
      }
    }

    return true;
  }

  /**
   * set the maximum duration after which the solver terminates
   *
   * \param d any \c std::chrono::duration
   */
  template <typename Duration>
  void set_cutoff_time(Duration d) {
    cutoff_time = std::chrono::duration_cast<duration_ms>(d);
  }

  /**
   * Set the size of the vertex cover at which the algorithm should stop
   */
  void set_optimal_size(int size) { optimal_size = size; }

  /**
   * set the seed used for the pseudo random number generator
   */
  void set_random_seed(unsigned int seed) { mt_rand.seed(seed); }

  /**
   * returns the current instance as a pair consisting of the
   * number of vertices and a vector of edges
   */
  std::pair<int, std::vector<Edge>> get_instance_as_edgelist() const {
    return std::make_pair(v_num, edge);
  }

  /**
   * return vertex indices of current best vertex cover
   */
  std::vector<int> get_cover(bool bias_by_one = true) const {
    std::vector<int> cover;
    for (int i = 0; i < v_num; i++) {
      if (best_v_in_c[i]) {
        cover.push_back(i + bias_by_one);
      }
    }
    return cover;
  }

  /**
   * returns a vector of flags, where a true-flag at position i denots
   * that vertex i is covered
   */
  std::vector<char> get_cover_as_flaglist() const { return best_v_in_c; }

  /**
   * return vertex indices of current best independent set
   */
  std::vector<int> get_independent_set(bool bias_by_one = true) const {
    std::vector<int> iset;
    for (int i = 0; i < v_num; i++) {
      if (!best_v_in_c[i]) {
        iset.push_back(i + bias_by_one);
      }
    }
    return iset;
  }

  /**
   * Number of vertices
   */
  inline int get_vertex_count() const { return v_num; }

  /**
   * Number of edges
   */
  inline int get_edge_count() const { return e_num; }

  /**
   * Size of current best vertex cover
   */
  inline int get_best_cover_size() const { return best_c_size; }

  /**
   * Tries necessary for current best vertex cover
   */
  inline long get_best_step() const { return best_step; }

  /**
   * duration for calculating current best vertex cover
   */
  inline std::chrono::milliseconds get_best_duration() const {
    return std::chrono::duration_cast<std::chrono::milliseconds>(
        best_comp_time);
  }

  /**
   * total duration since start of calculation
   */
  inline std::chrono::milliseconds get_total_duration() const {
    return std::chrono::duration_cast<std::chrono::milliseconds>(
        (std::chrono::system_clock::now() - start));
  }

  /**
   * Print statistics during calculation
   */
  static bool default_stats_printer(const NuMVC &solver,
                                    bool better_cover_found) {
    if (better_cover_found) {
      auto time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
          solver.get_best_duration());
      std::cout << "Better MVC found.\tSize: " << solver.get_best_cover_size()
                << "\tTime: " << std::fixed << std::setw(4)
                << std::setprecision(4) << time_ms.count() << "ms" << std::endl;
    }
    return false;
  }
};

}  // namespace libmvc

#endif
