/************************************************
** This is a local search solver adapter
************************************************/

#ifndef PARALLEL_SOLVER_ADAPTER_INCLUDED
#define PARALLEL_SOLVER_ADAPTER_INCLUDED

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <numeric>
#include <thread>
#include <utility>
#include <vector>
#include <memory>

namespace libmvc {

/**
 * Adapter to parallelize any MVC solver compatible with the LibMVC interface.
 *
 * The parallel solver itself fulfills the LibMVC interface as well.
 */
template <typename SOLVER>
class ParallelSolverAdapter {
 public:
  using Edge = std::pair<int, int>;

 private:
  using timepoint_t = std::chrono::time_point<std::chrono::system_clock>;
  using duration_ms = std::chrono::milliseconds;

 private:
  class SolverState {
   public:
    explicit SolverState(int opt_size){ optimal_cover = opt_size; }

    std::mutex mon_mx;
    std::atomic<int> best_cover{std::numeric_limits<int>::max()};
    std::atomic<bool> stop_solver{false};
    int best_solver = 0;
    int optimal_cover = 0;
  };

  bool verbose = false;
  // time limit
  duration_ms cutoff_time;
  // invariant master solver
  SOLVER master;
  // dynamic local solvers
  std::vector<std::unique_ptr<SOLVER>> solvers;
  // state object to coodinate solving
  SolverState global_state;

  unsigned long random_base_seed;
  unsigned int  num_instances = std::thread::hardware_concurrency();

 public:
  /**
   * Construct solver instance by importing a graph in DIMACS format
   */
  template <typename Is, typename Duration>
  ParallelSolverAdapter(
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
        master(str, optimal_size, cutoff_time, verbose),
        global_state(optimal_size),
        random_base_seed(rnd_seed) { }

  template <typename Duration>
  ParallelSolverAdapter(
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
        master(edges, num_vertices, optimal_size, cutoff_time),
        global_state(optimal_size),
        random_base_seed(rnd_seed) { }

 private:
  static bool monitor(
      const SOLVER & /* solver */,
      bool better_cover_found,
      unsigned int tid,
      ParallelSolverAdapter *self,
      std::function<bool(const ParallelSolverAdapter &, bool)> printer)
  {
    auto cover_size = self->solvers[tid]->get_best_cover_size();
    auto &state = self->global_state;


    // only lock if actual change present
    if (better_cover_found) {
      if (cover_size < state.best_cover) {
        {
          std::lock_guard<std::mutex> lock(state.mon_mx);
          state.best_cover = cover_size;
          state.best_solver = tid;
        }
        if (printer != nullptr && printer(*self, true)) {
          state.stop_solver = true;
          return true;
        }
      }
    } else {
      if (printer != nullptr && printer(*self, false)) {
        state.stop_solver = true;
        return true;
      }
    }
    // false -> continue if best cover not found yet
    // true  -> stop calculation
    return (state.best_cover == state.optimal_cover);
  }

  static void start_solver(
      ParallelSolverAdapter *self, unsigned int seed, unsigned int tid,
      std::function<bool(const ParallelSolverAdapter &, bool)> printer) {
    // copy master solver to local memory (first touch policy)
    self->solvers[tid] = std::make_unique<SOLVER>(self->master);
    self->solvers[tid]->set_random_seed(seed + tid);
    self->solvers[tid]->cover_LS(std::bind(monitor, std::placeholders::_1,
                                          std::placeholders::_2, tid, self,
                                          printer));
    if (self->verbose) {
      std::cout << "-- " << tid << " terminated" << std::endl;
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
  void cover_LS(const std::function<bool(const ParallelSolverAdapter &, bool)>
                    callback_on_update) {
    std::vector<std::thread> workers;
    solvers.resize(num_instances);
    global_state.stop_solver = false;

    if (verbose) {
      std::cout << "Using " << num_instances << " parallel instances"
                << std::endl;
    }

    // start threads
    for (unsigned int i = 1; i < num_instances; ++i) {
      workers.emplace_back(start_solver, this, random_base_seed, i,
                           callback_on_update);
    }
    start_solver(this, random_base_seed, 0, callback_on_update);

    for (auto &w : workers) {
      w.join();
    }
  }

  /**
   * Set the number of parallel solver instances to be started.
   */
  void set_concurrency(unsigned int instances) { num_instances = instances; }

  /**
   * Check if the calculated solution is a valid vertex cover
   */
  bool check_solution() const {
    auto best = global_state.best_solver;
    return solvers[best]->check_solution();
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
  void set_optimal_size(int size) { global_state.optimal_cover = size; }

  /**
   * set the base-seed used for the pseudo random number generator
   */
  void set_random_seed(unsigned int seed) { random_base_seed = seed; }

  /**
   * returns the current instance as a pair consisting of the
   * number of vertices and a vector of edges
   */
  std::pair<int, std::vector<Edge>> get_instance_as_edgelist() const {
    return master.get_instance_as_edgelist();
  }

  /**
   * returns the number of parallel solvers to be started
   */
  unsigned int get_concurrency() const { return num_instances; }

  /**
   * return vertex indices of current best vertex cover
   */
  std::vector<int> get_cover(bool bias_by_one = true) const {
    auto best = global_state.best_solver;
    return solvers[best]->get_cover(bias_by_one);
  }

  /**
   * returns a vector of flags, where a true-flag at position i denots
   * that vertex i is covered
   */
  std::vector<char> get_cover_as_flaglist() const {
    auto best = global_state.best_solver;
    return solvers[best]->get_cover_as_flaglist();
  }

  /**
   * return vertex indices of current best independent set
   */
  std::vector<int> get_independent_set(bool bias_by_one = true) const {
    auto best = global_state.best_solver;
    return solvers[best]->get_independent_set(bias_by_one);
  }

  /**
   * Number of vertices
   */
  inline int get_vertex_count() const { return master.get_vertex_count(); }

  /**
   * Number of edges
   */
  inline int get_edge_count() const { return master.get_edge_count(); }

  /**
   * Size of current best vertex cover
   */
  inline int get_best_cover_size() const {
    auto best = global_state.best_solver;
    return solvers[best]->get_best_cover_size();
  }

  /**
   * Tries necessary for current best vertex cover
   */
  inline long get_best_step() const {
    std::vector<long> all_best_steps(solvers.size());
    std::transform(solvers.begin(), solvers.end(), all_best_steps.begin(),
                   [](const auto &s) { return s->get_best_step(); });
    return std::accumulate(all_best_steps.begin(), all_best_steps.end(), 0);
  }

  /**
   * duration for calculating current best vertex cover
   */
  inline std::chrono::milliseconds get_best_duration() const {
    auto best = global_state.best_solver;
    return solvers[best]->get_best_duration();
  }

  /**
   * total duration since start of calculation
   */
  inline std::chrono::milliseconds get_total_duration() const {
    return master.get_total_duration();
  }

  /**
   * Print statistics during calculation
   */
  static bool default_stats_printer(const ParallelSolverAdapter &self,
                                    bool better_cover_found) {
    if (better_cover_found) {
      const int best_solver_id = self.global_state.best_solver;
      const auto &best_solver = *(self.solvers[best_solver_id]);
      auto time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
          best_solver.get_best_duration());
      std::cout << std::setw(2) << best_solver_id
                << ": Better MVC found.\tSize: "
                << best_solver.get_best_cover_size() << "\tTime: " << std::fixed
                << std::setw(4) << std::setprecision(4) << time_ms.count()
                << "ms" << std::endl;
    }
    return false;
  }
};

}  // namespace libmvc

#endif
