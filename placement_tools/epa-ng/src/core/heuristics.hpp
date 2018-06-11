#pragma once

#include <algorithm>

#ifdef __OMP
#include <omp.h>
#endif

#include "core/Work.hpp"
#include "sample/Sample.hpp"
#include "util/Options.hpp"
#include "set_manipulators.hpp"

static inline size_t get_num_threads(const Options& options)
{
  #ifdef __OMP
  const size_t num_threads  = options.num_threads
                            ? options.num_threads
                            : omp_get_max_threads();
  omp_set_num_threads(num_threads);
  #else
  (void) options;
  const size_t num_threads = 1;
  #endif
  return num_threads;
}

static inline size_t get_thread_id()
{
  #ifdef __OMP
  return omp_get_thread_num();
  #else
  return 0;
  #endif
}

inline Work dynamic_heuristic(Sample<Placement>& sample,
                              const Options& options)
{
  Work result;
  compute_and_set_lwr(sample);

  const auto num_threads = get_num_threads(options);

  std::vector<Work> workvec(num_threads);

  #ifdef __OMP
  #pragma omp parallel for schedule(dynamic)
  #endif
  for (size_t i = 0; i < sample.size(); ++i) {
    auto &pq = sample[i];
    const auto tid = get_thread_id();

    auto end = until_accumulated_reached( pq,
                                          options.prescoring_threshold,
                                          options.filter_min,
                                          options.filter_max);

    for (auto iter = pq.begin(); iter != end; ++iter) {
      workvec[tid].add(iter->branch_id(), pq.sequence_id());
    }
  }
  merge(result, workvec);
  return result;
}

inline Work fixed_heuristic(Sample<Placement>& sample,
                            const Options& options)
{
  Work result;
  compute_and_set_lwr(sample);

  const auto num_threads = get_num_threads(options);

  std::vector<Work> workvec(num_threads);

  #ifdef __OMP
  #pragma omp parallel for schedule(dynamic)
  #endif
  for (size_t i = 0; i < sample.size(); ++i) {
    auto &pq = sample[i];
    const auto tid = get_thread_id();

    auto end = until_top_percent(pq, options.prescoring_threshold);

    for (auto iter = pq.begin(); iter != end; ++iter) {
      workvec[tid].add(iter->branch_id(), pq.sequence_id());
    }
  }

  merge(result, workvec);
  return result;
}

inline Work baseball_heuristic( Sample<Placement>& sample,
                                const Options& options)
{
  Work result;

  const auto num_threads = get_num_threads(options);

  // strike_box: logl delta, keep placements within this many logl units from the best
  const double strike_box = 3;
  // max_strikes: number of additional branches to add after strike box is full
  const size_t max_strikes = 6;
  // max_pitches: absolute maximum of candidates to select
  const size_t max_pitches = 40;
  
  std::vector<Work> workvec(num_threads);
  #ifdef __OMP
  #pragma omp parallel for schedule(dynamic)
  #endif
  for (size_t i = 0; i < sample.size(); ++i) {
    auto &pq = sample[i];
    const auto tid = get_thread_id();

    assert(pq.size());
    // sort placements by likelihood (descending)
    sort_by_logl(pq);
    // keep any placements that are within strike box of the best
    const double best_logl = pq[0].likelihood();
    const double thresh = best_logl - strike_box;
    // get first element not within strike box
    auto keep_iter = std::find_if(pq.begin(), pq.end(),
      [thresh](const auto& p){
        return (p.likelihood() < thresh);
      }
    );

    const auto hits = std::distance(pq.begin(), keep_iter);

    // ensure we keep no more than max_pitches
    size_t to_add = std::min(max_pitches - hits, max_strikes);

    std::advance(keep_iter, to_add);

    for (auto iter = pq.begin(); iter != keep_iter; ++iter) {
      workvec[tid].add(iter->branch_id(), pq.sequence_id());
    }

  }
  merge(result, workvec);
  return result;
}

inline Work apply_heuristic(Sample<Placement>& sample,
                            const Options& options)
{
  if (options.baseball) {
    return baseball_heuristic(sample, options);
  } else if (options.prescoring_by_percentage) {
    return fixed_heuristic(sample, options);
  } else {
    return dynamic_heuristic(sample, options);
  }
}
