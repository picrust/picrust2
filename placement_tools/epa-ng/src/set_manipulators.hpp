#pragma once

#include <vector>
#include <unordered_map>
#include <limits>

#include "sample/Sample.hpp"
#include "util/Timer.hpp"
#include "util/Options.hpp"
#include "util/logging.hpp"
#include "seq/MSA.hpp"
#include "core/Work.hpp"

using pq_citer_t  = PQuery<Placement>::const_iterator;
using pq_iter_t   = PQuery<Placement>::iterator;

void merge(Work& dest, const Work& src);
void merge(Timer<>& dest, const Timer<>& src);

void sort_by_lwr(PQuery<Placement>& pq);
void sort_by_logl(PQuery<Placement>& pq);
void compute_and_set_lwr(Sample<Placement>& sample);
pq_iter_t until_top_percent( PQuery<Placement>& pq,
                              const double x);
void discard_bottom_x_percent(Sample<Placement>& sample, const double x);
void discard_by_support_threshold(Sample<Placement>& sample,
                                  const double thresh,
                                  const size_t min=1,
                                  const size_t max=std::numeric_limits<size_t>::max());

pq_iter_t until_accumulated_reached( PQuery<Placement>& pq,
                                      const double thresh,
                                      const size_t min,
                                      const size_t max);
void discard_by_accumulated_threshold(Sample<Placement>& sample,
                                      const double thresh,
                                      const size_t min=1,
                                      const size_t max=std::numeric_limits<size_t>::max());
void filter(Sample<Placement>& sample, const Options& options);
void find_collapse_equal_sequences(MSA& msa);

/**
 * collapses PQuerys with the same ID inside a Sample into one
 */
template<class T>
void collapse(Sample<T>& sample)
{
  const auto invalid = std::numeric_limits<
    typename PQuery<T>::seqid_type>::max();

  std::unordered_map< size_t, std::vector<size_t> > collapse_set;

  // build map of all pqueries
  for (size_t i = 0; i < sample.size(); ++i) {
    const auto& pq = sample[i];
    collapse_set[pq.sequence_id()].emplace_back(i);
  }

  // find all cases of puplicate entries and merge hem
  for (auto& pair : collapse_set) {
    auto pqlist = pair.second;
    // duplicate!
    if (pqlist.size() > 1) {
      // move entries from duplicate to original
      auto& dest = sample[pqlist[0]].data();
      for (size_t i = 1; i < pqlist.size(); ++i) {
        auto& src = sample[pqlist[i]].data();
        dest.reserve(dest.size() + src.size());
        std::move(std::begin(src), std::end(src), std::back_inserter(dest));
        // mark invalid in original sample
        sample[pqlist[i]].sequence_id(invalid);
      }
    }
  }

  // clear the original sample of invalid pqueries
  sample.erase(
    std::remove_if( std::begin(sample),
                    std::end(sample),
                    [invalid = invalid](auto& e){
                      return e.sequence_id() == invalid;}),
    std::end(sample) );

}

/**
 * dummy collapse function for default case of objects where
 * collapsing doesnt make sense
 */
template<class T>
void collapse(T&)
{ }

/**
 * special split function that Splits samples in buckets according to the global sequence ID
 * of their PQueries. The goal is to have them split such that each aggregate node gets their
 * correct set of sequence results (even if that part is empty, which constitutes a null-message)
 * 
 * @param
 * @param
 * @param
 */
template<class T>
void split( const Sample<T>& src, 
            std::vector<Sample<T>>& parts, 
            const unsigned int num_parts)
{
  parts.clear();
  // ensure that there are actually as many parts as specified. We want empty parts to enable null messages
  parts.resize(num_parts);

  for (auto& pq : src) {
    const auto bucket = pq.sequence_id() % num_parts;
    parts[bucket].push_back(pq);
  }

}

/**
 * Splits a <src> into <num_parts> number of equally sized <parts>.
 * <parts> is expected to be empty.
 *
 * REQUIRES T TO BE:
 *   Iterable
 *   has size()
 *   has insert()
 * 
 * @param src    Container to split
 * @param parts     resulting parts vector
 * @param num_parts number of parts
 */
template <class T>
void split( const T& src, 
            std::vector<T>& parts, 
            const unsigned int num_parts)
{
  parts.clear();
  unsigned int chunk_size = ceil(src.size() / static_cast<double>(num_parts));
  auto move_begin = src.begin();

  while (move_begin < src.end()) {
    auto move_end = std::min(move_begin + chunk_size, src.end());
    parts.emplace_back();
    parts.back().insert(move_begin, move_end);
    move_begin = move_end;
  }
}

void split( const Work& source, 
            std::vector<Work>& parts, 
            const unsigned int num_parts);

/**
  Merges a Sample <src> into a Sample <dest>. Leaves <src> intact.
*/
template<class T>
void merge(Sample<T>& dest, const Sample<T>& src)
{
  // merge in every source pquery...
  for (const auto& pquery : src) {
    // ... by checking if its sequence already exists in destination
    auto input_iter = find(dest.begin(), dest.end(), pquery);
    // if not, create a record
    if (input_iter == dest.end()) {
      dest.emplace_back(pquery.sequence_id(), pquery.header());
      input_iter = --(dest.end());
    }
    // then concat their vectors
    input_iter->insert(input_iter->end(), pquery.begin(), pquery.end());
  }
}

/**
  Merges a Sample <src> into a Sample <dest>. src here is an rvalue,
  and thus the elements are moved instead of copied
*/
template<class T>
void merge(Sample<T>& dest, Sample<T>&& src)
{
  for (auto& pquery : src) {
    // create new record
    dest.emplace_back(std::move(pquery));
  }
}

template<class T>
void merge(T& dest, std::vector<T>&& parts)
{
  for (auto& p : parts) {
    merge(dest, std::move(p));
  }
}

template <class T>
void merge(T& dest, std::vector<T>& parts)
{
  for (const auto& p : parts) {
    merge(dest, p);
  }
}

template <class T>
void merge(std::vector<T>& dest, const std::vector<T>& parts)
{
  dest.insert( std::end(dest), std::begin(parts), std::end(parts) );
}
