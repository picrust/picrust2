#pragma once

#include <chrono>
#include <vector>
#include <algorithm>
#include <cereal/types/vector.hpp>
#include <cereal/types/chrono.hpp>

template <class duration = std::chrono::seconds>
class Timer {
public:
  // Typedefs
  using clock           = std::chrono::high_resolution_clock;
  using iterator        = typename std::vector<duration>::iterator;
  using const_iterator  = typename std::vector<duration>::const_iterator;
  
  // Constructors/Destructors
  Timer(duration&& d)
  {
    ts_.push_back(d);
  }

  Timer()   = default;
  ~Timer()  = default;

  // Iterator Compatibility
  iterator begin() { return ts_.begin(); }
  iterator end() { return ts_.end(); }
  const_iterator begin() const { return ts_.cbegin(); }
  const_iterator end() const { return ts_.cend(); }
  const_iterator cbegin() { return ts_.cbegin(); }
  const_iterator cend() { return ts_.cend(); }

  // Methods
  void insert(iterator position, const_iterator first, const_iterator last)
  { 
    ts_.insert(position, first, last); 
  };

  void start() 
  {
    start_ = clock::now();
  }

  void pause()
  {
    pause_start_ = clock::now();
  }

  void resume()
  {
    auto end = clock::now();
    auto pause_time = std::chrono::duration_cast<duration>(end - pause_start_);
    pauses_.push_back(pause_time);
  }

  void stop()
  {
    auto end = clock::now();

    duration pause_total(this->sum_pauses()); 

    auto runtime = std::chrono::duration_cast<duration>(end - start_) - pause_total;

    ts_.push_back(runtime);
    pauses_.clear();
  }

  double sum() const
  {
    return this->sum_duration().count();
  }

  duration sum_duration() const
  {
    duration sum(0);
    for (auto p : ts_) {
      sum += p;
    }
    return sum;
  }

  duration sum_pauses() const
  {
    duration pause_total(0);
    for (auto p : pauses_) {
      pause_total += p;
    }
    return pause_total;
  }

  duration avg_duration() const
  {
    return this->sum_duration()/ts_.size(); 
  }

  double average() const
  {
    return this->sum()/ts_.size();
  }

  void clear() {ts_.clear();}

  // Serialization
  template <class Archive>
  void save(Archive & ar) const { ar( ts_ ); }

  template <class Archive>
  void load(Archive & ar) { ar( ts_ ); }
private:
  clock::time_point start_;
  std::vector<duration> ts_;
  clock::time_point pause_start_;
  std::vector<duration> pauses_;
};
