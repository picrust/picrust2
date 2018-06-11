#include "pipeline/schedule.hpp"

#include <stdexcept>
#include <cassert>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <iterator>
#include <vector>

void to_difficulty(std::vector<double>& perstage_avg)
{
  auto min = *std::min_element(perstage_avg.begin(), perstage_avg.end());
  for_each(perstage_avg.begin(), perstage_avg.end(), 
    [min](double& x){x /= min;}
  );
}

std::vector<unsigned int> solve(unsigned int stages, 
                                unsigned int nodes, 
                                const std::vector<double>& difficulty_per_stage)
{
  assert(difficulty_per_stage.size() == stages);
  if (nodes < stages) {
    throw std::runtime_error{"Must have more or equal number of nodes than stages"};
  }

  std::vector<unsigned int> nodes_per_stage(stages);

  // TODO for the record the following is ugly and I hate it
  auto constrained_begin = std::begin(difficulty_per_stage);
  // std::advance(constrained_begin, 1);
  auto constrained_end = std::end(difficulty_per_stage);
  // std::advance(constrained_end, -1);
  
  auto x1 = std::accumulate(constrained_begin, constrained_end, 0.0);
  x1 = static_cast<double>(nodes) / x1;

  for (size_t i = 0; i < stages; ++i) {
    if (i == 0 or i == stages - 1) {
      nodes_per_stage[i] = 1;
    } else {
      nodes_per_stage[i] = ceil(difficulty_per_stage[i] * x1);
    }
  }

  int off_by = 0;

  while ( (off_by = std::accumulate(nodes_per_stage.begin(), nodes_per_stage.end(), 0) - nodes) ) {
    auto max_stage = std::max_element(nodes_per_stage.begin(), nodes_per_stage.end());
    if (off_by < 0) {
      *max_stage += 1;
    } else {
      *max_stage -= 1;
    }
  }

  return nodes_per_stage;
}

void assign(const int local_rank,
            std::vector<unsigned int>& nodes_per_stage, 
            schedule_type& rank_assignm,
            int* local_stage)
{
  rank_assignm.clear();
  int rank = 0;
  for (size_t stage = 0; stage < nodes_per_stage.size(); ++stage) {
    auto nodes = nodes_per_stage[stage];
    rank_assignm.push_back(std::vector<int>(nodes));
    for (size_t j = 0; j < nodes; ++j) {
      rank_assignm.back()[j] = rank;
      if (local_rank == rank) {
        *local_stage = stage;
      }
      rank++;
    }
  }
}

void reassign(const int local_rank,
              std::vector<unsigned int>& nodes_per_stage, 
              schedule_type& rank_assignm,
              int* local_stage)
{
  assert(nodes_per_stage.size() == rank_assignm.size());
  // extract ranks from stages that have too many
  std::vector<int> cut_ranks;
  for (size_t i = 0; i < nodes_per_stage.size(); ++i) {
    auto& cur_stage = rank_assignm[i];
    int to_rm = (int)cur_stage.size() - nodes_per_stage[i] ;
    auto rm_iter = cur_stage.end();
    for (; to_rm > 0; --to_rm) {
      --rm_iter;
      cut_ranks.push_back(*rm_iter);
    }
    cur_stage.erase(rm_iter, cur_stage.end());
  }

  // reassign them where they are needed
  auto copy_iter = cut_ranks.begin();
  for (size_t i = 0; i < nodes_per_stage.size(); ++i) {
    auto& cur_stage = rank_assignm[i];
    int to_add = (int)nodes_per_stage[i] - cur_stage.size();
    if (to_add > 0)
    {
      auto old_copy_iter = copy_iter;
      std::advance(copy_iter, to_add);
      for (; old_copy_iter != copy_iter; ++old_copy_iter) {
        int rank = *old_copy_iter;
        if (rank == local_rank) {
          *local_stage = i;
        }
        cur_stage.push_back(rank);
      }
    }
  }
}
