#include "Epatest.hpp"

#include "pipeline/schedule.hpp"

#include <vector>
#include <numeric>

using namespace std;

TEST(schedule, solve)
{
  unsigned int stages, nodes;

  stages = 4;
  nodes = 32;

  vector<double> diff{1000.0, 1.0, 1000.0, 1.0};

  auto nps = solve(stages, nodes, diff);

  // for(const auto& n : nps)
  //   printf("%d, ", n);
    
  // printf("\nTotal: %d\n", accumulate(nps.begin(), nps.end(), 0));
  EXPECT_EQ(accumulate(nps.begin(), nps.end(), 0), nodes);
}

TEST(schedule, to_difficulty)
{
  std::vector<double> perstage_avg = {20.0, 2.0, 10.0, 3.0};
  std::vector<double> expected = {10.0, 1.0, 5.0, 1.5};

  to_difficulty(perstage_avg);

  for (unsigned int i = 0; i < perstage_avg.size(); ++i)
  {
    EXPECT_DOUBLE_EQ(perstage_avg[i], expected[i]);
  }   
  
}

TEST(schedule, assign)
{ 
  int local_stage;
  int local_rank = 0;

  std::vector<std::vector<unsigned int>> nps = { {15, 1, 15, 1}, {1, 1, 1, 1}, {2, 0, 1, 1} };
  std::vector<std::vector<int>> rank_assignm;

  for (auto& snps : nps)
  {
    assign(local_rank, snps, rank_assignm, &local_stage);
    for (unsigned int i = 0; i < rank_assignm.size(); ++i)
    {
      auto stage_assign = rank_assignm[i];
      EXPECT_EQ(stage_assign.size(), snps[i]);
    }

    rank_assignm[0].clear();
    rank_assignm[1].clear();
    rank_assignm[2].clear();
    rank_assignm[3].clear();
  }
  
}

TEST(schedule, reassign)
{ 
  int local_stage;
  int local_rank = 0;

  std::vector<unsigned int> init_nps = {8,8,8,8};
  std::vector<std::vector<unsigned int>> nps = { {15, 1, 15, 1}, {2, 14, 1, 15}, {30, 0, 1, 1} };
  std::vector<std::vector<int>> rank_assignm;
  assign(local_rank, init_nps, rank_assignm, &local_stage);

  for (auto& snps : nps)
  {
    reassign(local_rank, snps, rank_assignm, &local_stage);
    for (unsigned int i = 0; i < rank_assignm.size(); ++i)
    {
      auto stage_assign = rank_assignm[i];
      EXPECT_EQ(stage_assign.size(), snps[i]);
    }
  }  
}
