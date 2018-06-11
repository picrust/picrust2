#include "Epatest.hpp"

#include "core/pll/pllhead.hpp"
#include "core/raxml/Model.hpp"
#include "io/file_io.hpp"
#include "tree/Tree_Numbers.hpp"
#include "seq/MSA.hpp"

#include <string>
#include <tuple>

using namespace std;

TEST(file_io, build_MSA_from_file)
{
  auto msa = build_MSA_from_file(env->reference_file, MSA_Info(env->reference_file), true);
  int i = 0;
  for (auto s : msa) {
    i++;
  }

  EXPECT_EQ(i, 8);
  EXPECT_EQ(msa.num_sites(), 705);

}

TEST(file_io, build_partition_from_file)
{
  auto msa = build_MSA_from_file(env->reference_file, MSA_Info(env->reference_file), true);
  Tree_Numbers nums = Tree_Numbers();
  pll_partition_t * part;
  pll_utree_t * tree;

  tree = build_tree_from_file(env->tree_file, nums);
  part = build_partition_from_file( env->model, nums, msa.num_sites());

  EXPECT_EQ(nums.tip_nodes, 8);
  EXPECT_EQ(nums.nodes, 14);
  EXPECT_EQ(nums.inner_nodes, 6);
  EXPECT_EQ(nums.branches, 13);

  pll_partition_destroy(part);
  pll_utree_destroy(tree, nullptr);
}

TEST(file_io, file_check)
{
  EXPECT_ANY_THROW(file_check("asjbjibvi.hhs"));

  EXPECT_NO_THROW(file_check(env->combined_file));
  EXPECT_NO_THROW(file_check(env->query_file));
  EXPECT_NO_THROW(file_check(env->reference_file));
  EXPECT_NO_THROW(file_check(env->tree_file));
}
