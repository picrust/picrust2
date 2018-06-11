#include "Epatest.hpp"

#include "core/pll/pllhead.hpp"
#include "core/pll/pll_util.hpp"
#include "core/pll/epa_pll_util.hpp"
#include "io/file_io.hpp"
#include "tree/Tree_Numbers.hpp"
#include "core/raxml/Model.hpp"
#include "seq/MSA.hpp"

#include <string>
#include <tuple>

using namespace std;

TEST(pll_util, utree_query_branches)
{
  // buildup
  auto msa = build_MSA_from_file(env->reference_file, MSA_Info(env->reference_file), true);
  Tree_Numbers nums = Tree_Numbers();
  pll_partition_t * part;
  pll_utree_t * tree;

  tree = build_tree_from_file(env->tree_file, nums);
  part = build_partition_from_file( env->model, nums, msa.num_sites());

  // tests
  vector<pll_unode_t *> node_list(nums.branches);
  auto num_traversed = utree_query_branches(tree, &node_list[0]);

  EXPECT_EQ(nums.branches, num_traversed);

  for (auto x : node_list) {
    EXPECT_NE(x, nullptr);
  }

  // check for duplicates
  for (auto x : node_list) {
    int count = 0;
    for (auto y : node_list) {
      if (x == y || x == y->back)
        count++;
    }
    EXPECT_EQ(count, 1);
  }

  // teardown
  pll_partition_destroy(part);
  pll_utree_destroy(tree, nullptr);
}

TEST(pll_util, set_unique_clv_indices)
{
  // buildup
  auto msa = build_MSA_from_file(env->reference_file, MSA_Info(env->reference_file), true);
  Tree_Numbers nums = Tree_Numbers();

  auto tree = build_tree_from_file(env->tree_file, nums);
  auto part = build_partition_from_file( env->model, nums, msa.num_sites());

  set_unique_clv_indices(get_root(tree), nums.tip_nodes);

  // tests
  vector<pll_unode_t *> tip_nodes(nums.tip_nodes);
  tip_nodes.assign(tree->nodes, tree->nodes + nums.tip_nodes);

  vector<pll_unode_t *> inner_nodes(nums.inner_nodes);
  inner_nodes.assign(tree->nodes + nums.tip_nodes, tree->nodes + nums.nodes);

  // check for duplicate clv indices
  vector<unsigned int> index_list;
  for (auto x : tip_nodes)
    index_list.push_back(x->clv_index);

  for (auto x : inner_nodes) {
    index_list.push_back(x->clv_index);
    index_list.push_back(x->next->clv_index);
    index_list.push_back(x->next->next->clv_index);
  }

  for (auto x : index_list) {
    int count = 0;
    for (auto y : index_list) {
      if (x == y)
        count++;
    }
    EXPECT_EQ(count, 1);
  }

  // teardown
  pll_partition_destroy(part);
  pll_utree_destroy(tree, nullptr);
}

static int cb_set_branchlengths_one(pll_unode_t * node)
{
  node->length = 1;
  if(node->next)
  {
    node->next->length = 1;
    node->next->next->length = 1;
  }
  return 1;
};

TEST(pll_util, get_numbered_newick_string)
{
  // buildup
  auto msa = build_MSA_from_file(env->reference_file, MSA_Info(env->reference_file), true);
  Tree_Numbers nums = Tree_Numbers();
  pll_partition_t * part;
  pll_utree_t * tree;
  raxml::Model model;

  tree = build_tree_from_file(env->tree_file, nums);
  part = build_partition_from_file( env->model, nums, msa.num_sites());
  // auto valid_map = vector<Range>(nums.tip_nodes);
  link_tree_msa(tree, part, model, msa, nums.tip_nodes);

  // tests
  // valid output as returned by RAxML, with reset branch lengths, as we only want to test format
  string valid(
  "(Seal:1{0},(Whale:1{1},(Mouse:1{2},(Human:1{3},(Chicken:1{4},(Frog:1{5},Loach:1{6}):1{7}):1{8}):1{9}):1{10}):1{11},Cow:1{12});");

  vector<pll_unode_t*> travbuffer(nums.nodes);
  unsigned int traversal_size;
  pll_utree_traverse( get_root(tree), 
                      PLL_TREE_TRAVERSE_POSTORDER,
                      cb_set_branchlengths_one, 
                      &travbuffer[0], 
                      &traversal_size);

  auto ret_string =  get_numbered_newick_string(tree);

  EXPECT_STREQ(valid.c_str(), ret_string.c_str());

  // teardown

  pll_partition_destroy(part);
  pll_utree_destroy(tree, nullptr);
}

TEST(pll_util, sum_branch_lengths)
{
  auto msa = build_MSA_from_file(env->reference_file, MSA_Info(env->reference_file), true);
  Tree_Numbers nums = Tree_Numbers();
  pll_partition_t * part;
  pll_utree_t * tree;

  tree = build_tree_from_file(env->tree_file, nums);
  part = build_partition_from_file( env->model, nums, msa.num_sites());
  set_branch_lengths(tree, 1.0);
  auto total_length = sum_branch_lengths(tree);

  EXPECT_DOUBLE_EQ(nums.branches, total_length);

  pll_partition_destroy(part);
  pll_utree_destroy(tree, nullptr);
}

TEST(pll_util, shift_partition_focus_shifty)
{
  // buildup
  pll_partition_t * part;
  part =  pll_partition_create(
    3, // tips
    1, // extra clv's
    4, 10,
    1,
    3, // number of prob. matrices (one per possible unique branch length)
    4,
    1, // number of scale buffers (one per possible inner node)
    // pll_map_nt,
    PLL_ATTRIB_PATTERN_TIP | PLL_ATTRIB_ARCH_CPU);


  // tests

  int clv_size = part->states * part->rate_cats;
  int begin = 2;
  int span = 4;
  
  for (size_t i = 0; i < part->sites * clv_size; i++) {
    if ((int)i < begin * clv_size or (int)i >= (begin + span) * clv_size)
      part->clv[3][i] = 2.0;
    else
      part->clv[3][i] = 1.0;
  }

  auto seq = new char[part->sites];

  for (size_t i = 0; i < part->sites; i++) 
  {
    if ((int)i < begin or (int)i >= (begin + span))
    {
      seq[i] = 'A';
    }
    else
    {
      seq[i] = 'G';
    }
  }

  // printf("Sequence: %s\n", seq);

  pll_set_tip_states(part, 0, pll_map_nt, seq);
  ASSERT_NE(part->tipchars, nullptr);
  ASSERT_NE(part->tipchars[0], nullptr);

  shift_partition_focus(part, begin, span);

  // CLV test
  for (size_t i = 0; i < part->sites * clv_size; i++) {
    EXPECT_DOUBLE_EQ(part->clv[3][i], 1.0);
  }
  EXPECT_DOUBLE_EQ(part->clv[3][-1], 2.0);
  EXPECT_DOUBLE_EQ(part->clv[3][part->sites * clv_size], 2.0);

  // tipchars test
  for (size_t i = 0; i < part->sites; i++) {
    EXPECT_EQ(part->tipchars[0][i], '\x4');
  }
  EXPECT_EQ(part->tipchars[0][-1], '\x1');
  EXPECT_EQ(part->tipchars[0][part->sites], '\x1');


  // teardown
  shift_partition_focus(part, -begin, 10);
  delete[] seq;
  pll_partition_destroy(part);
}

TEST(pll_util, shift_partition_focus_logtest)
{
  // buildup
  pll_partition_t * part;
  part =  pll_partition_create(
    3, // tips
    1, // extra clv's
    4, 10,
    1,
    3, // number of prob. matrices (one per possible unique branch length)
    4,
    1, // number of scale buffers (one per possible inner node)
    // pll_map_nt,
    PLL_ATTRIB_PATTERN_TIP | PLL_ATTRIB_ARCH_CPU);

  double frequencies[4] = { 0.17, 0.19, 0.25, 0.39 };

  /* substitution rates for the GTR model */
  double subst_params[6] = {1,1,1,1,1,1};

  /* discretized category rates from a gamma distribution with alpha shape 1 */
  double rate_cats[4];
  pll_compute_gamma_cats(1.0, 4, rate_cats, PLL_GAMMA_RATES_MEAN);

  /* set frequencies */
  pll_set_frequencies(part, 0, frequencies);

  /* set substitution parameters */
  pll_set_subst_params(part, 0, subst_params);

  /* set rate categories */
  pll_set_category_rates(part, rate_cats);

  /* set the 5 tip CLVs, and use the pll_map_nt map for converting
     the sequences to CLVs */
  pll_set_tip_states(part, 0, pll_map_nt, "--AAAA----");
  pll_set_tip_states(part, 1, pll_map_nt, "--TTTT----");
  pll_set_tip_states(part, 2, pll_map_nt, "--GGGG----");

  double branch_lengths[3] = {0.123, 0.123, 0.123};
  unsigned int matrix_indices[3] = {0, 1, 2};
  unsigned int param_indices[4] = {0};
  pll_update_prob_matrices(part, param_indices, matrix_indices, branch_lengths, 3);

  pll_operation_t op;
  op.parent_clv_index = 3;
  op.child1_clv_index = 1;
  op.child1_scaler_index = PLL_SCALE_BUFFER_NONE;
  op.child2_clv_index = 2;
  op.child2_scaler_index = PLL_SCALE_BUFFER_NONE;
  op.parent_scaler_index = 0;
  op.child1_matrix_index = 1;
  op.child2_matrix_index = 2;

  pll_update_partials(part, &op, 1);

  // tests
  double full_logl = pll_compute_edge_loglikelihood(part,
                                        0,
                                        PLL_SCALE_BUFFER_NONE,
                                        3,
                                        0,
                                        0,
                                        param_indices, nullptr);

  shift_partition_focus(part, 2, 4);

  pll_update_partials(part, &op, 1);

  double ranged_logl = pll_compute_edge_loglikelihood(part,
                                        0,
                                        PLL_SCALE_BUFFER_NONE,
                                        3,
                                        0,
                                        0,
                                        param_indices, nullptr);

  EXPECT_DOUBLE_EQ(full_logl, ranged_logl);

  shift_partition_focus(part, 0, 3);

  pll_update_partials(part, &op, 1);

  double false_logl = pll_compute_edge_loglikelihood(part,
                                        0,
                                        PLL_SCALE_BUFFER_NONE,
                                        3,
                                        0,
                                        0,
                                        param_indices, nullptr);

  EXPECT_TRUE(abs(full_logl - false_logl) > 1.0 );
  // printf("Full:\t %.40f\n", full_logl);
  // printf("Ranged:\t %.40f\n", ranged_logl);
  // printf("False:\t %.40f\n", false_logl);


// extra tests
  shift_partition_focus(part, 3, 1);

  pll_update_partials(part, &op, 1);

  false_logl += pll_compute_edge_loglikelihood(part,
                                        0,
                                        PLL_SCALE_BUFFER_NONE,
                                        3,
                                        0,
                                        0,
                                        param_indices, nullptr);
  // printf("False+missing:\t %.40f\n", false_logl);

  shift_partition_focus(part, 1, 4);

  pll_update_partials(part, &op, 1);

  false_logl += pll_compute_edge_loglikelihood(part,
                                        0,
                                        PLL_SCALE_BUFFER_NONE,
                                        3,
                                        0,
                                        0,
                                        param_indices, nullptr);
  // printf("f+m+right:\t %.40f\n", false_logl);

  shift_partition_focus(part, -6, 2);

  pll_update_partials(part, &op, 1);

  false_logl += pll_compute_edge_loglikelihood(part,
                                        0,
                                        PLL_SCALE_BUFFER_NONE,
                                        3,
                                        0,
                                        0,
                                        param_indices, nullptr);
  // printf("fmr+left:\t %.40f\n", false_logl);

  shift_partition_focus(part, 0, 10);

  // teardown
  pll_partition_destroy(part);
}
