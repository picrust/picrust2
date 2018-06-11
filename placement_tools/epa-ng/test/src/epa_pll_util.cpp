#include "Epatest.hpp"

#include "core/pll/pllhead.hpp"
#include "core/pll/epa_pll_util.hpp"
#include "core/pll/pll_util.hpp"
#include "core/pll/optimize.hpp"
#include "io/file_io.hpp"
#include "tree/Tree_Numbers.hpp"
#include "tree/Tree.hpp"
#include "core/raxml/Model.hpp"
#include "util/Options.hpp"
#include "seq/MSA.hpp"

#include <string>

using namespace std;

TEST(epa_pll_util, link_tree_msa)
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
  vector<pll_unode_t *> tip_nodes(nums.tip_nodes);
  tip_nodes.assign(tree->nodes, tree->nodes + nums.tip_nodes);

  for (auto n : tip_nodes) {
    ASSERT_NE(n, nullptr);
    if (part->attributes & PLL_ATTRIB_PATTERN_TIP)
      EXPECT_NE(part->tipchars[n->clv_index][0], 0);
    else
    {
      EXPECT_NE(part->clv[n->clv_index][0], 0);
    }
  }

  // teardown
  pll_partition_destroy(part);
  pll_utree_destroy(tree, nullptr);

}

static void precompute_clvs_test(Options o)
{
  // auto tree_file = env->data_dir + "lucas/20k.newick";
  // auto reference_file = env->data_dir + "lucas/1k_reference.fasta";
  auto tree_file = env->tree_file;
  auto reference_file = env->reference_file;
  // buildup
  auto msa = build_MSA_from_file(env->reference_file, MSA_Info(env->reference_file), o.premasking);
  Tree_Numbers nums;
  raxml::Model model;


  auto tree = build_tree_from_file( tree_file, nums);
  auto part = build_partition_from_file( model, nums, msa.num_sites(), o.repeats);


  auto root = get_root(tree);
  set_unique_clv_indices(root, nums.tip_nodes);
  // auto valid_map = vector<Range>(nums.tip_nodes);
  link_tree_msa(tree, part, model, msa, nums.tip_nodes);

  optimize( model, 
            tree, 
            part, 
            nums, 
            o.opt_branches, 
            o.opt_model);

  precompute_clvs(tree, part, nums);

  // tests
  vector<pll_unode_t *> node_list(nums.branches, nullptr);
  utree_query_branches(tree, &node_list[0]);

  for (auto& p : node_list) {
    ASSERT_NE(p, nullptr);
  }

  // all edge logl should be the same
  auto first = true;
  double log_old = 0.0;
  double log_new;
  unsigned int param_indices[4] = {0};
  size_t id = 0;
  size_t success = 1;
  size_t failure = 0;
  for (auto x : node_list)
  {
    log_new = pll_compute_edge_loglikelihood(part,
                                         x->clv_index,
                                         x->scaler_index,
                                         x->back->clv_index,
                                         x->back->scaler_index,
                                         x->pmatrix_index,
                                         param_indices, nullptr);
    if (!first) {
      EXPECT_DOUBLE_EQ(log_old, log_new);
      if (fabs(log_old - log_new) > 1e-10) {
        failure++;

        printf("Failure!\n");
        printf("new logl: %f\n", log_new);
        printf("IDX is: %lu\n", id);
        printf("x->clv_index: \t\t%u\t%p\n", x->clv_index, part->clv[x->clv_index]);
        printf("x->scaler_index: \t%d\t%p\n", x->scaler_index,part->scale_buffer[x->scaler_index]);
        printf("x->back->clv_index: \t%u\t%p\n", x->back->clv_index, part->clv[x->back->clv_index]);
        printf("x->back->scaler_index: \t%d\t%p\n", x->back->scaler_index, part->scale_buffer[x->back->scaler_index]);
        printf("Is tip? %d\n", (x->next == nullptr));
        printf("Back tip? %d\n", (x->back->next == nullptr));
        printf("\n");
      } else {
        success++;
      }
    } else {
      log_old = log_new;
    }

    // log_old = log_new;
    first = false;
    id++;
  }

  // printf("success vs failue: %lu vs %lu\n", success, failure);

  EXPECT_NE(log_new, 0.0);

  // teardown
  utree_free_node_data(get_root(tree));
  pll_partition_destroy(part);
  pll_utree_destroy(tree, nullptr);
}

TEST(epa_pll_util, precompute_clvs)
{
  Options o;
  o.opt_model = o.opt_branches = false;
  o.repeats = false;
  precompute_clvs_test(o);  
}

TEST(epa_pll_util, split_combined_msa)
{
  // buildup
  auto combined_msa = build_MSA_from_file(env->combined_file, MSA_Info(env->combined_file), true);
  raxml::Model model;
  Options options;

  Tree tree(env->tree_file, combined_msa, model, options);

  // tests
  MSA query_msa;
  split_combined_msa(combined_msa, query_msa, tree);

  EXPECT_EQ(combined_msa.num_sites(), query_msa.num_sites());
  EXPECT_EQ(combined_msa.size(), 8);
  EXPECT_EQ(query_msa.size(), 2);

  for (auto x : combined_msa) {
    for (auto y : query_msa) {
      EXPECT_STRNE(x.header().c_str(), y.header().c_str());
    }
  }

  // teardown
}
