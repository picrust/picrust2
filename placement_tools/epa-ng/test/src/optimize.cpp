#include "Epatest.hpp"

#include <vector>

#include "core/pll/optimize.hpp"
#include "core/pll/epa_pll_util.hpp"
#include "io/file_io.hpp"
#include "util/Options.hpp"
#include "tree/Tree.hpp"
#include "tree/Tree_Numbers.hpp"

// TEST(optimize, repeats)
// {
//   // Options options;
//   Tree_Numbers nums;
//   raxml::Model model("GTR+G");
//   // auto tree_file = env->data_dir + "lucas/20k.newick";
//   // auto reference_file = env->data_dir + "lucas/1k_reference.fasta";
//   auto tree_file = env->tree_file;
//   auto reference_file = env->reference_file;

//   auto ref_msa = build_MSA_from_file(env->reference_file, MSA_Info(env->reference_file), false);
//   auto utree = build_tree_from_file(tree_file, nums);
//   auto part = build_partition_from_file(model, nums, ref_msa.num_sites(), true);
//   link_tree_msa(utree, 
//                 part, 
//                 model, 
//                 ref_msa, 
//                 nums.tip_nodes);

//   set_unique_clv_indices(get_root(utree), nums.tip_nodes);

//   optimize( model, 
//             utree, 
//             part, 
//             nums, 
//             true, 
//             true);

//   // Test: check all possible edge logl and throw if the variance is above some threshold
//   std::vector<pll_unode_t*> node_list(nums.branches);
//   utree_query_branches(utree, &node_list[0]);

//   std::vector<unsigned int> param_indices(part->rate_cats, 0);
//   std::vector<pll_unode_t*> travbuffer(nums.nodes);
//   std::vector<double> branch_lengths(nums.branches);
//   std::vector<unsigned int> matrix_indices(nums.branches);
//   std::vector<pll_operation_t> operations(nums.nodes);

//   std::vector<double> logl;

//   precompute_clvs(utree, part, nums);


//   for (auto& n : node_list) {

//   //   unsigned int traversal_size = 0;
//   //   unsigned int num_matrices = 0; 
//   //   unsigned int num_ops = 0;
//   //   pll_utree_traverse( n->back, 
//   //                       PLL_TREE_TRAVERSE_POSTORDER,
//   //                       cb_full_traversal,
//   //                       &travbuffer[0], 
//   //                       &traversal_size);

//   //   pll_utree_create_operations(&travbuffer[0],
//   //                               traversal_size,
//   //                               &branch_lengths[0],
//   //                               &matrix_indices[0],
//   //                               &operations[0],
//   //                               &num_matrices,
//   //                               &num_ops);

//   //   pll_update_prob_matrices(part,
//   //                            &param_indices[0],             // use model 0
//   //                            &matrix_indices[0],// matrices to update
//   //                            &branch_lengths[0],
//   //                            num_matrices); // how many should be updated

//   //   /* use the operations array to compute all num_ops inner CLVs. Operations
//   //      will be carried out sequentially starting from operation 0 towrds num_ops-1 */
//   //   pll_update_partials(part, &operations[0], num_ops);

//     logl.push_back(
//     pll_compute_edge_loglikelihood( part,
//                                     n->clv_index,
//                                     n->scaler_index,
//                                     n->back->clv_index,
//                                     n->back->scaler_index,
//                                     n->pmatrix_index,
//                                     &param_indices[0],
//                                     nullptr)
//     );
//   }

//   double previous = logl[0];
//   for (auto l : logl) {
//     EXPECT_DOUBLE_EQ(l, previous);
//     previous = l;
//     // printf("%f\n", l);
//   }

// }