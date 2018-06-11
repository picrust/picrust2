#include "core/pll/epa_pll_util.hpp"

#include <unordered_map>
#include <algorithm>

#include "core/pll/pll_util.hpp"
#include "set_manipulators.hpp"

void link_tree_msa( pll_utree_t * tree, 
                    pll_partition_t * partition, 
                    raxml::Model& model, 
                    const MSA& msa, 
                    const unsigned int num_tip_nodes)
{
  // associate the sequences from the MSA file with the correct tips
  /* create a hash table of size num_tip_nodes */
  std::unordered_map<std::string, unsigned int> map; // mapping labels to tip clv indices

  /* populate the hash table with tree tip labels */
  for (size_t i = 0; i < num_tip_nodes; ++i) {
    map[tree->nodes[i]->label] = i;
  }

  /* find sequences in hash table and link them with the corresponding taxa */
  for (auto const &s : msa) {
    auto map_value = map.find(s.header());

    // failure tolerance: the MSA may also contain query sequences
    if (map_value == map.end()) {
      continue;
      // throw runtime_error{std::string("Sequence with header does not appear in the tree: ") + s.header()};
    }

    auto clv_index = map_value->second;
    // associates the sequence with the tip by calculating the tips clv buffers
    pll_set_tip_states(partition, clv_index, model.charmap(), s.sequence().c_str());
  }
}

void precompute_clvs( pll_utree_t const * const tree, 
                      pll_partition_t * partition, 
                      const Tree_Numbers& nums)
{
  /* various buffers for creating a postorder traversal and operations structures */
  std::vector<unsigned int> param_indices(partition->rate_cats, 0);
  std::vector<pll_unode_t*> travbuffer(nums.nodes);
  std::vector<double> branch_lengths(nums.branches);
  std::vector<unsigned int> matrix_indices(nums.branches);
  std::vector<pll_operation_t> operations(nums.nodes);

  const auto root = get_root(tree);

  /* adjust clv indices such that every direction has its own */
  // set_unique_clv_indices(root, nums.tip_nodes);

  std::vector<std::vector<double>> persite(tree->tip_count);

  utree_free_node_data(root);

  for (size_t i = 0; i < tree->tip_count; ++i) {
    const auto node = tree->nodes[i];
    /* perform a partial postorder traversal of the unrooted tree  starting at the current tip
      and returning every node whose clv in the direction of the tip hasn't been calculated yet*/
    unsigned int traversal_size = 0;
    unsigned int num_matrices = 0;
    unsigned int num_ops = 0;
    if (pll_utree_traverse( node->back,
                            PLL_TREE_TRAVERSE_POSTORDER,
                            cb_partial_traversal,
                            &travbuffer[0], 
                            &traversal_size)
                != PLL_SUCCESS) {
      throw std::runtime_error{"Function pll_unode_traverse() requires inner nodes as parameters"};
    }

    /* given the computed traversal descriptor, generate the operations
       structure, and the corresponding probability matrix indices that
       may need recomputing */
    pll_utree_create_operations(&travbuffer[0],
                                traversal_size,
                                &branch_lengths[0],
                                &matrix_indices[0],
                                &operations[0],
                                &num_matrices,
                                &num_ops);

    pll_update_prob_matrices(partition,
                             &param_indices[0],             // use model 0
                             &matrix_indices[0],// matrices to update
                             &branch_lengths[0],
                             num_matrices); // how many should be updated

    /* use the operations array to compute all num_ops inner CLVs. Operations
       will be carried out sequentially starting from operation 0 towrds num_ops-1 */
    pll_update_partials(partition, &operations[0], num_ops);
  }
  utree_free_node_data(root);
}

void split_combined_msa(MSA& source, 
                        MSA& target, 
                        Tree& tree)
{
  std::vector<pll_unode_t*> tip_nodes(tree.nums().tip_nodes);
  tip_nodes.assign( tree.tree()->nodes, 
                    tree.tree()->nodes + tree.tree()->tip_count);

  auto falsegroup_begin = partition(source.begin(), source.end(),
    [tip_nodes = std::move(tip_nodes)](const Sequence& em) {
      return find(tip_nodes.begin(), tip_nodes.end(), em) != tip_nodes.end();
    });
  target.num_sites(source.num_sites());
  target.move_sequences(falsegroup_begin, source.end());
  source.erase(falsegroup_begin, source.end());
}

bool operator==(const pll_unode_t * node, const Sequence& s)
{
  return s.header().compare(node->label) == 0;
}

bool operator==(const Sequence& s, const pll_unode_t * node)
{
  return operator==(node, s);
}

raxml::Model get_model(pll_partition_t* partition)
{
  using namespace raxml;

  DataType seqtype = DataType::autodetect;

  if (partition->states == 4) {
    seqtype = DataType::dna;
  } else if (partition->states == 20) {
    seqtype = DataType::protein;
  } else {
    throw std::runtime_error{"Couldn't determine sequence type from partition"};
  }

  Model model(seqtype);

  assign(model, partition);

  return model;
}
