#include "tree/Tiny_Tree.hpp"

#ifdef __OMP
#include <omp.h>
#endif

#include <vector>
#include <numeric>

#include "tree/tiny_util.hpp"
#include "core/pll/pll_util.hpp"
#include "core/pll/optimize.hpp"
#include "core/raxml/Model.hpp"
#include "tree/Tree_Numbers.hpp"
#include "set_manipulators.hpp"
#include "util/logging.hpp"

static void precompute_sites_static(char nt,
                                    std::vector<double>& result,
                                    pll_partition_t * const partition,
                                    pll_utree_t const * const tree)
{
  const size_t sites  = partition->sites;
  const auto new_tip  = tree->nodes[2];
  const auto inner    = new_tip->back;
  result.clear();
  result.resize(sites);
  std::string seq(sites, nt);

  std::vector<unsigned int> param_indices(partition->rate_cats, 0);

  auto map = get_char_map(partition);

  auto err_check = pll_set_tip_states(partition, 
                                      new_tip->clv_index, 
                                      map,
                                      seq.c_str());

  if (err_check == PLL_FAILURE) {
    throw std::runtime_error{
      std::string("Set tip states during sites precompution failed! pll_errmsg: ")
      + pll_errmsg
    };
  }

  pll_compute_edge_loglikelihood( partition,
                                  new_tip->clv_index,
                                  PLL_SCALE_BUFFER_NONE, 
                                  inner->clv_index,
                                  inner->scaler_index,
                                  inner->pmatrix_index,
                                  &param_indices[0], 
                                  &result[0]);
}


Tiny_Tree::Tiny_Tree( pll_unode_t * edge_node, 
                      const unsigned int branch_id, 
                      Tree& reference_tree, 
                      const bool opt_branches, 
                      const Options& options, 
                      std::shared_ptr<Lookup_Store>& lookup_store)
  : partition_(nullptr, tiny_partition_destroy)
  , tree_(nullptr, utree_destroy)
  , opt_branches_(opt_branches)
  , premasking_(options.premasking)
  , sliding_blo_(options.sliding_blo)
  , branch_id_(branch_id)
  , lookup_(lookup_store)
{
  assert(edge_node);
  original_branch_length_ = edge_node->length;

  auto old_proximal = edge_node->back;
  auto old_distal = edge_node;

  // detect the tip-tip case. In the tip-tip case, the reference tip should
  // always be the DISTAL
  bool tip_tip_case = false;
  if (!old_distal->next) {
    tip_tip_case = true;
  } else if (!old_proximal->next) {
    tip_tip_case = true;
    // do the switcheroo
    old_distal = old_proximal;
    old_proximal = old_distal->back;
  }

  tree_ = std::unique_ptr<pll_utree_t, utree_deleter>(
      	                    make_tiny_tree_structure( old_proximal,
                                                      old_distal,
                                                      tip_tip_case),
                            utree_destroy);

  partition_ = std::unique_ptr<pll_partition_t, partition_deleter>(
                                make_tiny_partition(reference_tree,
                                                    tree_.get(),
                                                    old_proximal,
                                                    old_distal,
                                                    tip_tip_case),
                                tiny_partition_destroy);

  // operation for computing the clv toward the new tip (for initialization and logl in non-blo case)
  auto proximal = tree_->nodes[0];
  auto distal   = tree_->nodes[1];
  auto inner    = tree_->nodes[3];

  pll_operation_t op;
  op.parent_clv_index = inner->clv_index;
  op.child1_clv_index = distal->clv_index;
  op.child1_scaler_index = distal->scaler_index;
  op.child2_clv_index = proximal->clv_index;
  op.child2_scaler_index = proximal->scaler_index;
  op.parent_scaler_index = inner->scaler_index;
  op.child1_matrix_index = distal->pmatrix_index;
  op.child2_matrix_index = proximal->pmatrix_index;

  // wether heuristic is used or not, this is the initial branch length configuration
  double branch_lengths[3] = {proximal->length, distal->length, inner->length};
  unsigned int matrix_indices[3] = {proximal->pmatrix_index, distal->pmatrix_index, inner->pmatrix_index};

  // use branch lengths to compute the probability matrices
  std::vector<unsigned int> param_indices(reference_tree.partition()->rate_cats, 0);
  pll_update_prob_matrices( partition_.get(), 
                            &param_indices[0], 
                            matrix_indices, 
                            branch_lengths, 
                            3);

  // use update_partials to compute the clv pointing toward the new tip
  pll_update_partials(partition_.get(), &op, 1);

  if (not opt_branches) {
    const std::lock_guard<std::mutex> lock(lookup_store->get_mutex(branch_id));

    if (not lookup_store->has_branch(branch_id)) {
      const auto size = lookup_store->char_map_size();

      // precompute all possible site likelihoods
      std::vector<std::vector<double>> precomputed_sites(size);
      for (size_t i = 0; i < size; ++i) {
        precompute_sites_static(lookup_store->char_map(i),
                                precomputed_sites[i],
                                partition_.get(),
                                tree_.get());
      }
      lookup_store->init_branch(branch_id, precomputed_sites);
    }
  }
  
}

Placement Tiny_Tree::place(const Sequence &s) 
{
  assert(partition_);
  assert(tree_);

  const auto inner    = tree_->nodes[3];
  const auto distal   = tree_->nodes[1];
  const auto proximal = tree_->nodes[0];
  const auto new_tip  = inner->back;

  auto distal_length = distal->length;
  auto pendant_length = inner->length;
  double logl = 0.0;
  std::vector<unsigned int> param_indices(partition_->rate_cats, 0);

  if ( s.sequence().size() != partition_->sites ) {
    throw std::runtime_error{"Query sequence length not same as reference alignment!"};
  }

  Range range(0, partition_->sites);

  if (premasking_) {
    range = get_valid_range(s.sequence());
  }

  if (opt_branches_) {

    auto virtual_root = inner;

    // init the new tip with s.sequence(), branch length
    auto err_check = pll_set_tip_states(partition_.get(), 
                                        new_tip->clv_index, 
                                        get_char_map(partition_.get()),
                                        s.sequence().c_str());

    if (err_check == PLL_FAILURE) {
      throw std::runtime_error{"Set tip states during placement failed!"};
    }

    if (premasking_){
      logl = call_focused(optimize_branch_triplet, range, partition_.get(), virtual_root, sliding_blo_);
    } else {
      logl = optimize_branch_triplet(partition_.get(), virtual_root, sliding_blo_);
    }

    assert(inner->length >= 0);
    assert(inner->next->length >= 0);
    assert(inner->next->next->length >= 0);

    // rescale the distal length, as it has likely changed during optimization
    // done as in raxml
    const double new_total_branch_length = distal->length + proximal->length;
    distal_length = (original_branch_length_ / new_total_branch_length) * distal->length;
    pendant_length = inner->length;

    reset_triplet_lengths(inner,
                          partition_.get(), 
                          original_branch_length_);
    
    // re-update the partial
    auto child1 = virtual_root->next->back;
    auto child2 = virtual_root->next->next->back;

    pll_operation_t op;
    op.parent_clv_index = virtual_root->clv_index;
    op.parent_scaler_index = virtual_root->scaler_index;
    op.child1_clv_index = child1->clv_index;
    op.child1_scaler_index = child1->scaler_index;
    op.child1_matrix_index = child1->pmatrix_index;
    op.child2_clv_index = child2->clv_index;
    op.child2_scaler_index = child2->scaler_index;
    op.child2_matrix_index = child2->pmatrix_index;

    pll_update_partials(partition_.get(), &op, 1);

  } else {
    logl = lookup_->sum_precomputed_sitelk(branch_id_, s.sequence(), range);
  }

  assert(distal_length <= original_branch_length_);
  assert(distal_length >= 0.0);

  return Placement(branch_id_, logl, pendant_length, distal_length);
}
