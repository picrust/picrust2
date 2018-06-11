#include "tree/Tree.hpp"

#include <stdexcept>
#include <iostream>
#include <cstdio>

#include "core/pll/epa_pll_util.hpp"
#include "io/file_io.hpp"
#include "seq/Sequence.hpp"
#include "core/pll/optimize.hpp"
#include "set_manipulators.hpp"
#include "util/logging.hpp"
#include "util/stringify.hpp"

Tree::Tree( const std::string &tree_file, 
            const MSA &msa,
            raxml::Model &model, 
            const Options& options)
  : ref_msa_(msa)
  , model_(model)
  , options_(options)
{
  tree_ = utree_ptr(build_tree_from_file(tree_file, nums_), utree_destroy);
  partition_ = partition_ptr( build_partition_from_file(model_, 
                                                        nums_, 
                                                        ref_msa_.num_sites(),
                                                        options_.repeats ), 
                              pll_partition_destroy);

  locks_ = Mutex_List(partition_->tips + partition_->clv_buffers);

  link_tree_msa(tree_.get(), 
                partition_.get(), 
                model_, 
                ref_msa_, 
                nums_.tip_nodes);

  set_unique_clv_indices(get_root(tree_.get()), nums_.tip_nodes);

  LOG_DBG << model_;
  LOG_DBG << "Tree length: " << sum_branch_lengths(tree_.get());

  precompute_clvs(tree_.get(), partition_.get(), nums_);

  LOG_DBG << "Reference tree log-likelihood: "
          << std::to_string(this->ref_tree_logl());

}

/**
  Constructs the structures from binary file.
*/
Tree::Tree( const std::string& bin_file, 
            raxml::Model& model, 
            const Options& options)
  : model_(model)
  , options_(options)
  , binary_(bin_file)
{
  options_.load_binary_mode = true;
  partition_ = partition_ptr(binary_.load_partition(), pll_partition_destroy);
  nums_ = Tree_Numbers(partition_->tips);
  tree_ = utree_ptr(binary_.load_utree(partition_->tips), utree_destroy);
  locks_ = Mutex_List(partition_->tips + partition_->clv_buffers);

  raxml::assign(model_, partition_.get());
  LOG_DBG << model_;
  LOG_DBG << "Tree length: " << sum_branch_lengths(tree_.get());
  LOG_DBG << "Reference tree log-likelihood: "
        << std::to_string(this->ref_tree_logl());
}

/**
  Returns a pointer either to the CLV or tipchar buffer, depending on the index.
  If they are not currently in memory, fetches them from file.
  Ensures that associated scalers are allocated and ready on return.
*/
void* Tree::get_clv(const pll_unode_t* node)
{
  const auto i = node->clv_index;

  // prevent race condition from concurrent access to this function
  Scoped_Mutex lock_by_clv_id(locks_[i]);

  const auto scaler = node->scaler_index;
  const bool use_tipchars = partition_->attributes & PLL_ATTRIB_PATTERN_TIP;

  if (i >= partition_->tips + partition_->clv_buffers) {
    throw std::runtime_error{"Node index out of bounds"};
  }

  void* clv_ptr = nullptr;
  if (use_tipchars and i < partition_->tips) {
    clv_ptr = partition_->tipchars[i];
    // dynamically load from disk if not in memory
    if (options_.load_binary_mode 
        and clv_ptr == nullptr) {
      binary_.load_tipchars(partition_.get(), i);
      clv_ptr = partition_->tipchars[i];
    }
  } else {
    clv_ptr = partition_->clv[i];
    // dynamically load from disk if not in memory
    if (options_.load_binary_mode 
        and clv_ptr == nullptr) {
      binary_.load_clv(partition_.get(), i);
      clv_ptr = partition_->clv[i];
    }
  }

  // dynamically load the scaler if needed
  if (options_.load_binary_mode 
      and scaler != PLL_SCALE_BUFFER_NONE 
      and partition_->scale_buffer[scaler] == nullptr) {
    binary_.load_scaler(partition_.get(), scaler);
  }

  return clv_ptr;
}

double Tree::ref_tree_logl()
{
  std::vector<unsigned int> param_indices(partition_->rate_cats, 0);
  const auto root = get_root(tree_.get());
  // ensure clvs are there
  this->get_clv(root);
  this->get_clv(root->back);

  return pll_compute_edge_loglikelihood(partition_.get(),
                                        root->clv_index,
                                        root->scaler_index,
                                        root->back->clv_index,
                                        root->back->scaler_index,
                                        root->pmatrix_index,
                                        &param_indices[0],
                                        nullptr);
}
