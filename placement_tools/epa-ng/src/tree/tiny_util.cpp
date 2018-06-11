#include "tree/tiny_util.hpp"

#include <type_traits>

#include "core/pll/pll_util.hpp"


constexpr unsigned int proximal_clv_index         = 4;
constexpr unsigned int inner_clv_index            = 3;
constexpr unsigned int new_tip_clv_index          = 1;
constexpr unsigned int distal_clv_index_if_tip    = 2;
constexpr unsigned int distal_clv_index_if_inner  = 5;

template <class T, 
          typename = typename std::enable_if<std::is_pointer<T>::value>::type>
static void alloc_and_copy(T& dest, const T src, const size_t size)
{
  using base_t = std::remove_pointer_t<decltype(src)>;

  if (dest != nullptr) {
    free(dest);
  }

  dest = static_cast<T>(calloc(size, sizeof(base_t)));

  memcpy( dest,
          src,
          size * sizeof(base_t));
}

//TODO src_part should be const const
static void deep_copy_scaler( pll_partition_t* dest_part,
                              pll_unode_t* dest_node,
                              pll_partition_t const * const src_part,
                              pll_unode_t const * const src_node)
{
  if (src_node->scaler_index != PLL_SCALE_BUFFER_NONE
    and src_part->scale_buffer[src_node->scaler_index] != nullptr) {
    
    const auto scaler_size = 
    pll_get_sites_number( const_cast<pll_partition_t*>(src_part), 
                          src_node->clv_index);

    alloc_and_copy( dest_part->scale_buffer[dest_node->scaler_index], 
                    src_part->scale_buffer[src_node->scaler_index], 
                    scaler_size);
  }
}

static void deep_copy_repeats(pll_partition_t* dest_part,
                              pll_unode_t* dest_node,
                              pll_partition_t const * const src_part,
                              pll_unode_t const * const src_node)
{
  // copy size info
  if (src_node->scaler_index != PLL_SCALE_BUFFER_NONE) {
    dest_part->repeats->perscale_ids[dest_node->scaler_index]
      = src_part->repeats->perscale_ids[src_node->scaler_index];
  }
  
  dest_part->repeats->pernode_ids[dest_node->clv_index]
    = src_part->repeats->pernode_ids[src_node->clv_index];
  dest_part->repeats->pernode_allocated_clvs[dest_node->clv_index]
    = src_part->repeats->pernode_allocated_clvs[src_node->clv_index];

  if (src_part->repeats->pernode_ids[src_node->clv_index]) {
    const auto size = 
    pll_get_sites_number( const_cast<pll_partition_t*>(src_part), 
                          src_node->clv_index);
    
    alloc_and_copy( dest_part->repeats->pernode_site_id[dest_node->clv_index], 
                    src_part->repeats->pernode_site_id[src_node->clv_index], 
                    src_part->sites);
    alloc_and_copy( dest_part->repeats->pernode_id_site[dest_node->clv_index], 
                    src_part->repeats->pernode_id_site[src_node->clv_index], 
                    size);
  }

}


pll_partition_t * make_tiny_partition(Tree& reference_tree, 
                                      const pll_utree_t * tree, 
                                      pll_unode_t const * const old_proximal, 
                                      pll_unode_t const * const old_distal, 
                                      const bool tip_tip_case)
{
  /**
    As we work with PLL_PATTERN_TIP functionality, special care has to be taken in regards to the node and partition
    structure: PLL assumes that any node with clv index < number of tips is in fact a real tip, that is
    a tip that uses a character array instead of a real clv. Here we need to set up the node/partition to fool pll:
    the tips that actually contain CLVs copied over from the reference node have their index set to greater than
    number of tips. This results in a acceptable amount of wasted memory that is never used (num_sites * bytes
    * number of clv-tips)
  */
  pll_partition_t const * const old_partition = reference_tree.partition();
  assert(old_partition);

  bool use_tipchars = old_partition->attributes & PLL_ATTRIB_PATTERN_TIP;

  // tip_inner case: both reference nodes are inner nodes
  // tip tip case: one for the "proximal" clv tip
  const unsigned int num_clv_tips = tip_tip_case ? 1 : 2;

  auto proximal = tree->nodes[0];
  auto distal = tree->nodes[1];

  pll_partition_t * tiny = pll_partition_create(
    3, // tips
    1 + num_clv_tips, // extra clv's
    old_partition->states, old_partition->sites,
    old_partition->rate_matrices,
    3, // number of prob. matrices (one per possible unique branch length)
    old_partition->rate_cats,
    3, // number of scale buffers (one per possible inner node)
    old_partition->attributes);

  assert(tiny);

  unsigned int i;
  free(tiny->rates);
  tiny->rates = old_partition->rates;
  if (tiny->subst_params) {
    for (i = 0; i < tiny->rate_matrices; ++i) {
      pll_aligned_free(tiny->subst_params[i]);
    }
  }
  free(tiny->subst_params);
  tiny->subst_params = old_partition->subst_params;
  if (tiny->frequencies) {
    for (i = 0; i < tiny->rate_matrices; ++i) {
      pll_aligned_free(tiny->frequencies[i]);
    }
  }
  free(tiny->frequencies);
  tiny->frequencies = old_partition->frequencies;
  if (tiny->eigenvecs) {
    for (i = 0; i < tiny->rate_matrices; ++i) {
      pll_aligned_free(tiny->eigenvecs[i]);
    }
  }
  free(tiny->eigenvecs);
  tiny->eigenvecs = old_partition->eigenvecs;
  if (tiny->inv_eigenvecs) {
    for (i = 0; i < tiny->rate_matrices; ++i) {
      pll_aligned_free(tiny->inv_eigenvecs[i]);
    }
  }
  free(tiny->inv_eigenvecs);
  tiny->inv_eigenvecs = old_partition->inv_eigenvecs;
  if (tiny->eigenvals) {
    for (i = 0; i < tiny->rate_matrices; ++i) {
      pll_aligned_free(tiny->eigenvals[i]);
    }
  }
  free(tiny->eigenvals);
  tiny->eigenvals = old_partition->eigenvals;
  if (tiny->prop_invar) {
    free(tiny->prop_invar);
  }
  tiny->prop_invar = old_partition->prop_invar;
  free(tiny->eigen_decomp_valid);
  tiny->eigen_decomp_valid = old_partition->eigen_decomp_valid;
  if (tiny->pattern_weights) {
    free(tiny->pattern_weights);
  }
  tiny->pattern_weights = old_partition->pattern_weights;

  // shallow copy major buffers
  pll_aligned_free(tiny->clv[proximal->clv_index]);
  tiny->clv[proximal->clv_index] =
    static_cast<double*>(reference_tree.get_clv(old_proximal));


  if(tip_tip_case and use_tipchars) {
    std::string sequence(tiny->sites, 'A');
    if( pll_set_tip_states(tiny, distal->clv_index, get_char_map(old_partition), sequence.c_str())
        == PLL_FAILURE) {
      throw std::runtime_error{"Error setting tip state"};
    }
    pll_aligned_free(tiny->tipchars[distal->clv_index]);
    tiny->tipchars[distal->clv_index] = static_cast<unsigned char*>(reference_tree.get_clv(old_distal));
  } else {
    pll_aligned_free(tiny->clv[distal->clv_index]);
    tiny->clv[distal->clv_index] = static_cast<double*>(reference_tree.get_clv(old_distal));
  }


  // deep copy scalers
  deep_copy_scaler( tiny,
                    proximal,
                    old_partition,
                    old_proximal);

  deep_copy_scaler( tiny,
                    distal,
                    old_partition,
                    old_distal);

  // copy the repeats structures
  if (old_partition->repeats) {
    // then do the per-clv stuff, but only for the two relevant clv
    deep_copy_repeats(tiny,
                      proximal,
                      old_partition,
                      old_proximal);

    deep_copy_repeats(tiny,
                      distal,
                      old_partition,
                      old_distal);

    pll_resize_repeats_lookup(tiny, tiny->sites * tiny->states);
  }

  return tiny;
}

void tiny_partition_destroy(pll_partition_t * partition)
{
  if (partition) {
    // unset shallow copied things
    partition->rates              = nullptr;
    partition->subst_params       = nullptr;
    partition->frequencies        = nullptr;
    partition->eigenvecs          = nullptr;
    partition->inv_eigenvecs      = nullptr;
    partition->eigenvals          = nullptr;
    partition->prop_invar         = nullptr;
    partition->eigen_decomp_valid = nullptr;
    partition->pattern_weights    = nullptr;

    partition->clv[proximal_clv_index] = nullptr;

    const bool distal_is_tip    = partition->clv_buffers == 3 ? false : true;
    const bool pattern_tip_mode = partition->attributes & PLL_ATTRIB_PATTERN_TIP;

    if (distal_is_tip) {
      if (pattern_tip_mode) {
        partition->tipchars[distal_clv_index_if_tip] = nullptr;
      } else {
        partition->clv[distal_clv_index_if_tip] = nullptr;
      }
    } else {
      partition->clv[distal_clv_index_if_inner] = nullptr;
    }

    pll_partition_destroy(partition);
  }
}

pll_utree_t * make_tiny_tree_structure( const pll_unode_t * old_proximal, 
                                        const pll_unode_t * old_distal,
                                        const bool tip_tip_case)
{
  const unsigned int inner_scaler_index = 1;
  const unsigned int proximal_scaler_index = 0;
  const unsigned int distal_scaler_index = 2;

  /**
    As we work with PLL_PATTERN_TIP functionality, special care has to be taken in regards to the tree and partition
    structure: PLL assumes that any node with clv index < number of tips is in fact a real tip, that is
    a tip that uses a character array instead of a real clv. Here we need to set up the tree/partition to fool pll:
    the tips that actually contain CLVs copied over from the reference tree have their index set to greater than
    number of tips. This results in a acceptable amount of wasted memory that is never used (num_sites * bytes
    * number of clv-tips)
  */
  // if tip-inner case
  unsigned int distal_clv_index = (tip_tip_case) ? distal_clv_index_if_tip : distal_clv_index_if_inner;

  auto inner = static_cast<pll_unode_t *>(calloc(1,sizeof(pll_unode_t)));
  inner->next = static_cast<pll_unode_t *>(calloc(1,sizeof(pll_unode_t)));
  inner->next->next = static_cast<pll_unode_t *>(calloc(1,sizeof(pll_unode_t)));
  inner->next->next->next = inner;

  auto new_tip = static_cast<pll_unode_t *>(calloc(1,sizeof(pll_unode_t)));
  auto proximal = static_cast<pll_unode_t *>(calloc(1,sizeof(pll_unode_t)));
  auto distal = static_cast<pll_unode_t *>(calloc(1,sizeof(pll_unode_t)));

  // connect the nodes to each other
  inner->back = new_tip;
  new_tip->back = inner;
  inner->next->back = distal;
  distal->back = inner->next;
  inner->next->next->back = proximal;
  proximal->back = inner->next->next;

  // set up clv indices
  inner->clv_index = inner_clv_index;
  inner->next->clv_index = inner_clv_index;
  inner->next->next->clv_index = inner_clv_index;
  proximal->clv_index = proximal_clv_index;
  distal->clv_index = distal_clv_index;
  new_tip->clv_index = new_tip_clv_index;

  // set up scaler indices
  new_tip->scaler_index = PLL_SCALE_BUFFER_NONE;
  inner->scaler_index = inner_scaler_index;
  inner->next->scaler_index = inner_scaler_index;
  inner->next->next->scaler_index = inner_scaler_index;
  proximal->scaler_index = (old_proximal->scaler_index == PLL_SCALE_BUFFER_NONE) ?
    PLL_SCALE_BUFFER_NONE : proximal_scaler_index;
  distal->scaler_index = (old_distal->scaler_index == PLL_SCALE_BUFFER_NONE) ?
    PLL_SCALE_BUFFER_NONE : distal_scaler_index;


  reset_triplet_lengths(inner, nullptr, old_distal->length);

  auto tree = static_cast<pll_utree_t*>(calloc(1, sizeof(pll_utree_t)));

  tree->edge_count = 3;
  tree->tip_count = 3;
  tree->inner_count = 1;

  tree->nodes = static_cast<pll_unode_t **>(
    calloc( tree->inner_count + tree->tip_count,
            sizeof(pll_unode_t*)));

  tree->nodes[0] = proximal;
  tree->nodes[1] = distal;
  tree->nodes[2] = new_tip;
  tree->nodes[3] = inner;

  return tree;
}
