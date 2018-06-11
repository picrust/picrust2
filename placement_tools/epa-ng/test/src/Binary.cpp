#include "Epatest.hpp"

#include <vector>

#include "tree/Tree.hpp"
#include "io/Binary.hpp"
#include "io/file_io.hpp"
#include "io/msa_reader.hpp"
#include "util/Options.hpp"
#include "core/raxml/Model.hpp"

using namespace std;

static void check_equal_current(pll_unode_t const * const a, pll_unode_t const * const b)
{
  EXPECT_DOUBLE_EQ(a->length, b->length);
  EXPECT_EQ(a->node_index, b->node_index);
  EXPECT_EQ(a->clv_index, b->clv_index);
  EXPECT_EQ(a->scaler_index, b->scaler_index);
  EXPECT_EQ(a->pmatrix_index, b->pmatrix_index);
  EXPECT_STREQ(a->label, a->label);
}

static void check_equal(pll_unode_t const * const a, pll_unode_t const * const b)
{
  check_equal_current(a, b);

  if (a->next) {
    ASSERT_TRUE(b->next != nullptr);
    check_equal_current(a->next, b->next);
    check_equal_current(a->next->next, b->next->next);
  }
}

static void write_(const Options options)
{
  // setup
  auto msa = build_MSA_from_file(env->reference_file, MSA_Info(env->reference_file), options.premasking);
  raxml::Model model;

  Tree tree(env->tree_file, msa, model, options);

  // test
  dump_to_binary(tree, env->binary_file);
}

TEST(Binary, write)
{
  all_combinations(write_);
}

static int full_trav(pll_unode_t*)
{
  return 1;
}

static auto create_scaler_to_clv_map(Tree& tree)
{
  const auto num_scalers = tree.partition()->scale_buffers;
  std::vector<unsigned int> map(num_scalers);

  std::vector<pll_unode_t*> travbuffer(tree.nums().nodes);
  unsigned int trav_size = 0;
  pll_utree_traverse( get_root(tree.tree()),
                      PLL_TREE_TRAVERSE_POSTORDER,
                      full_trav, 
                      &travbuffer[0], 
                      &trav_size);

  for (auto& n : travbuffer) {
    if (n->scaler_index != PLL_SCALE_BUFFER_NONE) {
      map[n->scaler_index] = n->clv_index;
    }
  }

  return map;
}

static double loglh(pll_partition* partition, pll_unode_t* node)
{
  std::vector<unsigned int> param_indices(partition->rate_cats, 0);
  return pll_compute_edge_loglikelihood(partition,
                                        node->clv_index,
                                        node->scaler_index,
                                        node->back->clv_index,
                                        node->back->scaler_index,
                                        node->pmatrix_index,
                                        &param_indices[0],
                                        nullptr);
}

static void read_(Options options)
{
  // setup
  auto msa = build_MSA_from_file(env->reference_file, MSA_Info(env->reference_file), options.premasking);
  raxml::Model model;
  // double freqs[4] = {0.1,0.2,0.3,0.4};
  // double alpha = 42.42;
  // double subs[6] = {0.1,0.2,0.3,0.4,0.5,0.6};
  // int symm[6] = {1,2,3,4,5,6};
  // model.base_frequencies(freqs, 4);
  // model.substitution_rates(subs, 6);
  // model.symmetries(symm, 6);
  // model.alpha(alpha);
  Tree original_tree(env->tree_file, msa, model, options);
  dump_to_binary(original_tree, env->binary_file);

  // test
  Tree read_tree(env->binary_file, model, options);

  auto part = original_tree.partition();
  auto read_part = read_tree.partition();

  // compare numbered jplace strings
  string original_nns(get_numbered_newick_string(original_tree.tree()));
  string read_nns(get_numbered_newick_string(read_tree.tree()));

  EXPECT_STREQ(original_nns.c_str(), read_nns.c_str());
  // compare tree traversals
  ASSERT_EQ(original_tree.nums().nodes, read_tree.nums().nodes);
  vector<pll_unode_t *> original_nodes(original_tree.nums().nodes);
  vector<pll_unode_t *> read_nodes(read_tree.nums().nodes);
  unsigned int original_traversed, read_traversed;
  pll_utree_traverse( get_root(original_tree.tree()),
                      PLL_TREE_TRAVERSE_POSTORDER,
                      cb_full_traversal, 
                      &original_nodes[0], 
                      &original_traversed);
  pll_utree_traverse( get_root(read_tree.tree()),
                      PLL_TREE_TRAVERSE_POSTORDER,
                      cb_full_traversal,
                      &read_nodes[0],
                      &read_traversed);

  ASSERT_EQ(original_traversed, read_traversed);
  ASSERT_EQ(original_traversed, original_tree.nums().nodes);

  for (size_t i = 0; i < read_traversed; i++) {
    auto o = original_nodes[i];
    auto r = read_nodes[i];
    // printf("orig: %d back: %d\n", o->clv_index, o->back->clv_index);
    // printf("read: %d back: %d\n", r->clv_index, r->back->clv_index);
    check_equal(o, r);
  }


  ASSERT_EQ(part->sites, read_part->sites);
  ASSERT_EQ(part->states, read_part->states);
  ASSERT_EQ(part->states_padded, read_part->states_padded);
  ASSERT_EQ(part->rate_cats, read_part->rate_cats);
  ASSERT_EQ(part->tips, read_part->tips);
  ASSERT_EQ(part->clv_buffers, read_part->clv_buffers);
  ASSERT_EQ(part->attributes, read_part->attributes);

  auto read_freqs = read_part->frequencies[0];
  auto freqs = part->frequencies[0];
  for (size_t i = 0; i < 4; i++) {
    EXPECT_DOUBLE_EQ(freqs[i], read_freqs[i]);
  }

  auto read_subs = read_part->subst_params[0];
  auto subs = part->subst_params[0];
  for(size_t i = 0; i < 6; i++) {
    EXPECT_DOUBLE_EQ(subs[i], read_subs[i]);
  }

  auto read_rates = read_part->rates;
  auto rates = part->rates;
  for(size_t i = 0; i < part->rate_cats; i++) {
    EXPECT_DOUBLE_EQ(rates[i], read_rates[i]);
  }

  auto read_rate_weights = read_part->rate_weights;
  auto rate_weights = part->rate_weights;
  for(size_t i = 0; i < part->rate_cats; i++) {
    EXPECT_DOUBLE_EQ(rate_weights[i], read_rate_weights[i]);
  }



  // compare tips
  if (read_part->attributes & PLL_ATTRIB_PATTERN_TIP) {
    for (size_t i = 0; i < part->tips; i++) {
      pll_unode_t node;
      node.clv_index = i;
      node.scaler_index = 0;
      auto read_tipchars = static_cast<char*>(read_tree.get_clv(&node));
      for (size_t j = 0; j < part->sites; j++) {
        EXPECT_EQ(part->tipchars[i][j], read_tipchars[j]);
      }
    }
  }

  // compare clvs
  size_t start = (read_part->attributes & PLL_ATTRIB_PATTERN_TIP) ? part->tips : 0;

  for (size_t i = start; i < part->tips + part->clv_buffers; i++) {
    pll_unode_t node;
    node.clv_index = i;
    node.scaler_index = 0;
    const auto read_clv = static_cast<double*>(read_tree.get_clv(&node));
    const auto clv_size = pll_get_clv_size(part, i);
    for (size_t j = 0; j < clv_size; j++) {
      ASSERT_DOUBLE_EQ(part->clv[i][j], read_clv[j]);
    }
  }

  const auto stoc = create_scaler_to_clv_map(read_tree);

  // check scalers
  for (size_t i = 0; i < part->scale_buffers; i++) {
    pll_unode_t node;
    node.clv_index = 0;
    node.scaler_index = i;
    read_tree.get_clv(&node);
    const auto scaler_size = pll_get_sites_number(read_part, stoc[i]);
    for (size_t j = 0; j < scaler_size; j++) {
      // printf("%u v %u\n",part->scale_buffer[i][j], read_part->scale_buffer[i][j] );
      EXPECT_EQ(part->scale_buffer[i][j], read_part->scale_buffer[i][j]);
    }
  }

  // check repeats specific structures
  const auto rep = part->repeats;
  const auto read_rep = read_part->repeats;
  if (rep) {
    for (size_t i = 0; i < part->scale_buffers; ++i) {
      EXPECT_EQ(rep->perscale_ids[i], read_rep->perscale_ids[i]);
    }
    for (size_t i = 0; i < part->clv_buffers + part->tips; ++i) {
      EXPECT_EQ(rep->pernode_ids[i], read_rep->pernode_ids[i]);
      EXPECT_EQ(rep->pernode_allocated_clvs[i],
                read_rep->pernode_allocated_clvs[i]);

      for (size_t j = 0; j < part->sites; ++j) {
        EXPECT_EQ(rep->pernode_site_id[i][j], read_rep->pernode_site_id[i][j]);
      }

      for (size_t j = 0; j < pll_get_sites_number(part, i); ++j) {
        EXPECT_EQ(rep->pernode_id_site[i][j], read_rep->pernode_id_site[i][j]);
      }
    }
  }

  // finally, compare some likelihoods
  for (size_t i = 0; i < read_traversed; i++) {
    auto o = original_nodes[i];
    auto r = read_nodes[i];
    auto original_logl = loglh(part, o);
    auto read_logl = loglh(read_part, r);
    // printf("%f vs %f\n", original_logl, read_logl);
    EXPECT_DOUBLE_EQ(original_logl, read_logl);
  }
}

TEST(Binary, read)
{
  all_combinations(read_);
}
