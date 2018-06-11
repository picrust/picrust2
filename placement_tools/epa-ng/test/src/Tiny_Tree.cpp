#include "Epatest.hpp"

#include "core/pll/pllhead.hpp"
#include "core/pll/pll_util.hpp"
#include "core/pll/epa_pll_util.hpp"
#include "io/file_io.hpp"
#include "io/Binary.hpp"
#include "tree/Tree_Numbers.hpp"
#include "tree/Tiny_Tree.hpp"
#include "tree/Tree.hpp"
#include "sample/Sample.hpp"
#include "seq/MSA.hpp"
#include "set_manipulators.hpp"
#include "core/raxml/Model.hpp"
#include "core/Lookup_Store.hpp"

#include <tuple>
#include <limits>

using namespace std;

static void place_(const Options options) 
{
  // buildup
  auto msa = build_MSA_from_file(env->reference_file, MSA_Info(env->reference_file), options.premasking);
  auto queries = build_MSA_from_file(env->query_file, MSA_Info(env->query_file), options.premasking);

  auto ref_tree = Tree(env->tree_file, msa, env->model, options);

  shared_ptr<Lookup_Store> lu_ptr(new Lookup_Store(ref_tree.nums().branches, ref_tree.partition()->states));

  auto root = get_root(ref_tree.tree());

  // tests
  Tiny_Tree tt(root, 0, ref_tree, options.prescoring, options, lu_ptr);

  for (auto const &x : queries)
  {
    auto place = tt.place(x);
    auto brlen = root->length;
    EXPECT_NE(place.likelihood(), 0.0);
    EXPECT_NE(place.likelihood(), std::numeric_limits<double>::infinity());
    EXPECT_NE(place.likelihood(), -std::numeric_limits<double>::infinity());
    EXPECT_GT(place.distal_length(), 0.0);
    EXPECT_GT(brlen, place.distal_length());
    EXPECT_GT(place.pendant_length(), 0.0);
  }
  // teardown
}

TEST(Tiny_Tree, place)
{
  all_combinations(place_);
}

static void compare_samples(Sample<>& orig_samp, Sample<>& read_samp, bool verbose=false, unsigned int head=0)
{
  for (size_t seq_id = 0; seq_id < read_samp.size(); ++seq_id) {
    auto& orig_pq = orig_samp[seq_id];
    auto& read_pq = read_samp[seq_id];

    ASSERT_EQ(orig_pq.size(), read_pq.size());

    for (size_t branch_id = 0; branch_id < read_pq.size(); ++branch_id) {
      auto& orig_p = orig_pq[branch_id];
      auto& read_p = read_pq[branch_id];

      if (verbose) {
        auto limit = head ? head : 2000;
        if (branch_id <= limit and seq_id < 1) {
          printf("branch %lu: %f vs branch %lu: %f\n",orig_p.branch_id(), 
                                                    orig_p.likelihood(),
                                                    read_p.branch_id(),  
                                                    read_p.likelihood());
        }
      }

      EXPECT_DOUBLE_EQ(orig_p.likelihood(), read_p.likelihood());
      EXPECT_DOUBLE_EQ(orig_p.pendant_length(), read_p.pendant_length());
      EXPECT_DOUBLE_EQ(orig_p.distal_length(), read_p.distal_length());
      EXPECT_DOUBLE_EQ(orig_p.lwr(), read_p.lwr());
      EXPECT_EQ(orig_p.branch_id(), read_p.branch_id());

    }
  }
}

static void place_from_binary(const Options options)
{
  // setup
  auto tree_file = env->tree_file;
  auto reference_file = env->reference_file;
  auto msa = build_MSA_from_file(env->reference_file, MSA_Info(env->reference_file), options.premasking);
  auto queries = build_MSA_from_file(env->query_file, MSA_Info(env->query_file), options.premasking);
  
  raxml::Model model;

  Tree original_tree(tree_file, msa, model, options);

  dump_to_binary(original_tree, env->binary_file);
  Tree read_tree(env->binary_file, model, options);
  string invocation("./this --is -a test");

  auto orig_lup = std::make_shared<Lookup_Store>( original_tree.nums().branches, 
                                                  original_tree.partition()->states);
  auto read_lup = std::make_shared<Lookup_Store>( read_tree.nums().branches, 
                                                  read_tree.partition()->states);

  if (options.repeats) {
    ASSERT_TRUE(original_tree.partition()->attributes &
                PLL_ATTRIB_SITE_REPEATS);
    ASSERT_TRUE(read_tree.partition()->attributes & PLL_ATTRIB_SITE_REPEATS);
  }

  ASSERT_EQ(original_tree.nums().branches, read_tree.nums().branches);

  vector<pll_unode_t *> original_branches(original_tree.nums().branches);
  vector<pll_unode_t *> read_branches(read_tree.nums().branches);

  auto original_traversed = utree_query_branches( original_tree.tree(), 
                                                  &original_branches[0]);
  auto read_traversed = utree_query_branches( read_tree.tree(), 
                                              &read_branches[0]);

  ASSERT_EQ(original_traversed, read_traversed);
  ASSERT_EQ(original_traversed, original_tree.nums().branches);

  Sample<> orig_samp;
  Sample<> read_samp;

  // test
  for (size_t i = 0; i < original_traversed; i++) {

    Tiny_Tree original_tiny(original_branches[i], 
                            i, 
                            original_tree, 
                            !options.prescoring, 
                            options, 
                            orig_lup);
    Tiny_Tree read_tiny(read_branches[i], 
                        i, 
                        read_tree, 
                        !options.prescoring, 
                        options, 
                        read_lup);

    size_t seq_id = 0;
    for (auto &seq : queries) {
      auto orig_place = original_tiny.place(seq);
      auto read_place     = read_tiny.place(seq);

      ASSERT_DOUBLE_EQ(orig_place.likelihood(), read_place.likelihood());

      orig_samp.add_placement(seq_id, "", orig_place);
      read_samp.add_placement(seq_id, "",read_place);
      
      ++seq_id;
    }
  }

  ASSERT_EQ(orig_samp.size(), read_samp.size());

  compare_samples(orig_samp, read_samp);

  // effectively, do the candidate selection
  compute_and_set_lwr(orig_samp);
  compute_and_set_lwr(read_samp);
  compare_samples(orig_samp, read_samp);

  discard_by_accumulated_threshold(orig_samp, options.prescoring_threshold);
  discard_by_accumulated_threshold(read_samp, options.prescoring_threshold);
  compare_samples(orig_samp, read_samp);

  Work read_work(read_samp);
  Work orig_work(orig_samp);

  // printf("%lu vs %lu\n", orig_work.size(), read_work.size());

  EXPECT_EQ(orig_work.size(), read_work.size());

  // teardown
}

TEST(Tiny_Tree, place_from_binary)
{
  all_combinations(place_from_binary);
  // Options o;
  // o.prescoring = true;
  // o.opt_model = o.opt_branches = true;
  // o.repeats = true;
  // place_from_binary(o);
}
