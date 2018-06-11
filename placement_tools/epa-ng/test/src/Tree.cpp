#include "Epatest.hpp"

#include "io/file_io.hpp"
#include "io/Binary_Fasta.hpp"
#include "tree/Tree.hpp"
#include "seq/MSA.hpp"
#include "seq/MSA_Info.hpp"
#include "core/place.hpp"

#include <string>
#include <vector>
#include <limits>

using namespace std;

TEST(Tree, process_from_binary)
{
  // setup
  MSA_Info qry_info(env->query_file);
  MSA_Info ref_info(env->reference_file);

  MSA_Info::or_mask(qry_info, ref_info);

  auto queries = Binary_Fasta::fasta_to_bfast(env->query_file, env->out_dir);
  raxml::Model model;
  Options options;
  auto msa = build_MSA_from_file(env->reference_file, ref_info, options.premasking);
  Tree original_tree(env->tree_file, msa, model, options);
  dump_to_binary(original_tree, env->binary_file);
  Tree read_tree(env->binary_file, model, options);

  EXPECT_DOUBLE_EQ(original_tree.ref_tree_logl(), read_tree.ref_tree_logl());

  // test
  string invocation("./this --is -a test");

  simple_mpi(read_tree, queries, qry_info, env->out_dir, options, invocation);

  options.prescoring = true;
  simple_mpi(read_tree, queries, qry_info, env->out_dir, options, invocation);

  Tree mvstree;
  mvstree = Tree(env->binary_file, model, options);

  simple_mpi(mvstree, queries, qry_info, env->out_dir, options, invocation);

  // teardown
}

TEST(Tree, combined_input_file)
{
  auto combined_msa = build_MSA_from_file(env->combined_file, MSA_Info(env->combined_file), true);
  auto tree = Tree(env->tree_file, combined_msa, env->model, env->options);
}

TEST(Tree, rooted_input)
{
  auto msa = build_MSA_from_file(env->reference_file, MSA_Info(env->reference_file), true);
  auto tree = Tree(env->tree_file_rooted, msa, env->model, env->options);
}
