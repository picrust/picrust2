#pragma once

#include <string>
#include <vector>

#include "core/pll/pllhead.hpp"
#include "tree/Tree_Numbers.hpp"
#include "seq/MSA.hpp"
#include "seq/MSA_Stream.hpp"
#include "seq/Sequence.hpp"
#include "core/raxml/Model.hpp"
#include "tree/Tree.hpp"

void link_tree_msa( pll_utree_t * tree, 
                    pll_partition_t * partition, 
                    raxml::Model& model, 
                    const MSA& msa, 
                    const unsigned int num_tip_nodes);
void precompute_clvs( pll_utree_t const * const tree, 
                      pll_partition_t * partition, 
                      const Tree_Numbers& nums);
void split_combined_msa(MSA& source, 
                        MSA& target, 
                        Tree& tree);
raxml::Model get_model(pll_partition_t* partition);

// operator overloads
bool operator==(const pll_unode_t * node, const Sequence& s);
bool operator==(const Sequence& s, const pll_unode_t * node);
