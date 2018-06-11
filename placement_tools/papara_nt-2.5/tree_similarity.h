/*
 * Copyright (C) 2009-2012 Simon A. Berger
 * 
 * This file is part of papara.
 * 
 *  papara is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  papara is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with papara.  If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef __tree_similarity_h
#define __tree_similarity_h
#if 0

#include <vector>
#include <boost/dynamic_bitset_fwd.hpp>

namespace ivy_mike {
    namespace tree_parser_ms {
        class lnode;
    }
} 
   
   
typedef std::vector<boost::dynamic_bitset<> > split_set_t;

bool split_sets_equal( const split_set_t &s1, const split_set_t &s2 );
double compare_trees( ivy_mike::tree_parser_ms::lnode *t1, ivy_mike::tree_parser_ms::lnode *t2, split_set_t &t2_splits );

void get_all_splits( ivy_mike::tree_parser_ms::lnode *t, std::vector< std::pair< ivy_mike::tree_parser_ms::lnode*, ivy_mike::tree_parser_ms::lnode* > > &edges, std::vector<boost::dynamic_bitset<> > &splits, std::vector<ivy_mike::tree_parser_ms::lnode *> &sorted_tips );

#endif

#endif
