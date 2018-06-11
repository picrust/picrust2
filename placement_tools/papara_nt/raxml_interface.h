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


#ifndef __raxml_interface_h
#define __raxml_interface_h


#include "ivymike/tree_parser.h"
#include <string>
#include <vector>
#include <map>
#include <stdint.h>
#include <boost/array.hpp>
#include <boost/numeric/ublas/fwd.hpp>

#include <iosfwd>

void optimize_branch_lengths( ivy_mike::tree_parser_ms::lnode *tree, const std::map<std::string, const std::vector<uint8_t> * const> &name_to_seq );
ivy_mike::tree_parser_ms::lnode *optimize_branch_lengths2( ivy_mike::tree_parser_ms::lnode *tree, const std::map<std::string, const std::vector<uint8_t> * const> &name_to_seq, ivy_mike::tree_parser_ms::ln_pool &pool );



std::vector<boost::numeric::ublas::matrix<double> > read_binary_anc_probs( std::istream &pis );

//ivy_mike::tree_parser_ms::lnode *generate_marginal_ancestral_state_pvecs( ivy_mike::tree_parser_ms::ln_pool &pool, const std::string &tree_name, const std::string &ali_name, std::vector<boost::array<std::vector<double>, 4> > *pvecs );
ivy_mike::tree_parser_ms::lnode *generate_marginal_ancestral_state_pvecs( ivy_mike::tree_parser_ms::ln_pool &pool, const std::string &tree_name, const std::string &ali_name, std::vector<boost::numeric::ublas::matrix<double> > *pvecs );

#endif
