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


#include <iostream>

#include "ivymike/tree_traversal_utils.h"
#include "ivymike/tree_parser.h"
#include "ivymike/multiple_alignment.h"
#include "boost/dynamic_bitset.hpp"

using ivy_mike::tree_parser_ms::ln_pool;
using ivy_mike::tree_parser_ms::parser;
using ivy_mike::tree_parser_ms::lnode;

using ivy_mike::tip_collector_dumb;
using ivy_mike::visit_lnode;


bool is_gap( char c ) {
    return c == '-' || c == '?' || toupper(c) == 'N';
}

void pad( std::ostream &os, size_t n, char c ) {
    for( size_t i = 0; i < n; ++i ) {
        os << c;
    }
}

class tip_name_collector {

public:
	typedef ivy_mike::tree_parser_ms::lnode lnode;
	tip_name_collector( std::vector<std::string> *names ) : names_(names) {}
	void operator()(lnode *n) {
		if( n->m_data->isTip ) {
			names_->push_back(n->m_data->tipName);
		}
	}
private:
	std::vector<std::string> *names_;
};

int main( int argc, char *argv[] ) {
    if( argc != 3 ) {
        std::cerr << "usage: " << argv[0] << "<tree file> <phylip file>\n";
        return 0;
    }
    
    ln_pool pool;
    
    parser p(argv[1], pool );
    lnode *n = p.parse();
    
    std::vector<std::string> tip_names;
   
    tip_name_collector tnc(&tip_names); 
    visit_lnode( n, tnc );  
//    apply_lnode(n, [&](lnode *x) {
//        if( x->m_data->isTip ) {
//           tip_names.push_back( x->m_data->tipName );
//        }
//    });
    
    ivy_mike::multiple_alignment ma;
    ma.load_phylip(argv[2]);
    
    std::sort( tip_names.begin(), tip_names.end() );
    
    
    boost::dynamic_bitset<> bs_all(ma.data.front().size() );
    //bs_all.flip();
    
    size_t max_name_len = 0;
    for( size_t i = 0, e = ma.names.size(); i != e; ++i ) {
        const std::string &name = ma.names[i];
        max_name_len = std::max( max_name_len, name.size() );
        
        const std::vector<uint8_t> &data = ma.data[i];
        
        if( !std::binary_search( tip_names.begin(), tip_names.end(), name ) ) {
            boost::dynamic_bitset<> bs(data.size());
            
            for( size_t j = 0; j < data.size(); ++j ) {
                bs[j] = !is_gap(data[j]);
            }
            
            
            bs_all |= bs;
        }
    }
    
    std::cout << ma.names.size() << " " << bs_all.count() << "\n";
    
    std::vector<size_t> unmasked;
    {
        size_t i = bs_all.find_first();
        while( i != bs_all.npos ) {
            unmasked.push_back(i);
            i = bs_all.find_next(i);
        }
    }
    for( size_t i = 0, e = ma.names.size(); i != e; ++i ) {
        std::cout << ma.names[i];
        const std::vector<uint8_t> &data = ma.data[i];
        pad( std::cout, max_name_len - ma.names[i].size() + 1, ' ' );
        
        // excract 'unmasked' characters from current sequence
    //    std::transform(unmasked.begin(), unmasked.end(), std::ostream_iterator<char>(std::cout), [&](size_t u) { return data[u]; });
    //
    	for( std::vector<size_t>::iterator it = unmasked.begin(); it != unmasked.end(); ++it ) {
		assert( *it < data.size() );
		std::cout << data[*it];
	}
        std::cout << "\n";
    }
    
}
