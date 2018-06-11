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

#include <cassert>
#include <cstdlib>
#include <stdint.h>
#include <vector>


#include "align_utils.h"

void align_utils::trace_to_position_map( const std::vector< uint8_t >& gaps, std::vector< int > *map) {
    
    int seq_ptr = 0;
    
    for ( std::vector<uint8_t>::const_reverse_iterator git = gaps.rbegin(); git != gaps.rend(); ++git ) {
        
        if ( *git == 1) {
            ++seq_ptr;
        } else if ( *git == 0 ) {
            
            map->push_back(seq_ptr);
            
            ++seq_ptr;
        } else {
            map->push_back(seq_ptr);
            
        }
    }
}


uint8_t align_utils::decode_dna( int s ) {
    assert( s >= 0 && s < 4 );
    const static uint8_t map[4] = {'A','C','G','T'};
    
    return map[size_t(s)];
}


void align_utils::realize_trace( const std::vector<uint8_t> &seq, const std::vector<uint8_t> &tb, std::vector<uint8_t> *out ) {
    assert( out != 0 );
    
    std::vector<uint8_t>::const_reverse_iterator tb_it;
    std::vector<uint8_t>::const_iterator seq_it;
    for( tb_it = tb.rbegin(), seq_it = seq.begin(); tb_it != tb.rend(); ++tb_it ) {
        if( *tb_it == 0  ) {
            assert( seq_it != seq.end() );
            
            out->push_back(decode_dna(*seq_it));
            ++seq_it;
            
        } else if( *tb_it == 2 ) {
            assert( seq_it != seq.end() );
            ++seq_it;
        } else {
            out->push_back('-');
        }
    }
}
