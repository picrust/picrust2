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


#ifndef __align_utils_h__
#define __align_utils_h__


#include <vector>
#include <stdint.h>


namespace align_utils {
void trace_to_position_map( const std::vector< uint8_t >& gaps, std::vector< int > *map);
uint8_t decode_dna( int s );
void realize_trace( const std::vector<uint8_t> &seq, const std::vector<uint8_t> &tb, std::vector<uint8_t> *out );

}

#endif
