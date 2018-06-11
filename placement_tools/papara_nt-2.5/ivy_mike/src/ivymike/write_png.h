/*
 * Copyright (C) 2009-2012 Simon A. Berger
 * 
 * This file is part of ivy_mike.
 * 
 *  ivy_mike is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ivy_mike is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with ivy_mike.  If not, see <http://www.gnu.org/licenses/>.
 */


#include <iostream>
#include <algorithm>
#include "ivymike/tdmatrix.h"

#ifndef __ivy_mike__write_png_h
#define __ivy_mike__write_png_h
namespace ivy_mike {

template<typename ma_t>
static void write_png( ma_t &ma, std::ostream &os ) {
    
    
    
    typedef typename ma_t::value_type v_t;
    
    v_t max_value = (*std::max_element( ma.begin(), ma.end() ));
    v_t min_value = (*std::min_element( ma.begin(), ma.end() ));
    
    
    
    float rng = float(max_value-min_value);
    
    os << "P2\n";
    os << ma.size() << " " << ma[0].size() << "\n";
    os << 255 << "\n";
    for( typename ma_t::row_iterator row_it = ma.row_begin(); row_it !=  ma.row_end(); ++row_it ) {
        typename ma_t::row_type row = *row_it;
        
        for( typename ma_t::row_type::iterator col_it = row.begin(); col_it !=  row.end(); ++col_it ) {
            os << int((*col_it - min_value) / rng * 255) << " ";
        }
        os << "\n";
    }
    
}

}

#endif