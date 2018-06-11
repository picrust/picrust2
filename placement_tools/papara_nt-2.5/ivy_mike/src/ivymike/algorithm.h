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


#ifndef __ivy_mike__algorithm_h
#define __ivy_mike__algorithm_h

#include <algorithm>
#include <functional>
#include <vector>

namespace ivy_mike {
// the mighty twizzle algorithm (uhm, wouldn't the name binary_transform be more appropriate?)
// NOTE TO SELF: dumbass, std::transform works for binary ops already...
template<typename iiter1_, typename iiter2_, typename oiter_, typename function_>
function_ binary_twizzle( iiter1_ first1, iiter1_ last1, iiter2_ first2, oiter_ res, function_ func ) {
    for( ; first1 != last1; ++first1, ++first2, ++res ) {
        *res = func( *first1, *first2 );
    }
    return func;
}





template<typename iiter1_, typename iiter2_, typename pred_>
size_t binary_count_if( iiter1_ first1, iiter1_ last1, iiter2_ first2, pred_ pred ) {
    size_t count = 0;
    while( first1 != last1 ) {
        if( pred( *first1++, *first2++) ) {
            ++count;
        }
    }
    
    return count;
    
}


// template<typename iiter1_, typename iiter2_, typename oiter_, typename function_>
// size_t binary_transform( iiter1_ first1, iiter1_ last1, iiter2_ first2, oiter_ res, pred_ pred ) {
//     size_t count = 0;
//     while( first1 != last1 ) {
//         if( pred( *first1++, *first2++) ) {
//             ++count;
//         }
//     }
//     
//     return count;
//     
// }

template<typename iiter1_, typename iiter2_>
size_t count_equal( iiter1_ first1, iiter1_ last1, iiter2_ first2 ) {
    return binary_count_if( first1, last1, first2, std::equal_to<typename iiter1_::value_type>() );
}


template<typename T>
  struct scaler_clamp {
  	T m_s;
  	T m_hi, m_lo;
  	scaler_clamp( T s, T lo, T hi ) : m_s(s), m_hi(hi), m_lo(lo) {}

  	T operator()( T v ) {
  		return std::min( m_hi, std::max( m_lo, v * m_s ));;
  	}

  };

// this should basically emulate what push_back( T && ) is doing in C++11. right?
// WARNING: parameter v will be left empty.
template<typename T>
void push_back_swap( std::vector<T> &vec, T & v ) {
	vec.push_back(T());
	v.swap(vec.back());
}


// idea from http://www.dev102.com/2009/01/12/c-tip-how-to-get-array-length/
template<typename T, int size>
int arrlen(T(&)[size]){return size;}

}



#endif
