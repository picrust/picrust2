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

#ifndef __dtw_h
#define __dtw_h

#include <fstream>
#include <stdint.h>
#include <vector>
#include <cassert>

#include "vec_unit.h"
#include "aligned_buffer.h"
#include "ivymike/tdmatrix.h"

template<typename vec_iter_t>
static void traceback( ivy_mike::tdmatrix<uint8_t> &tbmat, vec_iter_t a_end, vec_iter_t b_end ) {
    ptrdiff_t a = tbmat.size() - 1;
    ptrdiff_t b = tbmat[0].size() - 1;
    int len = 0;
    std:: ofstream os( "tb.txt" );
    while( a > 0 || b > 0 ) {
        if( tbmat[a][b] == 0 ) {
            a--;
            b--;
            
            a_end--;
            b_end--;
            
            //std::cout << len << " " << 'm' << ' ' << 'm' << '\n'; 
            
        } else if( tbmat[a][b] == 1 ) {
            a--;
            a_end--;
            
        } else if( tbmat[a][b] == 2 ) {
            b--;
            b_end--;
//             std::cout << len << " " << 'm' << ' ' << '-' << '\n'; 
            
        } else {
            
            if( a > 0 ) {
                a--;
                a_end--;
//                 std::cout << len << " " << '-' << ' ' << 'm' << '\n';
                
            } else if( b > 0 ) {
                b--;
                b_end--;
//                 std::cout << len << " " << 'm' << ' ' << '-' << '\n';
                
            } else {
                throw std::runtime_error( "meeeeeeeeeeeeeeeeeep\n" );    
            }
            
        }
        os << len << " " << *a_end << ' ' << *b_end << '\n'; 
        len++;
        
    }
    
}

template<typename value_t, typename vec_iter_t, typename absfunc, typename min3func>
static value_t dtw_align( const vec_iter_t &a_begin, const vec_iter_t &a_end, const vec_iter_t &b_begin, const vec_iter_t &b_end, const value_t large_value, absfunc my_abs, min3func my_min3 ) {
    aligned_buffer<value_t> mat( b_end - b_begin + 1 );
    std::fill( mat.begin(), mat.end(), large_value );
    
#define DO_TB
#ifdef DO_TB    
    ivy_mike::tdmatrix<uint8_t> tbmat(a_end - a_begin + 1, b_end - b_begin + 1);
    std::fill( tbmat.begin(), tbmat.end(), 5 );
#endif    
    
    
    mat[0] = 0;
    
//     std:: ofstream os( "/tmp/yyy.txt" );
    value_t last_s = value_t(); // guru meditation: I thought that the default constructor of POD types is supposed to do nothing. Still this seems like a generic way to get rid of the 'last_s might be uninitialized' warning. Finally I should get a copy of TC++PL...
    assert( a_begin != a_end ); // just to make sure that last_s really is never used uninitialized...
    
    for( vec_iter_t ait = a_begin; ait != a_end; ++ait ) {
        const value_t a = *ait;
        last_s = large_value;

        // for the first row, the diagonal value for the first colum is zero (inialization above)
        value_t diag = mat[0];
        // for the remaining rows it will be infinity. note: mat[0] is never written to in the b-loop
        mat[0] = large_value;
        
        
        value_t *upper = mat(1);
        
        for( vec_iter_t bit = b_begin; bit != b_end; ++bit, ++upper ) {
            const value_t b = *bit;
            
            const value_t sd = diag;
            const value_t su = *upper;
            const value_t sl = last_s;
            
            value_t cost = my_abs( a - b );
            
            //std::cout << "cost: " << cost << "\n";
            
            diag = *upper;
            
            const value_t min = my_min3( sd, su, sl );
            
            if( int(min) + cost > large_value ) {
                std::cout << "overflow\n";
                cost = large_value - min;
            } 
            
            *upper = last_s = cost + min;
            
//             os << "x " << last_s << "\n";
#ifdef DO_TB
            if( sd <= su && sd <= sl ) {
                tbmat[ait - a_begin + 1][bit - b_begin + 1] = 0;
            } else {
                if( su <= sl ) {
                    tbmat[ait - a_begin + 1][bit - b_begin + 1] = 1;
                } else {
                    tbmat[ait - a_begin + 1][bit - b_begin + 1] = 2;
                }
            }
#endif            
        }
        
    }
    
    
#ifdef DO_TB
    traceback( tbmat, a_end, b_end );
#endif
    
    //return mat[mat.size() - 1];
    return last_s;
}

template<typename score_t>
struct dtw_align_ps {
    aligned_buffer<score_t> mat;
    aligned_buffer<score_t> out;
};

template <typename score_t, size_t W,typename vec_iter_t>
void dtw_align_vec( dtw_align_ps<score_t> &ps, aligned_buffer<score_t> &aprofile, size_t asize, vec_iter_t b_begin, vec_iter_t b_end, std::vector<score_t> &out ) {
    typedef vector_unit<score_t,W> vu;
    typedef typename vu::vec_t vec_t;
    
    aligned_buffer<score_t> &mat = ps.mat;
    if( mat.size() < (asize+1) * W ) {
        mat.resize( (asize + 1) * W );
    }
    const score_t LARGE_VALUE = vu::LARGE_VALUE;
    std::fill( mat.begin(), mat.end(), LARGE_VALUE );
    {
        vu::store( vu::set1(vu::BIAS), mat(0) );
    }
//     std:: ofstream os( "/tmp/yyy.txt" );
    vec_t last_s;
    vec_t diag_init = vu::set1(vu::BIAS);
    
    for( size_t ia = 0; ia < asize; ia++ ) {
        
        //const value_t a = *ait;
        const vec_t a = vu::load(aprofile(ia * W));
        
        last_s = vu::set1(vu::LARGE_VALUE);

        // for the first row, the diagonal value for the first colum is zero (inialization above)
        vec_t diag = diag_init;
        // for the remaining rows it will be infinity. note: mat[0] is never written to in the b-loop
        if( ia == 0 ) {
            diag_init = vu::set1(vu::LARGE_VALUE);
        }
        
        score_t *upper = mat(W);
        
        for( vec_iter_t bit = b_begin; bit != b_end; ++bit, upper += W ) {
            const vec_t b = vu::set1(*bit);
            
            const vec_t sd = diag;
            const vec_t su = vu::load(upper);
            const vec_t sl = last_s;
            
            
            
            const vec_t cost = vu::abs_diff(a, b);
//             const vec_t cost = vu::max( vu::sub(a, b), vu::sub(b, a) );
            //std::cout << "cost: " << cost << "\n";
            
            diag = su;
            last_s = vu::adds(cost, vu::min( sd, vu::min(su, sl) ) );
            
            vu::store( last_s, upper );
            

        }
        
    }
    
    
    if( ps.out.size() != W ) {
        ps.out.resize(W);
    }
    vu::store( last_s, ps.out.base() );
    
    out.assign( ps.out.begin(), ps.out.end() );
    //return mat[mat.size() - 1];
   
}
 


#endif
