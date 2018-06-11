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


#ifndef __align_pvec_vec_h
#define __align_pvec_vec_h

#include <iostream>
#include <vector>
#include <algorithm>

#include "ivymike/aligned_buffer.h"
#include "vec_unit.h"

template<typename score_t, size_t W>
class align_pvec_score{
    ivy_mike::aligned_buffer<score_t> m_a_prof;
    ivy_mike::aligned_buffer<score_t> m_cgap_prof;
    
    mutable ivy_mike::aligned_buffer<score_t> m_s;
    mutable ivy_mike::aligned_buffer<score_t> m_out_score;
    
    
    const score_t m_mismatch_score;
    const score_t m_match_cgap;
    const score_t m_gap_open;
    const score_t m_gap_extend;
    
    

    //aligned_buffer<score_t> m_a_aux_prof;
    const static score_t aux_cgap = 0x1;
//     score_t map_match_cgap( score_t aux ) {
//         if( aux == aux_cgap ) {
//             return m_match_cgap;
//         } else {
//             return 0;
//         }
//     }
//     
//     score_t map_gap_extend( score_t aux ) {
//         if( aux == aux_cgap ) {
//             return 0;
//         } else {
//             return m_gap_extend;
//         }
//     }
//     
//     score_t map_gap_oext( score_t aux ) {
//         if( aux == aux_cgap ) {
//             return 0;
//         } else {
//             return m_gap_open + m_gap_extend;
//         }
//     }
    
    score_t map_cgap( score_t aux ) {
        if( aux == aux_cgap ) {
            return -1; // FIXME: will this always yield a score_t with all bits set?
        } else {
            return 0;
        }
    }
public:
    align_pvec_score( const ivy_mike::aligned_buffer<score_t> &a_prof, const ivy_mike::aligned_buffer<score_t> &a_aux_prof, score_t mismatch_score, score_t match_cgap, score_t gap_open, score_t gap_extend )
     : m_mismatch_score(mismatch_score),m_match_cgap(match_cgap), m_gap_open(gap_open), m_gap_extend( gap_extend )
    {
        m_a_prof = a_prof;
       
        
        m_cgap_prof.reserve( a_aux_prof.size());
        for( typename ivy_mike::aligned_buffer<score_t>::const_iterator it = a_aux_prof.begin(); it != a_aux_prof.end(); ++it ) {
            
            m_cgap_prof.push_back( map_cgap(*it));
        }
        
        m_s.resize( a_prof.size() );
        m_out_score.resize(W);
    }
    
    
    align_pvec_score( const int **seqptrs, const unsigned int **auxptrs, size_t reflen, score_t mismatch_score, score_t match_cgap, score_t gap_open, score_t gap_extend )
     : m_mismatch_score(mismatch_score),m_match_cgap(match_cgap), m_gap_open(gap_open), m_gap_extend( gap_extend )
    {
        m_a_prof.reserve( reflen * W );
        m_cgap_prof.reserve( reflen * W );
        for( size_t i = 0; i < reflen; ++i ) {
            for( size_t j = 0; j < W; ++j ) {
                m_a_prof.push_back( seqptrs[j][i] );
                m_cgap_prof.push_back( map_cgap(auxptrs[j][i]));    
            }
            
            
        }
        
        assert( m_a_prof.size() == reflen * W );
        
        m_s.resize( m_a_prof.size() );
        m_out_score.resize(W);
    }
    

    inline void align( const std::vector<uint8_t> &b ) const {
        align(b.begin(), b.end() );
    }

    template<typename iiter>
    inline void align( iiter bstart, iiter bend ) const {
        typedef vector_unit<score_t, W> vu;
        typedef typename vu::vec_t vec_t;
        
        
        size_t asize = m_a_prof.size() / W;
        assert( ptrdiff_t(asize) > std::distance(bstart, bend) );
        
        const size_t band_width = asize - std::distance(bstart, bend);
        std::fill( m_s.begin(), m_s.end(), 0 );
        
        const score_t LARGE = 32000;
        //     std::fill( arr.si.begin(), arr.si.end(), LARGE );
        std::fill( m_s.begin() + (band_width + 1) * W, m_s.begin() + (band_width + 2) * W, LARGE );
        
        
        const vec_t gap_oext = vu::set1(m_gap_open + m_gap_extend);
        const vec_t gap_ext = vu::set1(m_gap_extend);
        const vec_t match_cgap = vu::set1( m_match_cgap );
        const vec_t mismatch_score = vu::set1(m_mismatch_score);
        const vec_t zero = vu::setzero();
        for( size_t ib = 0; bstart != bend; ++bstart, ++ib ) {
            const vec_t bc = vu::set1(*bstart);
            
            vec_t last_sl = vu::set1(LARGE);
            vec_t last_sc = vu::set1(LARGE);
            vec_t last_sdiag = vu::set1(0);
            
            size_t astart = ib;
            
            score_t * __restrict s_iter = &m_s[0];
            score_t * __restrict s_iter_next = &m_s[W];
            
            last_sdiag = vu::load(s_iter);
            
            
            for( size_t ia = astart; ia <= ib + band_width; ++ia, s_iter += W, s_iter_next += W ) {  
                const vec_t ac = vu::load( m_a_prof(ia * W) );
                const vec_t cgap = vu::load( m_cgap_prof(ia * W)); 
               
                const vec_t sm_match = vu::bit_and( mismatch_score, vu::cmp_eq( vu::bit_and( ac, bc ), zero) );
                const vec_t sm_cgap = vu::bit_and( match_cgap, cgap );
                const vec_t sm = vu::add( vu::add( last_sdiag, sm_match ), sm_cgap );
                
                
                
                const vec_t sl_ext = vu::add( last_sl, vu::bit_andnot( cgap, gap_ext));
                const vec_t sl_open = vu::add( last_sc, vu::bit_andnot( cgap, gap_oext));
                last_sl = vu::min( sl_ext, sl_open );
                
                const vec_t min_sm_sl = vu::min( sm, last_sl );
                last_sdiag = vu::load(s_iter_next);
                const vec_t su = vu::add( last_sdiag, mismatch_score );
                const vec_t sc = vu::min( min_sm_sl, su );
                
                last_sc = sc;
                
                vu::store( sc, s_iter );
                
                
            }
        }
        
        vec_t minscore = vu::set1(LARGE);
        for( size_t i = 0; i < band_width + 1; ++i ) {
            minscore = vu::min( minscore, vu::load( m_s(i * W)));
        }
        
        vu::store( minscore, m_out_score(0));
    }
    
    
    const score_t *get_scores() const {
        return m_out_score.data();
    }
    
};



template<size_t W>
class align_pvec_gapp_score{

    typedef int32_t pstate_t;
    typedef float score_t;


    ivy_mike::aligned_buffer<pstate_t> m_a_prof;
    ivy_mike::aligned_buffer<score_t> m_gapp_prof;

    mutable ivy_mike::aligned_buffer<score_t> m_s;
    mutable ivy_mike::aligned_buffer<score_t> m_out_score;


    const score_t m_mismatch_score;
    const score_t m_match_cgap;
    const score_t m_gap_open;
    const score_t m_gap_extend;




//    score_t map_cgap( score_t aux ) {
//        if( aux == aux_cgap ) {
//            return -1; // FIXME: will this always yield a score_t with all bits set?
//        } else {
//            return 0;
//        }
//    }
public:
    align_pvec_gapp_score( const ivy_mike::aligned_buffer<pstate_t> &a_prof, const ivy_mike::aligned_buffer<score_t> &a_gapp_prof, score_t mismatch_score, score_t match_cgap, score_t gap_open, score_t gap_extend )
     : m_mismatch_score(mismatch_score),m_match_cgap(match_cgap), m_gap_open(gap_open), m_gap_extend( gap_extend )
    {
        m_a_prof = a_prof;
        m_gapp_prof = a_gapp_prof;

        m_s.resize( a_prof.size() );
        m_out_score.resize(W);
    }


    align_pvec_gapp_score( const int **seqptrs, const double **gapp_ptrs, size_t reflen, score_t mismatch_score, score_t match_cgap, score_t gap_open, score_t gap_extend )
     : m_mismatch_score(mismatch_score),m_match_cgap(match_cgap), m_gap_open(gap_open), m_gap_extend( gap_extend )
    {
        m_a_prof.reserve( reflen * W );
        m_gapp_prof.reserve( reflen * W );
        for( size_t i = 0; i < reflen; ++i ) {
            for( size_t j = 0; j < W; ++j ) {
                m_a_prof.push_back( seqptrs[j][i] );
                m_gapp_prof.push_back( gapp_ptrs[j][i] );
            }


        }

        assert( m_a_prof.size() == reflen * W );

        m_s.resize( m_a_prof.size() );
        m_out_score.resize(W);
    }

    inline void align( const std::vector<uint8_t> &b ) const {
        typedef vector_unit<score_t, W> vu;
        typedef vector_unit<pstate_t, W> vui;

        typedef typename vu::vec_t vec_t;
        typedef typename vui::vec_t psvec_t;

        size_t asize = m_a_prof.size() / W;
        assert( asize > b.size() );

        const size_t band_width = asize - b.size();
        std::fill( m_s.begin(), m_s.end(), 0 );

        const score_t LARGE = 32000;
        //     std::fill( arr.si.begin(), arr.si.end(), LARGE );
        std::fill( m_s.begin() + (band_width + 1) * W, m_s.begin() + (band_width + 2) * W, LARGE );


        const vec_t gap_oext = vu::set1(m_gap_open + m_gap_extend);
        const vec_t gap_ext = vu::set1(m_gap_extend);
        const vec_t match_cgap = vu::set1( m_match_cgap );
        const vec_t mismatch_score = vu::set1(m_mismatch_score);
        const psvec_t zero = vui::setzero();
        for( size_t ib = 0; ib < b.size(); ib++ ) {
            const psvec_t bc = vui::set1(b[ib]);

            vec_t last_sl = vu::set1(LARGE);
            vec_t last_sc = vu::set1(LARGE);
            vec_t last_sdiag = vu::set1(0);

            size_t astart = ib;

            score_t * __restrict s_iter = &m_s[0];
            score_t * __restrict s_iter_next = &m_s[W];

            last_sdiag = vu::load(s_iter);


            for( size_t ia = astart; ia <= ib + band_width; ++ia, s_iter += W, s_iter_next += W ) {
                const psvec_t ac = vui::load( m_a_prof(ia * W) );
//                const vec_t cgap = vu::load( m_cgap_prof(ia * W));

                const psvec_t match_mask_int = vui::cmp_eq( vui::bit_and( ac, bc ), zero);
                const vec_t match_mask = vu::cast_from_int(match_mask_int);
                const vec_t sm_match = vu::bit_and( mismatch_score, match_mask );


                const vec_t p_nongap = vu::load( m_gapp_prof( ia * W ) );
                const vec_t p_gap = vu::sub( vu::set1( 1.0 ), p_nongap );

                const vec_t sm_cgap = vu::mul( match_cgap, p_gap );
                const vec_t sm = vu::add( vu::add( last_sdiag, sm_match ), sm_cgap );



                const vec_t sl_ext = vu::add( last_sl, vu::mul( p_nongap, gap_ext));
                const vec_t sl_open = vu::add( last_sc, vu::mul( p_nongap, gap_oext));
                last_sl = vu::min( sl_ext, sl_open );

                const vec_t min_sm_sl = vu::min( sm, last_sl );
                last_sdiag = vu::load(s_iter_next);
                const vec_t su = vu::add( last_sdiag, gap_oext );
                const vec_t sc = vu::min( min_sm_sl, su );

                last_sc = sc;

                vu::store( sc, s_iter );


            }
        }

        vec_t minscore = vu::set1(LARGE);
        for( size_t i = 0; i < band_width + 1; ++i ) {
            minscore = vu::min( minscore, vu::load( m_s(i * W)));
        }

        vu::store( minscore, m_out_score(0));
    }


    const score_t *get_scores() const {
        return m_out_score.data();
    }

};


#endif
