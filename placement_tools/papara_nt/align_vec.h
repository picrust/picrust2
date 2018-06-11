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

#ifndef __align_vec_h
#define __align_vec_h

#include "vec_unit.h"
#include "ivymike/aligned_buffer.h"
#include "ivymike/fasta.h"

using ivy_mike::aligned_buffer;
using ivy_mike::scoring_matrix;

template <class score_t>
struct persistent_state {
    aligned_buffer<score_t> out;
    aligned_buffer<score_t> s;
    aligned_buffer<score_t> si;  
};

//#define LIKELY(x) __builtin_expect((x),1)
// #define LIKELY(x) (x)
template <class score_t, class sscore_t, size_t W>
void align_vec( persistent_state<score_t> &ps, size_t asize, const std::vector<uint8_t> &b, const scoring_matrix &m, aligned_buffer<sscore_t> &qprofile, const sscore_t gap_open, const sscore_t gap_extend, std::vector<int> &out ) {
 
    typedef vector_unit<score_t,W> vu;
 
    typedef typename vu::vec_t vec_t;
    
//     assert( a.size() % W == 0 );
//     const size_t asize = a.size() / W;
    
//     std::vector<score_t> s( a.size());
//     std::vector<score_t> si( a.size());
    
    aligned_buffer<score_t> &s = ps.s;
    aligned_buffer<score_t> &si = ps.si;
    
    aligned_buffer<score_t> pb(W);
    
    if( s.size() < asize * W ) {
        s.resize( asize * W );
        si.resize( asize * W );
    }
    const score_t bias = vu::BIAS;
    std::fill( s.begin(), s.end(), bias );
    std::fill( si.begin(), si.end(), bias );
    const score_t SMALL = vu::SMALL_VALUE;
    vec_t max_vec = vu::set1(SMALL);
    
    const vec_t GAP_EXT_vec = vu::set1(score_t(-gap_extend)); // these values are _subtracted_ from the score. so use negative of user parameters!
    const vec_t GAP_OPEN_vec = vu::set1(score_t(-gap_open));
    const vec_t BIAS_vec = vu::set1( vu::BIAS );
    
//     vec_t len_vec = vu::load( len.m_ptr );
#define LOCAL_ALIGN
    for( size_t ib = 0; ib < b.size(); ib++ ) {
#ifndef LOCAL_ALIGN        
        const bool lastrow = ib == (b.size() - 1);
#endif        
        char bc = b[ib];
//         std::cout << "bc: " << bc << std::endl;
        assert( bc < char(qprofile.size() / (W * asize)));
        
        vec_t last_sl = vu::set1(SMALL);
        vec_t last_s = vu::set1(vu::BIAS);
        vec_t last_sdiag = vu::set1(vu::BIAS);
//         std::cout << "sbm: " << m.state_backmap(bc) << " " << (W * asize) << std::endl;
       // sscore_t *qpp_iter = qprofile(m.state_backmap(bc) * W * asize);
        sscore_t *qpp_iter = qprofile( bc * W * asize);
        score_t * __restrict s_iter = s.base();
        score_t * __restrict su_iter = si.base();
        score_t * __restrict s_end = s_iter + (asize * W);
        
        for( ; s_iter != s_end; s_iter += W, su_iter += W, qpp_iter += W ) {
            
            const vec_t match = vu::load( (score_t*) qpp_iter );
            
//             getchar();
            const vec_t sm = vu::add( last_sdiag, match );
            
            last_sdiag = vu::load( s_iter );
            
#ifdef LOCAL_ALIGN            
            const vec_t sm_zero = vu::max( sm, BIAS_vec);
#else
            const vec_t sm_zero = sm_vec;
#endif
            
            const vec_t sl_open = vu::sub( last_s, GAP_OPEN_vec ); 
            const vec_t sl_ext = vu::sub( last_sl, GAP_EXT_vec );
            const vec_t sl = vu::max( sl_ext, sl_open );
            
            
            last_sl = sl;
            

            const vec_t su_ext = vu::sub( vu::load(su_iter), GAP_EXT_vec );
            const vec_t su_open = vu::sub( last_sdiag, GAP_OPEN_vec );
            const vec_t su = vu::max( su_ext, su_open );
            
            vu::store( su, su_iter );
            
            const vec_t s = vu::max( sm_zero, vu::max( sl, su ));
            
//             vu::println( max_vec, pb.m_ptr );
            
            last_s = s;
            vu::store( s, s_iter ); 
#ifdef LOCAL_ALIGN  
            max_vec = vu::max( s, max_vec );
            
#else
            if( s_iter == s_end - W || lastrow ) {
                max_vec = vu::max( sc_vec, max_vec );
            }
#endif
        }
        
    }
    
    ps.out.resize(W);
    vu::store( max_vec, ps.out.base() );
    
    out.resize(0);
    out.reserve(W);
    for( size_t i = 0; i < W; i++ ) {
        out.push_back(int(ps.out[i]) - vu::BIAS );
//         std::cout << "x: " << int(ps.out.m_ptr[i]) << "\n";
    }
    //return max;
}


template<typename score_t, typename sscore_t>
struct block_cont {
    score_t * __restrict s_iter;
    score_t * __restrict su_iter;
    score_t * __restrict s_end;
    sscore_t * __restrict qpp_iter;
};

template <class score_t, class sscore_t>
struct persistent_state_blocked {
    aligned_buffer<score_t> out;
    aligned_buffer<score_t> s;
    aligned_buffer<score_t> si;  
    aligned_buffer<score_t> bcv_store;
    std::vector<block_cont<score_t, sscore_t> > bc_array;
};

template <class score_t, class sscore_t, size_t W, size_t BW>
void align_vec_blocked( persistent_state_blocked<score_t,sscore_t> &ps, size_t asize, const std::vector<uint8_t> &b, const scoring_matrix &m, aligned_buffer<sscore_t> &qprofile, const sscore_t gap_open, const sscore_t gap_extend, std::vector<int> &out ) {
 
    typedef vector_unit<score_t,W> vu;
 
    typedef typename vu::vec_t vec_t;
    
//     assert( a.size() % W == 0 );
//     const size_t asize = a.size() / W;
    
//     std::vector<score_t> s( a.size());
//     std::vector<score_t> si( a.size());
    
    aligned_buffer<score_t> &s = ps.s;
    aligned_buffer<score_t> &si = ps.si;
    
    
//     aligned_buffer<score_t> pb(W);
    
    if( s.size() < asize * W ) {
        s.resize( asize * W );
        si.resize( asize * W );
    }
    const score_t bias = vu::BIAS;
    std::fill( s.begin(), s.end(), bias );
    std::fill( si.begin(), si.end(), bias );
    const score_t SMALL = vu::SMALL_VALUE;
    vec_t max_vec = vu::set1(SMALL);
    
    const vec_t GAP_EXT_vec = vu::set1(score_t(-gap_extend)); // these values are _subtracted_ from the score. so use negative of user parameters!
    const vec_t GAP_OPEN_vec = vu::set1(score_t(-gap_open));
    const vec_t BIAS_vec = vu::set1( vu::BIAS );
    
    
    
    
    
    struct block_cont_vec {
        vec_t last_sl;
        vec_t last_s;
        vec_t last_sdiag;
        
        void init() {
            last_sl = vu::set1(vu::SMALL_VALUE);
            last_s = vu::set1(vu::BIAS);
            last_sdiag = vu::set1(vu::BIAS);
        }
        
        void load( score_t *ptr ) {
            last_sl = vu::load( ptr );
            last_s = vu::load( ptr + W ); 
            last_sdiag = vu::load( ptr + W + W );
            
        }
        
        void store( score_t *ptr ) {
            vu::store( last_sl, ptr );
            vu::store( last_s, ptr + W );
            vu::store( last_sdiag, ptr + W + W );
        }
        
        
    };
    
    size_t num_blocks = asize / BW;
    if( asize % BW != 0 ) {
        num_blocks++;
    }
    
    //aligned_buffer<score_t> bcv_store( b.size() * W * 3 );
    aligned_buffer<score_t> &bcv_store = ps.bcv_store;
    std::vector<block_cont<score_t, sscore_t> > &bc_array = ps.bc_array;
    
    if( bcv_store.size() < b.size() * W * 3 ) {
        bcv_store.resize(b.size() * W * 3);
        bc_array.resize( b.size() );
    }
    
    //std::vector<block_cont> bc_array(b.size());
    
    
    
//     vec_t len_vec = vu::load( len.m_ptr );
    
    for( size_t iblock = 0; iblock < num_blocks; ++iblock ) {
        block_cont_vec bcv_;
        
    
#define LOCAL_ALIGN
        for( size_t ib = 0; ib < b.size(); ib++ ) {
#ifndef LOCAL_ALIGN        
            const bool lastrow = ib == (b.size() - 1);
#endif        
            block_cont<score_t, sscore_t> &bc_ = bc_array[ib];
            
            
            if( iblock == 0 ) {
                const char bchar = b[ib];
                assert( bchar < char(qprofile.size() / (W * asize)));
                bcv_.init();
                bc_.s_iter = s.base();
                bc_.su_iter = si.base();
                bc_.s_end = bc_.s_iter + (asize * W);
                bc_.qpp_iter = qprofile( bchar * W * asize );
            } else {
                bcv_.load( bcv_store( ib * W * 3 ) );   
            }

            
            
            //         std::cout << "bc: " << bc << std::endl;
            
            
//             vec_t last_sl = vu::set1(SMALL);
//             vec_t last_s = vu::set1(vu::BIAS);
//             vec_t last_sdiag = vu::set1(vu::BIAS);
            //         std::cout << "sbm: " << m.state_backmap(bc) << " " << (W * asize) << std::endl;
            sscore_t *qpp_iter = bc_.qpp_iter;
            
            _mm_prefetch( (const char *) qpp_iter,_MM_HINT_T0);

            score_t * __restrict s_iter = bc_.s_iter;
            score_t * __restrict su_iter = bc_.su_iter;
            score_t * __restrict s_end = bc_.s_end;
            
            vec_t last_sdiag = bcv_.last_sdiag;
            vec_t last_sl = bcv_.last_sl;
            vec_t last_s = bcv_.last_s;
            
            score_t * __restrict s_block_end;// = std::min( bc.s_iter + BW * W, bc.s_end );
            if( s_iter + BW * W < s_end ) {
                s_block_end = s_iter + BW * W;
            } else {
                s_block_end = s_end;
            }
            
            
            assert( s_block_end != s_iter );
            
            for( ; s_iter != s_block_end; s_iter += W, su_iter += W, qpp_iter += W ) {
                
                const vec_t match = vu::load( (score_t*) qpp_iter );
                
                //             getchar();
                const vec_t sm = vu::add( last_sdiag, match );
                
                last_sdiag = vu::load( s_iter );
                
#ifdef LOCAL_ALIGN            
                const vec_t sm_zero = vu::max( sm, BIAS_vec);
#else
                const vec_t sm_zero = sm_vec;
#endif
                
                const vec_t sl_open = vu::sub( last_s, GAP_OPEN_vec ); 
                const vec_t sl_ext = vu::sub( last_sl, GAP_EXT_vec );
                const vec_t sl = vu::max( sl_ext, sl_open );
                
                
                last_sl = sl;
                
                
                const vec_t su_ext = vu::sub( vu::load(su_iter), GAP_EXT_vec );
                const vec_t su_open = vu::sub( last_sdiag, GAP_OPEN_vec );
                const vec_t su = vu::max( su_ext, su_open );
                
                vu::store( su, su_iter );
                
                const vec_t s = vu::max( sm_zero, vu::max( sl, su ));
                
                //             vu::println( max_vec, pb.m_ptr );
                
                last_s = s;
                vu::store( s, s_iter ); 
#ifdef LOCAL_ALIGN  
                max_vec = vu::max( s, max_vec );
                
#else
                if( s_iter == s_end - W || lastrow ) {
                    max_vec = vu::max( sc_vec, max_vec );
                }
#endif
            }
            
            bcv_.last_sdiag = last_sdiag;
            bcv_.last_sl = last_sl;
            bcv_.last_s = last_s;
            
            bcv_.store( bcv_store( ib * W * 3 ) );   
            bc_.qpp_iter = qpp_iter;
            bc_.s_iter = s_iter;
            bc_.su_iter = su_iter;
            
        }
    }
    
    
    
    ps.out.resize(W);
    vu::store( max_vec, ps.out.base() );
    
    out.resize(0);
    out.reserve(W);
    for( size_t i = 0; i < W; i++ ) {
        out.push_back(int(ps.out[i]) - vu::BIAS );
//         std::cout << "x: " << int(ps.out.m_ptr[i]) << "\n";
    }
    //return max;
}


#endif
