/*
 * Copyright (C) 2011 Simon A. Berger
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 */

#include <iostream>

#include "pars_align_vec.h"
#include "vec_unit.h"


// typedef vector_unit<short,pars_align_vec::WIDTH> vu_t;
typedef vector_unit<pars_align_vec::score_t,pars_align_vec::WIDTH> vu_t;

void pars_align_vec::align_freeshift_s11() {

    
    vu_t::vec_t z = vu_t::setzero(); 
   
    for( int ia = -1; ia <= int(m_na - m_nb - 1); ia++ ) {
        size_t p = svaddr(ia, -1);
        vu_t::store( z, m_arrays.m_score(p) );
        vu_t::store( z, m_arrays.m_score_l(p) );
    }
    
    for( size_t i = 0; i < m_na; i++ ) {
     
        vu_t::store( z, m_arrays.m_score( svaddr(i, -1) ));
    }
    
    
    vu_t::vec_t l = vu_t::set1(LARGE_VALUE); 
   
    for( size_t ia = int(m_na - m_nb); ia < m_na; ia++ ) {
        size_t p = svaddr(ia, -1);
        vu_t::store( l, m_arrays.m_score(p) );
    }
    
    for( size_t i = 0; i < m_nb; i++ ) {
    
        size_t p_1_0 = svaddr(i-1, i);
        size_t p_0_0 = svaddr(i, i);
        size_t p_x_1 = svaddr((m_na - m_nb) + i, i - 1 );
        
        vu_t::store( l, m_arrays.m_score(p_1_0) );
        vu_t::store( l, m_arrays.m_score_l(p_1_0) );
        vu_t::store( l, m_arrays.m_score(p_0_0) );
        vu_t::store( l, m_arrays.m_score(p_x_1) );
    }
    
    
}
void pars_align_vec::align_freeshift_s3() {
//     score_t best = LARGE_VALUE;
    const size_t band_width = m_na - m_nb;
    
    #define AUX_CGAP ( 0x1 )
    
    
    
    
    for ( size_t ib = 0; ib < m_nb; ib++ ) {
        pars_state_t cb = m_seq_b[ib];
        
        vu_t::vec_t cb_vec = vu_t::set1(cb);
        
        int astart = ib;
//         size_t p_as_b = svaddr(astart-1, ib);
        
//         vu_t::vec_t last_sp = vu_t::load( m_arrays.m_score(p_as_b));
//         vu_t::vec_t last_sLp = vu_t::load( m_arrays.m_score_l(p_as_b ));
        vu_t::vec_t last_sp = vu_t::set1(LARGE_VALUE);
        vu_t::vec_t last_sLp = vu_t::set1(LARGE_VALUE);
        
        int ibl2 = -1;
        int ibl = 0;
        vu_t::vec_t last_diag = vu_t::load(m_arrays.m_score(svaddr(astart-1,ibl-1)));         
//         std::cout << _mm_extract_epi16( last_diag, 0 ) << "\n";
        for ( size_t ia = astart; ia <= ib + band_width; ia++ /*, saddr_0_0++, saddr_1_1++, sp_0_0++, sp_1_1++*/ ) {
//             const pars_state_t ca = get_seq_a( ia );
//             vu_t::vec_t ca_vec = vu_t::set1(ca);
            
            vu_t::vec_t ca_vec = vu_t::load( m_seq_a(WIDTH * ia) );
            
            
            // CGAP mask: aux_cgap_mask_vec[i] == 0xffff if cgap column
            
#if 0
            vu_t::vec_t aux_vec = vu_t::load( m_aux_a.m_ptr + (WIDTH * ia ));
            
            vu_t::vec_t aux_cgap_mask_vec = vu_t::cmp_eq( aux_vec, vu_t::set1(AUX_CGAP));
            
            vu_t::vec_t aux_match_cgap_vec = vu_t::bit_and( aux_cgap_mask_vec, vu_t::set1( MATCH_CGAP ));
#else
/*            vu_t::vec_t aux_cgap_mask_vec = vu_t::load( m_cgap_mask_pre(WIDTH * ia));
            vu_t::vec_t aux_match_cgap_vec = vu_t::bit_and( aux_cgap_mask_vec, vu_t::set1( MATCH_CGAP ));*/
            const vu_t::vec_t aux_match_cgap_vec = vu_t::load( m_match_cgap_pre(WIDTH * ia));
#endif
            const vu_t::vec_t mismatch_vec = vu_t::set1(MISMATCH_PENALTY);
            // match comparison: bitwise and between parsimony states
            // RELUSLT: sd_add[i] == MISMATCH_PENALTY iff match[i] == 0
            
            const vu_t::vec_t match = vu_t::bit_and( ca_vec, cb_vec );
            const vu_t::vec_t sd_add = vu_t::bit_and(vu_t::cmp_zero( match ), mismatch_vec );
            
            
            // load upper and diagonal cells.
            // add MISMATCH_PENALTY to upper value
            
            
            
//             vu_t::vec_t su_vec = vu_t::add( mismatch_vec, vu_t::load( m_arrays.m_score(svaddr(ia, -1)) ));
            
//             vu_t::vec_t sd_vec = vu_t::load( m_arrays.m_score(svaddr(ia-1, ib-1)) );
            vu_t::vec_t sd_vec = last_diag;
            
            last_diag = vu_t::load( m_arrays.m_score(svaddr(ia, ibl-1)));
            const vu_t::vec_t su_vec = vu_t::add( mismatch_vec, last_diag );
            
            //vu_t::vec_t sd_vec = vu_t::load( m_arrays.m_score(svaddr(ia-1, -1)) );
//             vu_t::vec_t sd_vec = last_diag;
//             last_diag = su_vec;
            sd_vec = vu_t::add( sd_vec, vu_t::add( sd_add, aux_match_cgap_vec ) );
           
            // calculate gap open/extend scores
            // influnece of cgap: only add penalties if not CGAP:
            // bitwise andnot with the cgap_mask vector does the trick...
#if 0
            const vu_t::vec_t gap_extend_vec = vu_t::bit_andnot( aux_cgap_mask_vec, vu_t::set1( GAP_EXTEND ));
            const vu_t::vec_t gap_open_ext_vec = vu_t::bit_andnot( aux_cgap_mask_vec, vu_t::set1( GAP_OPEN_EXTEND ));
#else
            const vu_t::vec_t gap_extend_vec = vu_t::load( m_extend_pre( ia * WIDTH ));;
            const vu_t::vec_t gap_open_ext_vec = vu_t::load( m_open_ext_pre( ia * WIDTH ));
#endif
            const vu_t::vec_t score_open = vu_t::add( last_sp, gap_open_ext_vec );
            const vu_t::vec_t score_extend = vu_t::add( last_sLp, gap_extend_vec );
            
            // decide on QS gap open vs. extension
            
            last_sLp = vu_t::min( score_open, score_extend );
            
            // decide on match vs. RS gap
            sd_vec = vu_t::min( sd_vec, su_vec );
            
          
            last_sp = vu_t::min( last_sLp, sd_vec );
            const size_t addr = svaddr( ia, ibl2 );
            
            vu_t::store( last_sp, m_arrays.m_score(addr) );
//             vu_t::store( last_sp, m_arrays.m_score(svaddr( ia, -1 )));
//             vu_t::store( vu_t::set1(LARGE_VALUE), m_arrays.m_score(svaddr( ia + 1, -1 )));
//             if( ia == astart ) {
//              
//                 last_diag = last_sp;
//             }
        }
        
        
        
    }
    
//     throw std::runtime_error("exit");
    // look for minimum values in the last row
    const int tb_start_b = -1;
    //const int tb_start_b = -1;
    vu_t::vec_t best_score = vu_t::set1( LARGE_VALUE );
//         std::cout << "find best\n";
    for ( size_t a = m_na - 1; a >= m_nb - 1; a-- ) {
        vu_t::vec_t cand = vu_t::load(m_arrays.m_score(svaddr ( a, tb_start_b )));
        
//         std::cout << " " << a << ":" << _mm_extract_epi16(cand, 0) << "\n";
        best_score = vu_t::min(best_score, cand);
    }
    
    
//     throw std::runtime_error("exit");
    vu_t::store(best_score, m_best.base());
    
    m_ncups = band_width * m_nb * WIDTH;
//     std::cout << "best: " << m_best.m_ptr[0] << "\n";
}

pars_align_vec::pars_align_vec(const int** seqA, unsigned char* seqB, size_t n_a, size_t n_b, size_t aStride, const unsigned int** aAux,
                               size_t aAuxStride, pars_align_vec::arrays<WIDTH>& arr, const unsigned int* bvtrans,
                               pars_align_vec::score_t gapOpen, pars_align_vec::score_t gap_extend,
                               pars_align_vec::score_t mismatch, pars_align_vec::score_t match_cgap) 

    : 
    GAP_OPEN(gapOpen),
    GAP_EXTEND(gap_extend),
    GAP_OPEN_EXTEND( GAP_OPEN + GAP_EXTEND ),
    MISMATCH_PENALTY(mismatch),
    MATCH_CGAP(match_cgap),
    m_na( n_a ), m_nb( n_b ),
    m_seqlist_a( seqA ), m_auxlist_a(aAux), m_seq_b(seqB),
    m_stride_a( aStride ), m_stride_aux_a(aAuxStride),
    m_best(WIDTH),
    m_seq_a(WIDTH * n_a ),
    m_aux_a(WIDTH * n_a ),
//     m_cgap_mask_pre( WIDTH * n_a ),
    m_match_cgap_pre(WIDTH * n_a),
    m_extend_pre(WIDTH * n_a),
    m_open_ext_pre(WIDTH * n_a),
    
    m_arrays(arr),
    m_ncups(0)
{    
    m_arrays.size(ma() * mb());
    
    score_t *sap = m_seq_a.base();
    score_t *aap = m_aux_a.base();
    score_t *mcp = m_match_cgap_pre.base();
    score_t *ep = m_extend_pre.base();
    score_t *oep = m_open_ext_pre.base();
//     score_t *cmp = m_cgap_mask_pre.m_ptr;
    
    for( size_t i = 0; i < m_na; i++ ) {
        score_t *aap_start = aap;
        
        for( size_t j = 0; j < WIDTH; j++ ) {
        
            *sap = get_seq_a(j, i);
            *aap = get_aux_a(j, i);
            
            sap++;
            aap++;
        }
        vu_t::vec_t aux_vec = vu_t::load( aap_start);
            
        vu_t::vec_t aux_cgap_mask_vec = vu_t::cmp_eq( aux_vec, vu_t::set1(AUX_CGAP));
//         vu_t::store( aux_cgap_mask_vec, m_cgap_mask_pre(i * WIDTH));
        
        vu_t::vec_t aux_match_cgap_vec = vu_t::bit_and( aux_cgap_mask_vec, vu_t::set1( MATCH_CGAP ));
        
        vu_t::store( aux_match_cgap_vec, mcp );
        mcp += WIDTH;
        
        const vu_t::vec_t gap_extend_vec = vu_t::bit_andnot( aux_cgap_mask_vec, vu_t::set1( GAP_EXTEND ));
        vu_t::store( gap_extend_vec, ep );
        ep += WIDTH;
        
        const vu_t::vec_t gap_open_ext_vec = vu_t::bit_andnot( aux_cgap_mask_vec, vu_t::set1( GAP_OPEN_EXTEND ));
        vu_t::store( gap_open_ext_vec, oep );
        oep += WIDTH;
        
    }
    
}
