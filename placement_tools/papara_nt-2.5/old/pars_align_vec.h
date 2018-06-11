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

#ifndef __pars_align_vec_h
#define __pars_align_vec_h

#include <cstdlib>
#include <cstddef>
#include <stdexcept>
#include "aligned_buffer.h"

class pars_align_vec {
    
    
public:    
    //typedef short score_t;
#if 0
#error BOGUS!
    typedef unsigned char score_t;
    const static size_t WIDTH = 16;
    const static score_t LARGE_VALUE = 110; //  score_t must be able to keep LARGE_VALUE + GAP_OPEN + GAP_EXTEND without overflowing!. waiting for c++0x to be able to use std::numeric_limits at compile-time...
    
    
#else
    typedef short score_t;
    const static size_t WIDTH = 8;
    const static score_t LARGE_VALUE = 32000; //  score_t must be able to keep LARGE_VALUE + GAP_OPEN + GAP_EXTEND without overflowing!. waiting for c++0x to be able to use std::numeric_limits at compile-time...
    
#endif
    
    
    
    typedef unsigned int pars_state_t;
    
    template<size_t N>
    struct arrays {
    
        aligned_buffer<score_t> m_score;
        aligned_buffer<score_t> m_score_l;

        void size( const size_t s ) {
            const size_t vs = N * s;
            
            if( vs > m_score.size() ) {
                m_score.resize(vs);
                m_score_l.resize(vs);
            }
        }
        
    };
    
    
private:
    const score_t GAP_OPEN; //= 1;
    const score_t GAP_EXTEND; // = 1;
    const score_t GAP_OPEN_EXTEND;
    const score_t MISMATCH_PENALTY;// = 4;
    const score_t MATCH_CGAP;
    const size_t m_na;
    const size_t m_nb;
    
    const int **m_seqlist_a;
    const unsigned int **m_auxlist_a;
    const unsigned char *m_seq_b;
    
    
    
    
    const size_t m_stride_a;
    const size_t m_stride_aux_a;
    
    aligned_buffer<score_t> m_best;
    
    aligned_buffer<score_t> m_seq_a;
    aligned_buffer<score_t> m_aux_a;
    
//     aligned_buffer<score_t> m_cgap_mask_pre;
    aligned_buffer<score_t> m_match_cgap_pre;
    aligned_buffer<score_t> m_extend_pre;
    aligned_buffer<score_t> m_open_ext_pre;
    
    
    
    arrays<WIDTH> &m_arrays;
    
    size_t m_ncups;
    
    
    
    inline size_t ma() { return m_na + 1; }
    inline size_t mb() { return m_nb + 1; }
    
    
    
    
    
    
    inline size_t addr( size_t a, size_t b ) {
        return a + b * ma();
    }
    

    inline size_t saddr( ptrdiff_t a, ptrdiff_t b ) {
        return addr( a + 1, (b + 1));
    }
    
    inline size_t vaddr( size_t a, size_t b ) {
        return (a + b * ma()) * WIDTH;
    }
    

    inline size_t svaddr( ptrdiff_t a, ptrdiff_t b ) {
        return vaddr( a + 1, (b + 1));
    }
    
    inline pars_state_t get_seq_a ( size_t idx, size_t a ) {
        return pars_state_t(m_seqlist_a[idx][a * m_stride_a]);
    }
    
    inline int get_aux_a( size_t idx, size_t a ) {
        return m_auxlist_a[idx][a * m_stride_aux_a];
    }
    
    void align_freeshift_s11();
    void align_freeshift_s3() ;

    
public:
    pars_align_vec( const int** seqA, unsigned char* seqB, size_t n_a, size_t n_b, size_t aStride, 
                    const unsigned int **aAux, size_t aAuxStride, arrays<WIDTH> &arr, const unsigned int *bvtrans,
                    score_t gapOpen = 1, score_t gapExtend = 1, score_t mismatch = 3, score_t matchCGap = 10 );
    
    ~pars_align_vec() {}
    
    inline score_t *align_freeshift() {
        align_freeshift_s11();
        align_freeshift_s3();
        return get_best();
    }
        
    inline score_t *get_best() {
        return m_best.base();
    }
    inline size_t get_ncups() {
     
        return m_ncups;
    }
    
};


#endif
