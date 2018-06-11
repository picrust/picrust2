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

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <cassert>
#include <stdint.h>
#include <cstdlib>
#include <cstddef>
#include "align_vec.h"
#include "vec_unit.h"
#include "fasta.h"





struct seqs {
    std::vector<std::string> names;
    std::vector<std::vector<int> > data;
  
};

template <class score_t>
score_t align( std::string &a, std::string &b, scoring_matrix &m ) {
 
    std::vector<score_t> s( a.size());
    std::vector<score_t> si( a.size());
    
    std::fill( s.begin(), s.end(), 0 );
    std::fill( si.begin(), si.end(), 0 );
    const score_t SMALL = -32000;
    score_t max = SMALL;
    
    for( int ib = 0; ib < b.size(); ib++ ) {
        char bc = b[ib];
        
        score_t last_sl = SMALL;
        score_t last_sc = 0;
        score_t last_sdiag = 0;
        
       
        for( int ia = 0; ia < a.size(); ia++ ) {
            char ac = a[ia];
            score_t match = m.get_score( ac, bc );
            
            score_t sm = last_sdiag + match;
            const score_t GAP_EXT = -1;
            const score_t GAP_OPEN = -2;
            
            score_t sl = last_sl + GAP_EXT;
            sl = std::max( sl, score_t(last_sc + GAP_OPEN) );
            last_sl = sl;
            
            last_sdiag = s[ia];
            score_t su = si[ia] + GAP_EXT;
            su = std::max( su, score_t(last_sdiag + GAP_OPEN) );
            si[ia] = su;
            
            score_t sc = std::max( sm , std::max( std::max( sl, su ), score_t(0) ) );
            
            s[ia] = sc;
            last_sc = sc;
            max = std::max( sc, max );
        }
        
    }
    
    return max;
}



// void test() {
//  
//     aligned_buffer<short> a(16);
//     aligned_buffer<short> b(16);
//     aligned_buffer<short> c(16);
//     
//     std::fill( a.begin(), a.end(), 1 );
//     std::fill( b.begin(), b.end(), 2 );
//     a.m_ptr[0] = 2;
//     
//     typedef vector_unit<short,16> vu;
//     
//     vu::vec_t av = vu::load(a.m_ptr);
//     vu::vec_t bv = vu::load(b.m_ptr);
//     
//     vu::store( vu::cmp_lt(av,bv), c.m_ptr );
//     
//     for( int i = 0; i < c.size(); i++ ) {
//      
//         std::cout << i << ": " << c.m_ptr[i] << "\n";
//     }
//     
// }
int main( int argc, char *argv[] ) {
//     test();
//     return 0;
    
    
    seqs sq;
    seqs sd;
    
   // typedef mapped_file file_input;
    typedef std::ifstream file_input;
    
    if( argc != 6 ) {
     
        throw std::runtime_error( "missing parameters. expect: <open> <ext> <query.fa> <db.fa> <matrix>" );
        
    }
    
    const int gap_open = atoi( argv[1] );
    const int gap_extend = atoi( argv[2] );
    
    //std::cout << "gap: " <<

    std::ifstream ism( argv[5] );
    scoring_matrix sm( ism );
    
    file_input qfi( argv[3] );
    inc_fasta<file_input, scoring_matrix> qfasta( qfi, sm );
    
    file_input dfi( argv[4] );
    inc_fasta<file_input, scoring_matrix> dfasta( dfi, sm );
        
        
    
    
    
    
//     for( int i = 0 ; i < sd.names.size(); i++ ) {
//         int sc = align<int>( sd.data[i], sq.data[0], sm );
//         std::cout << sc << "\t" << sd.names[i] << "\n";
//     }
    
    
#if 0
    const int W = 8;
    typedef short score_t;
    typedef short sscore_t;
#else
    const int W = 16;
    typedef unsigned char score_t;
    typedef char sscore_t;
#endif
    persistent_state<score_t> ps;
//     size_t db_size = (sd.names.size() / W ) * W;
    
#ifndef LOCAL_ALIGN
#define PAD_FRONT  
#endif
    
    typedef uint8_t seq_char_t;
    
    std::string dname[W];
    std::vector<seq_char_t> ddata[W];
//     size_t dpad[W];
    aligned_buffer<seq_char_t> ddata_int;
    
//     std::vector<score_t> dmask[W];
    
    std::string qname;
    std::vector<seq_char_t> qdata;
    
    bool have_input = true;
    size_t n_qseq = 0;
    size_t n_qchar = 0;
    
    size_t n_dseq = 0;
    size_t n_dchar = 0;
    
    bool first_block = true;
    while( have_input ) {
        
        size_t maxlen = 0;
        
        int lj = -1;
        
        // determine db sequences for the current block
        
        for( int j = 0; j < W; j++ ) {
            dname[j].resize(0);
            ddata[j].resize(0);
            have_input = dfasta.next_seq( dname[j], ddata[j] );
            
            
            // if there aren't enough db sequences left to fill the block, pad with last db sequence
            if( !have_input ) {
                
                // break immediately if there are no db sequences left (means #db-seqs % W == 0, or otherwise have_input would have been == false from last iteration)
                if( j == 0 ) {
                    break;
                } else {
                    dname[j] = dname[lj];
                    ddata[j] = ddata[lj];
                }
            } else {
                n_dseq++;
                n_dchar += ddata[j].size();
                lj = j; // store largest valid 'j'
                
//                 for( int i = 0; i < ddata[j].length(); i++ ) {
//                         
//                     ddata[j][i] = sm.state_backmap(ddata[j][i]);
//                 }
            }
            
//             dmask[j].clear();
//             dmask[j].resize(ddata[j].length(), 0xffff );
            
            
            maxlen = std::max( maxlen, ddata[j].size() );
        }
        
        // jf == -1 at this point means that the block is empty (#db-seqs % W == 0)
        if( lj == -1 ) {
            break;
        }
        
        
        // WARNING: front-padding is currently defunct
#ifdef PAD_FRONT
        for( int j = 0; j < W; j++ ) {
            dpad[j] = maxlen - ddata[j].size();
        }
#endif

        // lj == -1 means that there are no remaining sequences
        
        
        
//         std::cout << "sdis: " << sdi.size() << "\n";
        aligned_buffer<sscore_t> qprofile( maxlen * W * sm.num_states());
        aligned_buffer<sscore_t>::iterator qpi = qprofile.begin();
#if 0
        for( int j = 0; j < sm.num_states(); j++ ) {
            const int jstate = sm.get_state(j);
//             const scoring_matrix::score_t *slice = sm.get_slice(jstate);
            const scoring_matrix::score_t *cslice = sm.get_cslice(j);
            for( int k = 0; k < maxlen; k++ ) {
                for( int l = 0; l < W; l++ ) {
                    
                    
                    std::vector<seq_char_t> &sdi = ddata[l];
                    
                    
#ifdef PAD_FRONT
                    const size_t frontpad = dpad[l];//maxlen - sdi.size();

                    //if( k < sdi.size() ) {
                    if( k >= frontpad )  {
                        //*qpi = sm.get_score( sdi[k], jstate);
                        *qpi = cslice[sdi[k - frontpad]];
#else
                    if( k < sdi.size() ) {
//                         *qpi = slice[sdi[k]];
                        
                        *qpi = cslice[sdi[k]];
#endif
                    } else {
                        *qpi = 0;   
                    }
//                     std::cout << "prof: " << (qpi - qprofile.begin() ) << " " << *qpi << " " << int(j) << " " << char(sm.get_state(j)) << " " << sdi[k] << "\n";
                    qpi++;
                }
            }
        }
#else
        // setup the qprofile (= lookup table for match penalties along the db-sequences in the current block)
        // this is the faster (at least on core i5) two-step version, using interleaved db-sequences
        
        // setup buffer for interleaved db sequences
        if( ddata_int.size() < maxlen * W ) {
            ddata_int.resize(maxlen * W);
        }
 
        // copy individual db sequences into interleaved buffer (padding the shorter sequnences 
        aligned_buffer<seq_char_t>::iterator dint_iter = ddata_int.begin();
        const int zero_state = sm.get_zero_state();
        for( size_t i = 0; i < maxlen; i++ ) {
            for( int j = 0; j < W; j++ ) {
                std::vector<seq_char_t> &sdi = ddata[j];
                if( i < sdi.size() ) {
                    *dint_iter = sdi[i];

                } else {
                    *dint_iter = zero_state;   
                }
                ++dint_iter;
            }
        }
        for( size_t j = 0; j < sm.num_states(); j++ ) {
            dint_iter = ddata_int.begin();
            const char *cslice = sm.get_cslice(j);
            for( size_t k = 0; k < maxlen; k++ ) {
                for( int l = 0; l < W; l++ ) {
//                     if( *dint_iter == zero_state ) {
//                         std::cout << int(cslice[*dint_iter]) << "\n";
//                        
//                     }
                    
                    *qpi = cslice[*dint_iter];
                    ++dint_iter;
                    ++qpi;
                }
            }
        }
#endif
//         std::cout << "sdis2: " << sdi.size() << std::endl;
//         aligned_buffer<short> dv(sdi.size() * W);
//         short *dv_iter = dv.begin();
//         
//         for( int j = 0; j < sdi.size(); j++ ) {
//                 
//             for( int k = 0; k < W; k++ ) {
//                 *dv_iter = sdi[j];
//                 dv_iter++;
//             }
//             
//         }
        
        
//         aligned_buffer<score_t> len(W);
//         for( int j = 0; j < W; j++ ) {
//             len.m_ptr[j] = sd.data[i + j].size();
//         }
        
        std::vector<int> out(W);
        
        qfasta.reset();
        
        while( qfasta.next_seq( qname, qdata )) {
            if( first_block ) {
                n_qseq++;
                n_qchar+=qdata.size();
            }
            
            align_vec<score_t,sscore_t,W>( ps, maxlen, qdata, sm, qprofile, gap_open, gap_extend, out );
        
            for( int j = 0; j <= lj; j++ ) {
//                 std::cout << out[j] << "\t" << dname[j] << " " << qname << " " << ddata[j].size() << "\n";
                std::cout << out[j] << "\t" << dname[j] << "\n";
            }
            qname.resize(0);
            qdata.resize(0);
            
        }
        first_block = false;
    }
    
    std::cerr << n_qseq << "[" << n_qchar << "] x " << n_dseq << "[" << n_dchar << "]\n";
}