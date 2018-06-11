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

#include <iostream>
#include <fstream>


#include <algorithm>
#include <deque>


//#include <boost/bind.hpp>
#include "align_vec.h"
// #define USE_BOOST_THREADS
#ifdef USE_BOOST_THREADS
#define BOOST_LIB_DIAGNOSTIC
#include <boost/thread.hpp>

namespace timpl = boost;

#else
#include "ivymike/thread.h"
namespace timpl = ivy_mike;
#endif


#include "ivymike/time.h"
#include "ivymike/write_png.h"
#include "ivymike/tdmatrix.h"
#include "ivymike/fasta.h"
#include "ivymike/aligned_buffer.h"

#ifndef PWDIST_INLINE 
// this means this file is not included by pairwise_seq_distance.h itself...
#include "pairwise_seq_distance.h"
#endif




#ifndef PSD_DECLARE_INLINE
#define PSD_DECLARE_INLINE
#endif


using ivy_mike::scoring_matrix;
using ivy_mike::aligned_buffer;

//typedef boost::multi_array<int,2> pw_score_matrix;
typedef ivy_mike::tdmatrix<int> pw_score_matrix;

static double read_temp() {
    std::ifstream is("/sys/class/hwmon/hwmon0/temp1_input" );
    
    if( is.good() ) {
        int temp;
        is >> temp;
        
        return temp * 1e-3;
    } else {
        return -1;
    }
}

template <size_t W, typename seq_char_t>
struct db_block {
    size_t didx[W];
//     std::vector<seq_char_t> *ddata[W];
    size_t dpad[W];    
    size_t maxlen;
    int lj;
};

template <typename block>
struct block_queue {
    std::deque<block> m_blocks;
    timpl::mutex m_mtx;
    volatile size_t m_ncup;
    volatile size_t m_ok_flags;
    block_queue() : m_ncup(0), m_ok_flags(0) {}
};

template <typename block_t>
struct worker {
    block_queue<block_t> &m_queue;
    
    worker( block_queue<block_t> &q ) : m_queue(q) {}
    void operator()() {
        
    }
};


// alignment worker thread. consumes block objects from the block-queue and writes results to the 2d matrix (m_outscore)

template <size_t W, typename seq_char_t, typename score_t, typename sscore_t>
struct lworker {
    typedef db_block<W, seq_char_t> block_t;
    const size_t m_nthreads;
    const size_t m_rank;
    block_queue<block_t> &m_queue;
    const scoring_matrix &m_sm;
    const std::vector< std::vector<uint8_t> > &m_seq1;
    const std::vector< std::vector<uint8_t> > &m_seq2;
    const sscore_t gap_open;
    const sscore_t gap_extend;
    const size_t m_block_size;
    
    pw_score_matrix &m_outscore;
    const bool m_half_matrix;
    lworker( size_t nthreads, size_t rank, block_queue<block_t>&q, const scoring_matrix &sm, const std::vector< std::vector<uint8_t> > &seq1_, const std::vector< std::vector<uint8_t> > &seq2_, const sscore_t gap_open_, const sscore_t gap_extend_,pw_score_matrix &outscore, bool half_matrix, size_t block_size ) 
    : m_nthreads(nthreads), m_rank(rank), m_queue(q), m_sm(sm), m_seq1(seq1_), m_seq2(seq2_), gap_open(gap_open_), gap_extend(gap_extend_), m_block_size(block_size), m_outscore(outscore), m_half_matrix( half_matrix ) 
    {
        if( m_half_matrix ) {
            if( m_seq1.size() != m_seq2.size() ) {
                throw std::runtime_error( "half_matrix mode set with m_seq1.size() != m_seq2.size()." );
            }
        }
    }
    
    void operator()() {
        
        // thread entry point
        
        size_t n_qseq = 0;
        size_t n_qchar = 0;
        
        size_t n_dseq = 0;
        size_t n_dchar = 0;
        
        bool first_block = true;
        aligned_buffer<seq_char_t> ddata_int;
        persistent_state<score_t> ps;
        persistent_state_blocked<score_t, sscore_t> ps_blocked;
    
//         {
//             cpu_set_t cs;
//             CPU_ZERO( &cs );
//             CPU_SET( 0, &cs );
//             if(pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cs) != 0)
//             {
//                 printf("\n\nThere was a problem finding a physical core for thread number %d to run on.\n", 0);
//                 
//                 assert(0);
//             }
//             
//             
//         }
        
//         for( int i = 0; i < m_queue.m_blocks.size(); i++ ) {
//             
//             if( (i % m_nthreads) != m_rank ) {
//                 continue;
//             }
//         
//         
//             block_t block = m_queue.m_blocks[i];
        size_t ncups = 0;
        
        ivy_mike::timer t1;
        ivy_mike::timer t2;
        size_t ncups_last = 0;
        while(true) {
            // get next block from the queue
            
            block_t block;
            {
                timpl::lock_guard<timpl::mutex> lock( m_queue.m_mtx );
                
                if ( m_queue.m_blocks.empty() ) {
                    break;
                }
                block = m_queue.m_blocks.front();
                m_queue.m_blocks.pop_front();
            }
            
            

            aligned_buffer<sscore_t> qprofile( block.maxlen * W * m_sm.num_states());
            typename aligned_buffer<sscore_t>::iterator qpi = qprofile.begin();
            
            // setup the qprofile (= lookup table for match penalties along the db-sequences in the current block)
            // this is the faster (at least on core i5) two-step version, using interleaved db-sequences
            
            // setup buffer for interleaved db sequences
            if ( ddata_int.size() < block.maxlen * W ) {
                ddata_int.resize(block.maxlen * W);
            }
            
            // copy individual db sequences into interleaved buffer (padding the shorter sequnences
            typename aligned_buffer<seq_char_t>::iterator dint_iter = ddata_int.begin();
            const int zero_state = m_sm.get_zero_state();
            for ( size_t i = 0; i < block.maxlen; i++ ) {
                for ( size_t j = 0; j < W; j++ ) {
                    const std::vector<seq_char_t> &sdi = m_seq1.at(block.didx[j]);//*(block.ddata[j]);
                    if ( i < sdi.size() ) {
                        *dint_iter = sdi[i];
                        // the aligner will catch this later if assertions are enabled
// #ifdef DEBUG
//                         if( *dint_iter >= m_sm.num_states() ) {
//                             throw std::runtime_error( "meeeep. illegal character in input sequences\n" );
//                         }
// #endif                   
                    } else {
                        *dint_iter = zero_state;
                    }
                    
                    
                    
                    //                 std::cout << j << " " << int(*dint_iter) << " " << (i < sdi.size()) << "\n";
                    ++dint_iter;
                }
            }
            
            //copy interleaved scoring-matrix
            for ( size_t j = 0; j < m_sm.num_states(); j++ ) {
                dint_iter = ddata_int.begin();
                const char *cslice = m_sm.get_cslice(j);
                for ( size_t k = 0; k < block.maxlen; k++ ) {
                    for ( size_t l = 0; l < W; l++ ) {
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
            
            std::vector<int> out(W);
            size_t i_max;
            if( m_half_matrix ) {
                i_max = block.didx[block.lj];
            } else {
                i_max = m_seq2.size() - 1;
            }
//             std::cout << "i_max: " << i_max << " " << block.maxlen << "\n";
            
//             const size_t i_max = m_seq.size() - 1;
            
            // loop over all sequences and align them against the current profile
            
            for ( size_t i_seq2 = 0; i_seq2 <= i_max; ++i_seq2 ) {
//             for ( size_t i_seq2 = 0; i_seq2 < m_seq.size(); ++i_seq2 ) {
                const std::vector<uint8_t> &qdata = m_seq2.at(i_seq2);
                
                if ( first_block ) {
                    n_qseq++;
                    n_qchar+=qdata.size();
                }
                
                
                // call the alignment kernel 
                
                if( m_block_size == 0 ) {
                    align_vec<score_t,sscore_t,W>( ps, block.maxlen, qdata, m_sm, qprofile, gap_open, gap_extend, out );
                } else if( m_block_size == 32 ) {
                    align_vec_blocked<score_t,sscore_t,W,32>( ps_blocked, block.maxlen, qdata, m_sm, qprofile, gap_open, gap_extend, out );
                } else if( m_block_size == 64 ) {
                    align_vec_blocked<score_t,sscore_t,W,64>( ps_blocked, block.maxlen, qdata, m_sm, qprofile, gap_open, gap_extend, out );
                } else if( m_block_size == 128 ) {
                    align_vec_blocked<score_t,sscore_t,W,128>( ps_blocked, block.maxlen, qdata, m_sm, qprofile, gap_open, gap_extend, out );
                } else if( m_block_size == 256 ) {
                    align_vec_blocked<score_t,sscore_t,W,256>( ps_blocked, block.maxlen, qdata, m_sm, qprofile, gap_open, gap_extend, out );
                } else {
                    std::cerr << "worker thread abort: unsupported block size: " << m_block_size << "\n";
                    return;
                }
                
                // write output scores to the output matrix. no lock necessary, as writes are independent.
                
                for ( int j = 0; j <= block.lj; j++ ) {
                    //                 std::cout << out[j] << "\t" << dname[j] << " " << qname << " " << ddata[j].size() << "\n";
                    //                     std::cout << out[j] << "\t" << block.didx[j] << " " << i_seq2 << "\n";
//                     m_outscore[block.didx[j]][i_seq2] = out[j];
                    
                    
                    m_outscore[block.didx[j]][i_seq2] = out[j];
                    if( m_half_matrix ) { 
                        m_outscore[i_seq2][block.didx[j]] = out[j];
                    }
                    ncups += m_seq2[i_seq2].size() * m_seq1[block.didx[j]].size();
                }
                
                
                
            }
            
            for ( int j = 0; j <= block.lj; j++ ) {
                    //                 std::cout << out[j] << "\t" << dname[j] << " " << qname << " " << ddata[j].size() << "\n";
//                     std::cout << out[j] << "\t" << block.didx[j] << " " << i_seq2 << "\n";
                    
                n_dseq++;
                n_dchar += m_seq1[block.didx[j]].size();
            }
            first_block = false;
            if( m_rank == 0 && t1.elapsed() > 2 ) {
                size_t dncup = ncups - ncups_last;
                
                std::cerr << t2.elapsed() << " " << dncup << " in " << t1.elapsed() << " s " << dncup / (t1.elapsed() * 1e6) << " " << read_temp() << std::endl;
                t1 = ivy_mike::timer();
                ncups_last = ncups;
            }
        }
        
        {
            std::cerr << n_qchar << " x " << n_dchar << "\n";
            timpl::lock_guard<timpl::mutex> lock( m_queue.m_mtx );
            m_queue.m_ncup += ncups;
            m_queue.m_ok_flags++;
        }
        
        
        
    }
};


// WARNING: the sequences are expected to be transformed to 'compressed states' (= 0, 1, 2 ...) rather than characters.
// The state mapping must be consistent with the supplied scoring matrix and its compressed form.
// Sequences containing numbers >= sm.num_states() will likely blow up the aligner, as there are no checks after this point!
PSD_DECLARE_INLINE bool pairwise_seq_distance( const std::vector< std::vector<uint8_t> > &seq1, const std::vector< std::vector<uint8_t> > &seq2, bool identical, pw_score_matrix &out_scores, const scoring_matrix &sm, const int gap_open, const int gap_extend, const size_t n_thread, const size_t block_size ) {
#if 1
    const int W = 8;
    typedef short score_t;
    typedef short sscore_t;
#else
    const int W = 16;
    typedef unsigned char score_t;
    typedef char sscore_t;
#endif
    
//     size_t db_size = (sd.names.size() / W ) * W;
    

    ivy_mike::timer t1;
    
    
//     std::vector< std::vector<uint8_t> > seq( seq_raw.size() );
//     seq.resize(400);
//     for( int i = 0; i < seq.size(); i++ ) {
//         std::for_each( seq_raw[i].begin(), seq_raw[i].end(), scoring_matrix::valid_state_appender<std::vector<uint8_t> >(sm, seq[i]) );
//     }
    
    if( seq1.size() != out_scores.size() || seq2.size() != out_scores[0].size() ) {
        throw std::runtime_error( "out_scores matrix is too small" );
    }
    
//     const sscore_t gap_open = -5;
//     const sscore_t gap_extend = -2;
    
    typedef uint8_t seq_char_t;
    
    
   
    //std::string dname[W];
    
    
    
//     std::vector<score_t> dmask[W];
    
    
    
    bool have_input = true;
    
    
    size_t i_seq1 = 0;
    
    // TODO: update comment for seq1 * seq2 alignment!
    // the following code basically consists of two nested loops which align all elements in seq against each other (N*N alignments).
    
    // It is a bit hard to recognize, though as the alignments operations are distributed to 'blocks' (=independent work units)
    // consumed by the worker threads.

    // each block normally consists of W (=vector unit width) sequences to be aligned agains all other sequences
    block_queue<db_block<W, seq_char_t> > q;
    std::deque<db_block<W, seq_char_t> > &blocks = q.m_blocks;
    
    // generate the block objects and put them in the queue.
    while( have_input ) {
        
        // determine db sequences for the current block
        db_block<W, seq_char_t> block;
        block.maxlen = 0;
        
        block.lj = -1;
        for( int j = 0; j < W; j++ ) {
//             dname[j].resize(0);
            //ddata[j].resize(0);
//             have_input = (i_seq1 != seq.size());                       
           have_input = (i_seq1 != seq1.size());
//             have_input = i_seq1 < 30;
           // std::cout << "have_input " << have_input << " " << seq.size() << "\n";
            
            
            
            // if there aren't enough db sequences left to fill the block, pad with last db sequence
            if( !have_input ) {
                
                // break immediately if there are no db sequences left (means #db-seqs % W == 0, or otherwise have_input would have been == false from last iteration)
                if( j == 0 ) {
                    break;
                } else {
//                     block.ddata[j] = block.ddata[block.lj];
                    block.didx[j] = block.didx[block.lj];
                }
            } else {
                block.didx[j] = i_seq1;
//                 block.ddata[j] = &seq[i_seq1];
                ++i_seq1;
                
                
                block.lj = j; // store largest valid 'j'
                
//                 for( int i = 0; i < ddata[j].length(); i++ ) {
//                         
//                     ddata[j][i] = sm.state_backmap(ddata[j][i]);
//                 }
            }
            
//             dmask[j].clear();
//             dmask[j].resize(ddata[j].length(), 0xffff );
            
            
            block.maxlen = std::max( block.maxlen, seq1[block.didx[j]].size() );
        }
//         std::cout << "maxlen; " << block.maxlen << "\n";
        // jf == -1 at this point means that the block is empty (#db-seqs % W == 0)
        if( block.lj == -1 ) {
            break;
        }
        
//         std::cout << "block: " << block.didx[block.lj] << "\n";
        
        blocks.push_back(block);
        
    }
    
    std::cerr << "blocks: " << blocks.size() << "\n";
    //throw std::runtime_error( "exit" );
    
 
//     return;
    
    
    //pw_score_matrix out_scores(boost::extents[seq.size()][seq.size()]) ;
    
        
    
    // spawn the worker threads. Each of them will consume blocks from the block-queue until it is empty.
    // the results are concurrently written to the 2d matrix out_scores
    timpl::thread_group tg;
    
    for( size_t i = 0; i < n_thread; ++i ) {
        lworker<W, seq_char_t, score_t, sscore_t> lw( n_thread, i, q, sm, seq1, seq2, gap_open, gap_extend, out_scores, identical, block_size );
        
        std::cerr << "thread " << i << "\n";
        
        tg.create_thread( lw );
    }
    
        
    tg.join_all();

    if( q.m_ok_flags != n_thread ) {
        std::cerr << n_thread - q.m_ok_flags << " threads did not exit properly.\n";
        return false;
    }
    
    
    std::cerr << "aligned " << seq1.size() << " x " << seq2.size() << " sequences. " << q.m_ncup << " " << (q.m_ncup / (t1.elapsed() * 1.0e9)) << " GCup/s\n";
    return true;

}
