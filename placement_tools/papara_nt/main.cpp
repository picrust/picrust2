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


#include "ivymike/multiple_alignment.h"
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <vector>
#include <deque>
#include <map>
#include <functional>
#include <cstring>
#include <boost/io/ios_state.hpp>

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/bind.hpp>




#include "parsimony.h"
#include "pvec.h"
#include "align_utils.h"

#include "pars_align_seq.h"
#include "pars_align_gapp_seq.h"
#include "fasta.h"
#include "vec_unit.h"
#include "align_pvec_vec.h"

#include "ivymike/tree_parser.h"
#include "ivymike/time.h"
#include "ivymike/getopt.h"
#include "ivymike/thread.h"
#include "ivymike/demangle.h"
#include "ivymike/stupid_ptr.h"
#include "ivymike/algorithm.h"
#include "ivymike/smart_ptr.h"

using namespace ivy_mike;
using namespace ivy_mike::tree_parser_ms;

using namespace boost::numeric;



// template<typename score_t>
// struct align_arrays_vec {
//     aligned_buffer<score_t> s;
//     
// };
// 
// template<typename score_t, size_t W>
// static void align_pvec_score( aligned_buffer<score_t> &a_prof, aligned_buffer<score_t> &a_aux_prof, const std::vector<uint8_t> &b, score_t mismatch_score, score_t match_cgap, score_t gap_open, score_t gap_extend, align_arrays_vec<score_t> &arr, aligned_buffer<score_t> &out ) {
// 
//     typedef vector_unit<score_t, W> vu;
//     typedef typename vu::vec_t vec_t;
//     
//     if( arr.s.size() < a_prof.size() ) {
//         arr.s.resize( a_prof.size() );
//         
//     }
//     
//     const score_t aux_cgap = 0x1;
//     
//     assert( a_prof.size() % W == 0 );
//     size_t asize = a_prof.size() / W;
//     assert( asize > b.size() );
//     
//     const size_t band_width = asize - b.size();
//     std::fill( arr.s.begin(), arr.s.end(), 0 );
//     
//     const score_t LARGE = 32000;
// //     std::fill( arr.si.begin(), arr.si.end(), LARGE );
//     std::fill( arr.s.begin() + (band_width + 1) * W, arr.s.begin() + (band_width + 2) * W, LARGE );
//     //arr.s.at(band_width+1) = LARGE;
// 
// //     score_t min_score = LARGE;
//     
//     for( size_t ib = 0; ib < b.size(); ib++ ) {
//         const vec_t bc = vu::set1(b[ib]);
//         
//         vec_t last_sl = vu::set1(LARGE);
//         vec_t last_sc = vu::set1(LARGE);
//         vec_t last_sdiag = vu::set1(0);
//         
//         size_t astart = ib;
//                 
//         score_t * __restrict s_iter = &arr.s[0];
// 
//         
//         
//         
//        
//         last_sdiag = vu::load(s_iter);
//         
//         
//         for( size_t ia = astart; ia <= ib + band_width; ++ia, s_iter += W ) {  
//             vec_t ac = vu::load( a_prof(ia * W) );
//             vec_t aaux = vu::load(a_aux_prof(ia * W));
//             
//             vec_t cgap_mask = vu::cmp_eq( aaux, vu::set1(aux_cgap));
//             
//             
//             // determine match or mis-match according to parsimony bits coming from the tree.
//             
//             vec_t mismatch_mask = vu::cmp_eq( vu::bit_and( ac, bc ), vu::setzero() );
//             
//             vec_t sm = vu::add( last_sdiag, vu::bit_and( vu::set1(mismatch_score), mismatch_mask ));
//             sm = vu::add( sm, vu::bit_and( vu::set1( match_cgap ), cgap_mask ));
//             
//             last_sdiag = vu::load(s_iter+W);
// 
// //            cgap_mask = vu::bit_invert( cgap_mask );
//             
//             vec_t last_sc_OPEN = vu::add( last_sc, vu::bit_andnot( cgap_mask, vu::set1(gap_open + gap_extend)));
//             vec_t sl_score_stay = vu::add( last_sl, vu::bit_andnot( cgap_mask, vu::set1(gap_extend)));
//             
//                        
//             vec_t sl = vu::min( sl_score_stay, last_sc_OPEN );
//                         
//             last_sl = sl;
//             
// 
//             vec_t su = vu::add( last_sdiag, vu::set1(mismatch_score));
//             
// //             std::cout << "x: " << sm << " " << su << " " << sl << "  ";
//             vec_t sc = vu::min( sm, vu::min( su, sl ) );
//             
//             last_sc = sc;
//             
//             vu::store( sc, s_iter );
//             
//            
//         }
//     }
//     
//     vec_t minscore = vu::set1(LARGE);
//     for( size_t i = 0; i < band_width + 1; ++i ) {
//         minscore = vu::min( minscore, vu::load( arr.s(i * W)));
//     }
//     out.resize(W);
//     vu::store( minscore, out(0));
//     
// }

namespace {
    typedef boost::iostreams::tee_device<std::ostream, std::ofstream> log_device;
    typedef boost::iostreams::stream<log_device> log_stream;
    
    log_stream lout;
    
    template<typename stream_, typename device_>
    class bios_open_guard {
        stream_ &m_stream;
    public:
        bios_open_guard( stream_ &stream, device_ &device ) : m_stream(stream) {
            m_stream.open( device );
        }
        ~bios_open_guard() {
            m_stream.close();
        }
    };
    
    typedef bios_open_guard<log_stream, log_device> log_stream_guard;
}

class ostream_test {
    std::ostream &m_os;

public:
    ostream_test( std::ostream &os ) : m_os(os) {}
    void operator()(int i) {
        m_os << i;
    }
};



// static inline parsimony_state dna_to_parsimony_state( uint8_t c ) {
//     switch( c ) {
//         case 'A':
//         case 'a':
//             return 0x1;
//
//         case 'C':
//         case 'c':
//             return 0x2;
//
//         case 'G':
//         case 'g':
//             return 0x4;
//
//         case 'U':
//         case 'u':
//         case 'T':
//         case 't':
//             return 0x8;
//
//         default:
//             return 0xf;
//     };
//
// }
//
//
// static inline int dna_to_cgap( uint8_t c ) {
//     switch( c ) {
//     case 'A':
//     case 'a':
//     case 'C':
//     case 'c':
//     case 'G':
//     case 'g':
//     case 'T':
//     case 't':
//     case 'U':
//     case 'u':
//         return 0x0;
//
//     default:
//         return AUX_CGAP;
//     };
// }

static bool g_dump_aux = false;


// probgap_model *pvec_pgap::pgap_model = 0;

// class pvec_ugly {
//     parsimony_state *m_v;
//     size_t m_size;
//
// public:
//     pvec_ugly() : m_v(0), m_size(0) {}
//     ~pvec_ugly() {
//         delete[] m_v;
//     }
//     void init( const std::vector<uint8_t> &seq ) {
//         assert( m_v == 0 );
//         delete[] m_v;
//         m_size = seq.size(0);
//         m_v = new parsimony_state[m_size];
//
//         std::transform( seq.begin(), seq.end(), m_v, dna_to_parsimony_state );
//     }
//
//     static void newview( pvec_ugly &p, pvec_ugly &c1, pvec_ugly &c2 ) {
//         assert( c1.m_size() == c2.m_size() );
//
// //         p.v.resize(0);
//         //p.v.resize(c1.v.size());
// #error continue here!
//
//         for( size_t i = 0; i < c1.v.size(); i++ ) {
//             parsimony_state ps = c1.v[i] & c2.v[i];
//
//             if( ps == 0 ) {
//                 ps = c1.v[i] | c2.v[i];
//             }
//
//             //p.v.push_back( ps );
//             p.v[i] = ps;
//
//         }
//     }
//
// };


template<class pvec_t>
class my_adata_gen : public ivy_mike::tree_parser_ms::adata {
//     static int ct;
    //std::vector<parsimony_state> m_pvec;
    pvec_t m_pvec;

public:
//     int m_ct;
    my_adata_gen() {

//         std::cout << "my_adata\n";

    }

    virtual ~my_adata_gen() {

//         std::cout << "~my_adata\n";

    }

    virtual void visit() {
//         std::cout << "tr: " << m_ct << "\n";
    }
    void init_pvec(const std::vector< uint8_t >& seq) {


        m_pvec.init( seq );
//         std::cout << "init_pvec: " << m_pvec.size() << "\n";
//                 m_pvec.reserve(seq.size());
//         for( std::vector< uint8_t >::const_iterator it = seq.begin(); it != seq.end(); ++it ) {
//             m_pvec.push_back(dna_to_parsimony_state(*it));
//
//         }
    }
    pvec_t &get_pvec() {
        return m_pvec;
    }


};

// inline void newview_parsimony( std::vector<parsimony_state> &p, const std::vector<parsimony_state> &c1, const std::vector<parsimony_state> &c2 ) {
//
// }



// inline std::ostream &operator<<( std::ostream &os, const my_adata &rb ) {
//
//     os << "my_adata: " << rb.m_ct;
// }

template<class ndata_t>
class my_fact : public ivy_mike::tree_parser_ms::node_data_factory {

    virtual ndata_t *alloc_adata() {

        return new ndata_t;
    }

};



template<class lnode>
void traverse_rec( lnode *n ) {

    n->m_data->visit();

    if( n->next->back != 0 ) {
        traverse_rec(n->next->back);
    }

    if( n->next->next->back != 0 ) {
        traverse_rec(n->next->next->back);
    }
}

// template<class lnode>
// void traverse( lnode *n ) {
//     n->m_data->visit();
//
//     if( n->back != 0 ) {
//         traverse_rec(n->back);
//     }
//
//     if( n->next->back != 0 ) {
//         traverse_rec(n->next->back);
//     }
//
//     if( n->next->next->back != 0 ) {
//         traverse_rec(n->next->next->back);
//     }
//
//
// }


uint8_t to_hex( double v ) {
    int vi = int(fabs(v));

    vi = std::min( vi, 15 );

    if( vi <= 9 ) {
        return '0' + vi;
    } else {
        return 'a' + (vi - 10);
    }

}


template<class pvec_t>
void do_newview( pvec_t &root_pvec, lnode *n1, lnode *n2, bool incremental ) {
    typedef my_adata_gen<pvec_t> my_adata;

    std::deque<rooted_bifurcation<lnode> > trav_order;

    //std::cout << "traversal for branch: " << *(n1->m_data) << " " << *(n2->m_data) << "\n";

    rooted_traversal_order( n1, n2, trav_order, incremental );
//     std::cout << "traversal: " << trav_order.size() << "\n";

    for( std::deque< rooted_bifurcation< ivy_mike::tree_parser_ms::lnode > >::iterator it = trav_order.begin(); it != trav_order.end(); ++it ) {
//         std::cout << *it << "\n";

        my_adata *p = dynamic_cast<my_adata *>( it->parent->m_data.get());
        my_adata *c1 = dynamic_cast<my_adata *>( it->child1->m_data.get());
        my_adata *c2 = dynamic_cast<my_adata *>( it->child2->m_data.get());
//         rooted_bifurcation<ivy_mike::tree_parser_ms::lnode>::tip_case tc = it->tc;

//         std::cout << "tip case: " << (*it) << "\n";
        pvec_t::newview(p->get_pvec(), c1->get_pvec(), c2->get_pvec(), it->child1->backLen, it->child2->backLen, it->tc);

    }





    {
        my_adata *c1 = dynamic_cast<my_adata *>( n1->m_data.get());
        my_adata *c2 = dynamic_cast<my_adata *>( n2->m_data.get());

//         tip_case tc;

        if( c1->isTip && c2->isTip ) {
//                 std::cout << "root: TIP TIP\n";
            pvec_t::newview(root_pvec, c1->get_pvec(), c2->get_pvec(), n1->backLen, n2->backLen, TIP_TIP );
        } else if( c1->isTip && !c2->isTip ) {
//                 std::cout << "root: TIP INNER\n";
            pvec_t::newview(root_pvec, c1->get_pvec(), c2->get_pvec(), n1->backLen, n2->backLen, TIP_INNER );
//             root_pvec = c2->get_pvec();
        } else if( !c1->isTip && c2->isTip ) {
//                 std::cout << "root: INNER TIP\n";
            pvec_t::newview(root_pvec, c2->get_pvec(), c1->get_pvec(), n1->backLen, n2->backLen, TIP_INNER );
//             root_pvec = c1->get_pvec();
        } else {
//                 std::cout << "root: INNER INNER\n";
            pvec_t::newview(root_pvec, c1->get_pvec(), c2->get_pvec(), n1->backLen, n2->backLen, INNER_INNER );
        }


    }
//     std::cout << std::hex;
//     for( std::vector< parsimony_state >::const_iterator it = root_pvec.begin(); it != root_pvec.end(); ++it ) {
//         std::cout << *it;
//     }
//
//     std::cout << std::dec << std::endl;

}


static void seq_to_nongappy_pvec( std::vector<uint8_t> &seq, std::vector<uint8_t> &pvec ) {
    pvec.resize( 0 );

    for( unsigned int i = 0; i < seq.size(); i++ ) {
        uint8_t ps = dna_parsimony_mapping::d2p(seq[i]);

        if( ps == 0x1 || ps == 0x2 || ps == 0x4 || ps == 0x8 ) {
            pvec.push_back(ps);
        }

    }

}

void pairwise_seq_distance( std::vector< std::vector<uint8_t> > &seq );

class papara_nt_i : boost::noncopyable{
public:
    virtual ~papara_nt_i() {}

    virtual void calc_scores( size_t ) = 0;
    virtual void print_best_scores( std::ostream & ) = 0;
    virtual void write_result_phylip( std::ostream &, std::ostream &) = 0;
};
struct scoring_result {
	size_t edge;
	size_t qs;

	int res;

	scoring_result( size_t edge_, size_t qs_, int res_ ) : edge( edge_), qs(qs_), res(res_) {}
};

template<typename pvec_t>
class papara_nt : public papara_nt_i {

    const static int score_gap_open = 3;
    const static int score_gap_extend = 1;
    const static int score_mismatch = 1;
    const static int score_match_cgap = 3;


//    const static int score_gap_open = 1;
//    const static int score_gap_extend = 1;
//    const static int score_mismatch = 3;
//    const static int score_match_cgap = 10;

    //typedef pvec_pgap pvec_t;
    typedef my_adata_gen<pvec_t> my_adata;


    const static size_t VW = 8;

    struct block_t {
        block_t() {
            memset( this, 0, sizeof( block_t )); // FIXME: hmm, this is still legal?
        }

        // WARNING: these are pointers into m_ref_pvecs and m_ref_aux
        // make sure they stay valid!
        const int *seqptrs[VW];
        const unsigned int *auxptrs[VW];
        const double *gapp_ptrs[VW];
        size_t ref_len;
        size_t edges[VW];
        int num_valid;
    };

//    papara_nt( const papara_nt &other );
//    papara_nt & operator=( const papara_nt &other );



    //multiple_alignment m_ref_ma;

    std::vector <std::string > m_ref_names;
    std::vector <std::vector<uint8_t> > m_ref_seqs;
    std::auto_ptr<ivy_mike::tree_parser_ms::ln_pool> m_ln_pool;
    edge_collector<lnode> m_ec;

    std::vector <std::string> m_qs_names;
    std::vector <std::vector<uint8_t> > m_qs_seqs;

    std::vector<std::vector <uint8_t> > m_qs_pvecs;

    std::vector<std::vector <int> > m_ref_pvecs;
    std::vector<std::vector <unsigned int> > m_ref_aux;
    std::vector<std::vector <double> > m_ref_gapp;

    ivy_mike::mutex m_qmtx; // mutex for the block queue and the qs best score/edge arrays
    std::deque<block_t> m_blockqueue;
    std::vector <int> m_qs_bestscore;
    std::vector <int> m_qs_bestedge;




    class worker {
        papara_nt &m_pnt;
        size_t m_rank;
    public:
        worker( papara_nt & pnt, size_t rank ) : m_pnt(pnt), m_rank(rank) {}
        void operator()() {

            pars_align_seq<>::arrays seq_arrays(true);
            pars_align_gapp_seq::arrays seq_arrays_gapp(true);
            
            ivy_mike::timer tstatus;
            double last_tstatus = 0.0;
   
            uint64_t ncup = 0;
            while( true ) {
                block_t block;

                {
                    ivy_mike::lock_guard<ivy_mike::mutex> lock( m_pnt.m_qmtx );
                    if( m_pnt.m_blockqueue.empty() ) {
                        break;
                    }
                    block = m_pnt.m_blockqueue.front();
                    m_pnt.m_blockqueue.pop_front();

                    if( m_rank == 0 && (tstatus.elapsed() - last_tstatus) > 10.0 ) {
                        std::cerr << tstatus.elapsed() << " " << m_pnt.m_blockqueue.size() << " blocks remaining. " << ncup / (tstatus.elapsed() * 1e9) << " gncup/s\n";
                        last_tstatus = tstatus.elapsed();
                    }
                }
#if 1
                assert( VW == 8 );
                const align_pvec_score<short,8> aligner( block.seqptrs, block.auxptrs, block.ref_len, score_mismatch, score_match_cgap, score_gap_open, score_gap_extend );
                for( unsigned int i = 0; i < m_pnt.m_qs_names.size(); i++ ) {

                    size_t stride = 1;
                    size_t aux_stride = 1;
                    
                    aligner.align(m_pnt.m_qs_pvecs[i]);
                    const short *score_vec = aligner.get_scores();

                    ncup += block.num_valid * block.ref_len * m_pnt.m_qs_pvecs[i].size();
                    {
                        ivy_mike::lock_guard<ivy_mike::mutex> lock( m_pnt.m_qmtx );

                        for( int k = 0; k < block.num_valid; k++ ) {



                            if( score_vec[k] < m_pnt.m_qs_bestscore[i] || (score_vec[k] == m_pnt.m_qs_bestscore[i] && block.edges[k] < size_t(m_pnt.m_qs_bestedge[i]) )) {
                                const bool validate = false;
                                if( validate ) {
                                    const int *seqptr = block.seqptrs[k];
                                    const unsigned int *auxptr = block.auxptrs[k];

                                    pars_align_seq<> pas( seqptr, m_pnt.m_qs_pvecs[i].data(), block.ref_len, m_pnt.m_qs_pvecs[i].size(), stride, auxptr, aux_stride, seq_arrays, 0, score_gap_open, score_gap_extend, score_mismatch, score_match_cgap );
                                    int res = pas.alignFreeshift(INT_MAX);

                                    if( res != score_vec[k] ) {


                                        std::cout << "meeeeeeep! score: " << score_vec[k] << " " << res << "\n";
                                    }
                                }

                                m_pnt.m_qs_bestscore[i] = score_vec[k];

								assert( block.edges[k] <= size_t(std::numeric_limits<int>::max()) );
                                m_pnt.m_qs_bestedge[i] = int(block.edges[k]); // TODO: dto, maybe change to size_t
                            }
                        }
                    }
                }
            }
#else

            assert( block.gapp_ptrs[0] != 0 );
            assert( VW == 4 );
            const align_pvec_gapp_score<4> aligner( block.seqptrs, block.gapp_ptrs, block.ref_len, score_mismatch, score_match_cgap, score_gap_open, score_gap_extend );
            for( unsigned int i = 0; i < m_pnt.m_qs_names.size(); i++ ) {

                size_t stride = 1;
                size_t aux_stride = 1;

                aligner.align(m_pnt.m_qs_pvecs[i]);
                const float *score_vec = aligner.get_scores();

                ncup += block.num_valid * block.ref_len * m_pnt.m_qs_pvecs[i].size();
                {
                    ivy_mike::lock_guard<ivy_mike::mutex> lock( m_pnt.m_qmtx );

                    for( int k = 0; k < block.num_valid; k++ ) {



                        if( score_vec[k] < m_pnt.m_qs_bestscore[i] || (score_vec[k] == m_pnt.m_qs_bestscore[i] && block.edges[k] < m_pnt.m_qs_bestedge[i] )) {
                            const bool validate = false;
                            if( validate ) {
                                const int *seqptr = block.seqptrs[k];
                                const double *gapp_ptr = block.gapp_ptrs[k];

//                                std::vector<double> gapp_tmp(gapp_ptr, gapp_ptr + block.ref_len);


                                pars_align_gapp_seq pas( seqptr, m_pnt.m_qs_pvecs[i].data(), block.ref_len, m_pnt.m_qs_pvecs[i].size(), stride, gapp_ptr, aux_stride, seq_arrays_gapp, 0, score_gap_open, score_gap_extend, score_mismatch, score_match_cgap );
                                int res = pas.alignFreeshift(INT_MAX);

                                if( res != score_vec[k] ) {


                                    std::cout << "meeeeeeep! score: " << score_vec[k] << " " << res << "\n";
                                }
                            }

                            m_pnt.m_qs_bestscore[i] = score_vec[k];
                            m_pnt.m_qs_bestedge[i] = block.edges[k];
                        }
                    }
                }
            }
        }

#endif
            {
                ivy_mike::lock_guard<ivy_mike::mutex> lock( m_pnt.m_qmtx );
                std::cout << "thread " << m_rank << ": " << ncup / (tstatus.elapsed() * 1e9) << " gncup/s\n";
            }
        }
    };

    void build_block_queue() {
        // creates the list of ref-block to be consumed by the worker threads.  A ref-block onsists of N ancestral state sequences, where N='width of the vector unit'.
        // The vectorized alignment implementation will align a QS against a whole ref-block at a time, rather than a single ancestral state sequence as in the
        // sequencial algorithm.

        assert( m_blockqueue.empty() );
        size_t n_groups = (m_ec.m_edges.size() / VW);
        if( (m_ec.m_edges.size() % VW) != 0 ) {
            n_groups++;
        }


//         std::vector<int> seqlist[VW];
//         const int *seqptrs[VW];
//         std::vector<unsigned int> auxlist[VW];
//         const unsigned int *auxptrs[VW];



        for ( size_t j = 0; j < n_groups; j++ ) {
            int num_valid = 0;



            block_t block;

            for( unsigned int i = 0; i < VW; i++ ) {

                size_t edge = j * VW + i;
                if( edge < m_ec.m_edges.size()) {
                    block.edges[i] = edge;
                    block.num_valid++;

                    block.seqptrs[i] = m_ref_pvecs[edge].data();
                    block.auxptrs[i] = m_ref_aux[edge].data();

                    if( !m_ref_gapp[edge].empty() ) {
                    	block.gapp_ptrs[i] = m_ref_gapp[edge].data();
                    } else {
                    	block.gapp_ptrs[i] = 0;
                    }

                    block.ref_len = m_ref_pvecs[edge].size();
                    //                     do_newview( root_pvec, m_ec.m_edges[edge].first, m_ec.m_edges[edge].second, true );
//                     root_pvec.to_int_vec(seqlist[i]);
//                     root_pvec.to_aux_vec(auxlist[i]);
//
//                     seqptrs[i] = seqlist[i].data();
//                     auxptrs[i] = auxlist[i].data();

                    num_valid++;
                } else {
                    if( i < 1 ) {
                        std::cout << "edge: " << edge << " " << m_ec.m_edges.size() << std::endl;

                        throw std::runtime_error( "bad integer mathematics" );
                    }
                    block.edges[i] = block.edges[i-1];

                    block.seqptrs[i] = block.seqptrs[i-1];
                    block.auxptrs[i] = block.auxptrs[i-1];
                    block.gapp_ptrs[i] = block.gapp_ptrs[i-1];
                }

            }
            m_blockqueue.push_back(block);
        }
    }

    void build_ref_vecs() {
        // pre-create the ancestral state vectors. This step is necessary for the threaded version, because otherwise, each
        // thread would need an independent copy of the tree to do concurrent newviews. Anyway, having a copy of the tree
        // in each thread will most likely use more memory than storing the pre-calculated vectors.

        // TODO: maybe try lazy create/cache of the asv's in the threads
        
        ivy_mike::timer t1;



        assert( m_ref_aux.empty() && m_ref_pvecs.empty() );

        m_ref_pvecs.reserve( m_ec.m_edges.size() );
        m_ref_aux.reserve( m_ec.m_edges.size() );


        for( size_t i = 0; i < m_ec.m_edges.size(); i++ ) {
            pvec_t root_pvec;

//             std::cout << "newview for branch " << i << ": " << *(m_ec.m_edges[i].first->m_data) << " " << *(m_ec.m_edges[i].second->m_data) << "\n";

            if( i == 340 ) {
                g_dump_aux = true;
            }

            do_newview( root_pvec, m_ec.m_edges[i].first, m_ec.m_edges[i].second, true );

            g_dump_aux = false;
            // TODO: try something fancy with rvalue refs...

            m_ref_pvecs.push_back( std::vector<int>() );
            m_ref_aux.push_back( std::vector<unsigned int>() );

            root_pvec.to_int_vec(m_ref_pvecs.back());
            root_pvec.to_aux_vec(m_ref_aux.back());


            m_ref_gapp.push_back( std::vector<double>() );

            if( ivy_mike::same_type<pvec_t,pvec_pgap>::result ) {
            	// WTF: this is why mixing static and dynamic polymorphism is a BAD idea!
            	pvec_pgap *rvp = reinterpret_cast<pvec_pgap *>(&root_pvec);
            	rvp->to_gap_post_vec(m_ref_gapp.back());

//            	std::transform( m_ref_gapp.back().begin(), m_ref_gapp.back().end(), std::ostream_iterator<int>(std::cout), ivy_mike::scaler_clamp<double>(10,0,9) );
//
//            	std::cout << "\n";
            }


        }
        
        std::cout << "pvecs created: " << t1.elapsed() << "\n";
        
    }
    void seq_to_position_map(const std::vector< uint8_t >& seq, std::vector< int > &map) {
        for( size_t i = 0; i < seq.size(); ++i ) {
            uint8_t ps = dna_parsimony_mapping::d2p(seq[i]);

            if( ps == 0x1 || ps == 0x2 || ps == 0x4 || ps == 0x8 ) {
                map.push_back(int(i));
            }
        }
    }
    
    void gapstream_to_position_map( const std::vector< uint8_t >& gaps, std::vector< int > &map) { 
        align_utils::trace_to_position_map( gaps, &map );
//        int seq_ptr = 0;
//
//        for ( std::vector<uint8_t>::const_reverse_iterator git = gaps.rbegin(); git != gaps.rend(); ++git ) {
//
//            if ( *git == 1) {
//                ++seq_ptr;
//            } else if ( *git == 0 ) {
//
//                map.push_back(seq_ptr);
//
//                ++seq_ptr;
//            } else {
//                map.push_back(seq_ptr);
//
//            }
//        }
    }
    
    void write_qs_pvecs( const char * name ) {
        std::ofstream os( name );
        
        os << m_qs_pvecs.size();
        for( std::vector< std::vector< uint8_t > >::iterator it = m_qs_pvecs.begin(); it != m_qs_pvecs.end(); ++it ) {
            os << " " << it->size() << " ";
            os.write( (char *)it->data(), it->size() );
            
        }
    }
    void write_ref_pvecs( const char * name ) {
        std::ofstream os( name );
        
        os << m_ref_pvecs.size();
        for( size_t i = 0; i < m_ref_pvecs.size(); ++i ) {
            os << " " << m_ref_pvecs[i].size() << " ";
            os.write( (char *)m_ref_pvecs[i].data(), m_ref_pvecs[i].size() * sizeof(int));
            os.write( (char *)m_ref_aux[i].data(), m_ref_aux[i].size() * sizeof(unsigned int));
        }
    }


public:



    papara_nt( const char* opt_tree_name, const char *opt_alignment_name, const char *opt_qs_name, bool write_testbench )
      : m_ln_pool(new ln_pool( std::auto_ptr<node_data_factory>(new my_fact<my_adata>) ))
    {

            //std::cerr << "papara_nt instantiated as: " << typeid(*this).name() << "\n";
        lout << "papara_nt instantiated as: " << ivy_mike::demangle(typeid(*this).name()) << "\n";
        lout << "scores: " << score_gap_open << " " << score_gap_extend << " " << score_mismatch << " " << score_match_cgap << "\n";
        
        
        
        std::cerr << ivy_mike::isa<papara_nt<pvec_cgap> >(*this) << " " << ivy_mike::isa<papara_nt<pvec_pgap> >(*this) << "\n";
        // load input data: ref-tree, ref-alignment and query sequences

        //
        // parse the reference tree
        //


        ln_pool &pool = *m_ln_pool;
        tree_parser_ms::parser tp( opt_tree_name, pool );
        tree_parser_ms::lnode * n = tp.parse();

        n = towards_tree( n );
        //
        // create map from tip names to tip nodes
        //
        typedef tip_collector<lnode> tc_t;
        tc_t tc;

        visit_lnode( n, tc );

        std::map<std::string, std::shared_ptr<lnode> > name_to_lnode;

        for( std::vector< std::shared_ptr<lnode> >::iterator it = tc.m_nodes.begin(); it != tc.m_nodes.end(); ++it ) {
//             std::cout << (*it)->m_data->tipName << "\n";
            name_to_lnode[(*it)->m_data->tipName] = *it;
        }


        {
            //
            // read reference alignment: store the ref-seqs in the tips of the ref-tree
            //
            multiple_alignment ref_ma;
            ref_ma.load_phylip( opt_alignment_name );

            

            for( unsigned int i = 0; i < ref_ma.names.size(); i++ ) {

                std::map< std::string, std::shared_ptr<lnode> >::iterator it = name_to_lnode.find(ref_ma.names[i]);

                // process sequences fomr the ref_ma depending on, if they are contained in the tree.
                // if they are, they are 'swapped' into m_ref_seqs
                // if they are not, into m_qs_seqs. (gaps in the QS are removed later)
                
                if( it != name_to_lnode.end() ) {
                    std::shared_ptr< lnode > ln = it->second;
                    //      adata *ad = ln->m_data.get();

                    assert( ivy_mike::isa<my_adata>(ln->m_data.get()) ); //typeid(*ln->m_data.get()) == typeid(my_adata ) );
                    my_adata *adata = static_cast<my_adata *> (ln->m_data.get());

                    m_ref_names.push_back(std::string() );
                    m_ref_seqs.push_back(std::vector<uint8_t>() );

                    m_ref_names.back().swap( ref_ma.names[i] );
                    m_ref_seqs.back().swap( ref_ma.data[i] );

                    // WARNING: make sure not to keep references to elements of m_ref_seqs at this point!
                    adata->init_pvec( m_ref_seqs.back() );
                } else {
                    
                    //std::cout << "transfer ref -> qs: " << ref_ma.names[i] << "\n";
                    m_qs_names.push_back(std::string() );
                    m_qs_seqs.push_back(std::vector<uint8_t>() );

                    m_qs_names.back().swap( ref_ma.names[i] );
                    m_qs_seqs.back().swap( ref_ma.data[i] );
                }
            }
        }
        probgap_model pm( m_ref_seqs );

        std::cout << "p: " << pm.setup_pmatrix(0.1) << "\n";

        stupid_ptr_guard<probgap_model> spg( pvec_pgap::pgap_model, &pm );

        //
        // collect list of edges
        //

        visit_edges( n, m_ec );

        lout << "edges: " << m_ec.m_edges.size() << "\n";

//         std::vector< pvec_t > m_parsvecs;
//         m_parsvecs.resize( m_ec.m_edges.size() );


        //
        // read query sequences
        //

        if( opt_qs_name != 0 ) {
            std::ifstream qsf( opt_qs_name );

            if( qsf.bad() ) {
            	throw std::runtime_error( "cannot open qs file");
            }

            // mix them with the qs from the ref alignment
            read_fasta( qsf, m_qs_names, m_qs_seqs);
        }

        if( m_qs_names.empty() ) {
            throw std::runtime_error( "no qs" );
        }

        //
        // setup qs best-score/best-edge lists
        //


        m_qs_pvecs.resize( m_qs_names.size() );

        m_qs_bestscore.resize(m_qs_names.size());
        std::fill( m_qs_bestscore.begin(), m_qs_bestscore.end(), 32000);
        m_qs_bestedge.resize(m_qs_names.size());

        //
        // pre-create the reference pvecs/auxvecs
        //

        build_ref_vecs();


        //
        // preprocess query sequences
        //

        for( size_t i = 0; i < m_qs_seqs.size(); i++ ) {
            seq_to_nongappy_pvec( m_qs_seqs[i], m_qs_pvecs[i] );
        }
        
        if( write_testbench ) {

        	write_qs_pvecs( "qs.bin" );
        	write_ref_pvecs( "ref.bin" );
        }

    }



    ivy_mike::mutex m_score_mtx;

    struct tb_thread_entry {
    	papara_nt *pnt;
    	size_t rank;
    	size_t num;

    	tb_thread_entry (papara_nt *pnt_, size_t rank_, size_t num_ ) : pnt(pnt_), rank(rank_), num(num_) {}

    	void operator()() {
    		pnt->calc_scores_testbench(rank, num);
    	}
    };




    void calc_scores_testbench( size_t rank = 0, size_t num_threads = 0 ) {


    	pars_align_gapp_seq::arrays arrays;


    	std::vector<scoring_result>results;

    	const size_t num_edges = m_ref_pvecs.size();

    	for( size_t j = 0; j < num_edges; ++j ) {

    		if( num_threads != 0 && (j % num_threads) != rank ) {
    			continue;
    		}

//    		std::cout << "thread " << rank << " " << j << "\n";

    		int *seqptr = m_ref_pvecs[j].data();
    		double *gappptr = m_ref_gapp[j].data();

    		assert( !m_ref_gapp[j].empty());


    		size_t ref_len = m_ref_pvecs[j].size();

			for( unsigned int i = 0; i < m_qs_names.size(); i++ ) {



				size_t stride = 1;
				size_t aux_stride = 1;


				pars_align_gapp_seq pas( seqptr, m_qs_pvecs[i].data(), ref_len, m_qs_pvecs[i].size(), stride, gappptr, aux_stride, arrays, 0, score_gap_open, score_gap_extend, score_mismatch, score_match_cgap );
				int res = pas.alignFreeshift(INT_MAX);
				results.push_back(scoring_result(j, i, res));

			}
    	}

    	{
			ivy_mike::lock_guard<ivy_mike::mutex> lock(m_score_mtx);

			for( std::vector<scoring_result>::iterator it = results.begin(); it != results.end(); ++it ) {
				if( it->res < m_qs_bestscore[it->qs] || (it->res == m_qs_bestscore[it->qs] && int(it->edge) < m_qs_bestedge[it->qs] )) {

					m_qs_bestscore[it->qs] = it->res;

					assert( it->edge <= size_t(std::numeric_limits<int>::max()) );
					m_qs_bestedge[it->qs] = int(it->edge); // TODO: review: can m_qs_bestedge become size_t? 
				}

			}



		}


    }

    void calc_scores( size_t n_threads ) {


    	if( true ) {
			//
			// build the alignment blocks
			//


			build_block_queue();

			//
			// work
			//
			ivy_mike::timer t1;
			ivy_mike::thread_group tg;
			lout << "start scoring, using " << n_threads <<  " threads\n";


			while( tg.size() < n_threads ) {
				tg.create_thread(worker(*this, tg.size()));
			}
			tg.join_all();

			lout << "scoring finished: " << t1.elapsed() << "\n";
    	} else {
    		lout << "testbench mode\n";

    		if( !false ) {
				ivy_mike::thread_group tg;


				for( size_t i = 0; i < n_threads; ++i ) {
					tg.create_thread( tb_thread_entry( this, i, n_threads));
				}

				tg.join_all();
				std::cout << "joined threads\n";
    		} else {
    			calc_scores_testbench();
    		}

    	}
    }

    void print_best_scores( std::ostream &os ) {
        boost::io::ios_all_saver ioss(os);
        os << std::setfill ('0');
        for( unsigned int i = 0; i < m_qs_names.size(); i++ ) {
            os << m_qs_names[i] << " "  << std::setw (4) << m_qs_bestedge[i] << " " << std::setw(5) << m_qs_bestscore[i] << "\n";

        }
    }


    void gapstream_to_alignment( const std::vector<uint8_t> &gaps, const std::vector<uint8_t> &raw, std::vector<uint8_t> &out, uint8_t gap_char ) {

        std::vector<uint8_t>::const_reverse_iterator rit = raw.rbegin();

        for ( std::vector<uint8_t>::const_iterator git = gaps.begin(); git != gaps.end(); ++git ) {

            if ( *git == 1) {
                out.push_back(gap_char);
            } else if ( *git == 0 ) {
            	assert( rit < raw.rend() );
                out.push_back(*rit);
                ++rit;
            } else {
                ++rit; // just consume one QS character
            }
        }

        std::reverse(out.begin(), out.end());
    }


    static uint8_t normalize_dna( uint8_t c ) {
        c = std::toupper(c);

        if( c == 'U' ) {
            c = 'T';
        }

        return c;
    }

    static char num_to_ascii( int n ) {
        if( n >= 0 && n <= 9 ) {
            return '0' + n;
        } else if( n >= 0xa && n <= 0xf ) {
            return 'a' + n;
        } else {
            throw std::runtime_error( "not a single digit (hex) number" );
        }
    }

    void align_best_scores( std::ostream &os, std::ostream &os_quality ) {
        // create the actual alignments for the best scoring insertion position (=do the traceback)
        
        lout << "generating best scoring alignments\n";
        ivy_mike::timer t1;
        pars_align_seq<>::arrays seq_arrays(true);
        pars_align_gapp_seq::arrays seq_arrays_gapp(true);

        double mean_quality = 0.0;
        double n_quality = 0.0;

        const size_t pad = max_name_len() + 1;

        for( unsigned int i = 0; i < m_qs_names.size(); i++ ) {
            int best_edge = m_qs_bestedge[i];

            assert( best_edge >= 0 && size_t(best_edge) < m_ref_pvecs.size() );

            int res = -1;
            std::vector<uint8_t> tbv;

            if( true ) {
				const int *seqptr = m_ref_pvecs[best_edge].data();
				const unsigned int *auxptr = m_ref_aux[best_edge].data();

				const size_t ref_len = m_ref_pvecs[best_edge].size();

				const size_t stride = 1;
				const size_t aux_stride = 1;
				pars_align_seq<> pas( seqptr, m_qs_pvecs[i].data(), ref_len, m_qs_pvecs[i].size(), stride, auxptr, aux_stride, seq_arrays, 0, score_gap_open, score_gap_extend, score_mismatch, score_match_cgap );
				res = pas.alignFreeshift(INT_MAX);
				pas.tracebackCompressed(tbv);
            } else {
            	const int *seqptr = m_ref_pvecs[best_edge].data();

            	assert( !m_ref_gapp[best_edge].empty() );

				const double *gappptr = m_ref_gapp[best_edge].data();

				const size_t ref_len = m_ref_pvecs[best_edge].size();

				const size_t stride = 1;
				const size_t aux_stride = 1;
				pars_align_gapp_seq pas( seqptr, m_qs_pvecs[i].data(), ref_len, m_qs_pvecs[i].size(), stride, gappptr, aux_stride, seq_arrays_gapp, 0, score_gap_open, score_gap_extend, score_mismatch, score_match_cgap );
				res = pas.alignFreeshift(INT_MAX);
				pas.tracebackCompressed(tbv);
            }
//            std::copy( tbv.begin(), tbv.end(), std::ostream_iterator<int>(std::cout, " " ));
//            std::cout << "\n";



            std::vector<uint8_t> out_qs;


            gapstream_to_alignment(tbv, m_qs_pvecs[i], out_qs, 0xf);

            //
            // in place transform out_qs from pstate to dna form. WATCH OUT!
            //

            std::transform( out_qs.begin(), out_qs.end(), out_qs.begin(), dna_parsimony_mapping::p2d );

            os << std::setw(pad) << std::left << m_qs_names[i];
            std::copy( out_qs.begin(), out_qs.end(), std::ostream_iterator<char>(os));
            os << "\n";

            const bool dump_auxvec = false;
            if( dump_auxvec )
            {
                std::string auxv;
                auxv.resize(m_ref_aux[best_edge].size());

                std::transform( m_ref_aux[best_edge].begin(), m_ref_aux[best_edge].end(), auxv.begin(), num_to_ascii );
                os << m_qs_names[i] << "\t" << auxv << "\n";
            }


            if( res != m_qs_bestscore[i] ) {
                std::cout << "meeeeeeep! score: " << m_qs_bestscore[i] << " " << res << "\n";
            }

            if( os_quality.good() && out_qs.size() == m_qs_seqs[i].size() ) {
                

                std::vector<int> map_ref;
                std::vector<int> map_aligned;
                seq_to_position_map( m_qs_seqs[i], map_ref );
                gapstream_to_position_map( tbv, map_aligned );
                

                if( map_ref.size() != map_aligned.size() ) {
                    throw std::runtime_error( "alignment quirk: map_ref.size() != map_aligned.size()" );
                }
                
                size_t num_equal = ivy_mike::count_equal( map_ref.begin(), map_ref.end(), map_aligned.begin() );
                
                //std::cout << "size: " << map_ref.size() << " " << map_aligned.size() << " " << m_qs_seqs[i].size() << "\n";
                //std::cout << num_equal << " equal of " << map_ref.size() << "\n";
                
                double score = num_equal / double(map_ref.size());
                //double score = alignment_quality( out_qs, m_qs_seqs[i], debug );

                os_quality << m_qs_names[i] << " " << score << "\n";

                mean_quality += score;
                n_quality += 1;
            }


        }
        lout << "alignment finished: " << t1.elapsed() << "\n"; 
        lout << "mean quality: " << mean_quality / n_quality << "\n";

    }

    ~papara_nt() {

        ivy_mike::timer t2;
        m_ln_pool->clear();
        //   pool.mark(n);
        m_ln_pool->sweep();

        std::cout << t2.elapsed() << std::endl;

    }

    size_t max_name_len() {
        size_t ml = 0;

        for( size_t i = 0; i < m_ref_names.size(); i++ ) {
            ml = std::max( ml, m_ref_names[i].size() );
        }


        for( size_t i = 0; i < m_qs_names.size(); i++ ) {
            ml = std::max( ml, m_qs_names[i].size() );
        }

        return ml;
    }

    void dump_ref_seqs ( std::ostream &os ) {
        const size_t pad = max_name_len() + 1;

        for( size_t i = 0; i < m_ref_seqs.size(); i++ ) {
            std::string outs;
            outs.resize(m_ref_seqs[i].size() );

            std::transform( m_ref_seqs[i].begin(), m_ref_seqs[i].end(), outs.begin(), normalize_dna );

            os << std::setw(pad) << std::left << m_ref_names[i] << outs << "\n";
        }
    }


    void write_result_phylip( std::ostream &os, std::ostream &os_quality ) {
        os << " " << m_ref_seqs.size() + m_qs_names.size() << " " << m_ref_seqs.at(0).size() << "\n";
        dump_ref_seqs(os);
        align_best_scores(os, os_quality);

    }
    double alignment_quality_very_strict ( const std::vector< uint8_t > &s1, const std::vector< uint8_t >& s2, bool debug = false ) {
        size_t nident = 0;
        size_t ngap1 = 0;
        size_t ngap2 = 0;


        for( std::vector< uint8_t >::const_iterator it1 = s1.begin(), it2 = s2.begin(); it1 != s1.end(); ++it1, ++it2 ) {

            if( dna_parsimony_mapping::is_gap( *it1 ) ) {
                ngap1++;
            }

            if( dna_parsimony_mapping::is_gap( *it2 ) ) {
                ngap2++;
            }
            if( debug ) {
                std::cerr << ngap1 << " " << ngap2 << " " << *it1 << " " << *it2 << "\n";
            }

            if( ngap1 == ngap2 ) {
                nident++;
            }
        }

        return double(nident) / s1.size();

    }
    double alignment_quality ( const std::vector< uint8_t > &s1, const std::vector< uint8_t >& s2, bool debug = false ) {
        size_t nident = 0;

//         size_t nident_nongap = 0;
//         size_t n_nongap = 0;

        for( std::vector< uint8_t >::const_iterator it1 = s1.begin(), it2 = s2.begin(); it1 != s1.end(); ++it1, ++it2 ) {
            if( dna_parsimony_mapping::d2p(*it1) == dna_parsimony_mapping::d2p(*it2) ) {
                nident++;
            }
        }

        return double(nident) / s1.size();

    }

};

std::string filename( const std::string &run_name, const char *type ) {
    std::stringstream ss;
    
    ss << "papara_" << type << "." << run_name;
    
    return ss.str();
}

bool file_exists(const char *filename)
{
  std::ifstream is(filename);
  return is.good();
}


int main( int argc, char *argv[] ) {

//     aligned_buffer<int> xxx(1024);
    
    

    
    namespace igo = ivy_mike::getopt;

    ivy_mike::getopt::parser igp;

    std::string opt_tree_name;
    std::string opt_alignment_name;
    std::string opt_qs_name;
    bool opt_use_cgap;
    int opt_num_threads;
    std::string opt_run_name;
    bool opt_write_testbench;

    
    igp.add_opt( 't', igo::value<std::string>(opt_tree_name) );
    igp.add_opt( 's', igo::value<std::string>(opt_alignment_name) );
    igp.add_opt( 'q', igo::value<std::string>(opt_qs_name) );
    igp.add_opt( 'c', igo::value<bool>(opt_use_cgap, true).set_default(false) );
    igp.add_opt( 'j', igo::value<int>(opt_num_threads).set_default(1) );
    igp.add_opt( 'n', igo::value<std::string>(opt_run_name).set_default("default") );
    igp.add_opt( 'b', igo::value<bool>(opt_write_testbench, true).set_default(false) );
    
    igp.parse(argc,argv);

    if( igp.opt_count('t') != 1 || igp.opt_count('s') != 1  ) {
        std::cerr << "missing options -t and/or -s (-q is optional)\n";
        return 0;
    }
    ivy_mike::timer t;

    const char *qs_name = 0;
    if( !opt_qs_name.empty() ) {
        qs_name = opt_qs_name.c_str();
    }
    
    std::string log_filename = filename( opt_run_name, "log" );

    if( opt_run_name != "default" && file_exists(log_filename.c_str()) ) {
        std::cout << "log file already exists for run '" << opt_run_name << "'\n";
        return 0;
    }
    
    std::ofstream logs( log_filename.c_str());
    if( !logs ) {
        std::cout << "could not open logfile for writing: " << log_filename << std::endl;
        return 0;
    }
    
    log_device ldev( std::cout, logs );
    log_stream_guard lout_guard( lout, ldev );
    
    
    
    std::auto_ptr<papara_nt_i> pnt_ptr;

    if( !opt_use_cgap ) {
        pnt_ptr.reset( new papara_nt<pvec_pgap>( opt_tree_name.c_str(), opt_alignment_name.c_str(), qs_name, opt_write_testbench ));
    } else {
        pnt_ptr.reset( new papara_nt<pvec_cgap>( opt_tree_name.c_str(), opt_alignment_name.c_str(), qs_name, opt_write_testbench ));
    }
    
//     return 0;
    
    std::cerr << "using " << opt_num_threads << " threads\n";


    papara_nt_i &pnt = *pnt_ptr;
    pnt.calc_scores( opt_num_threads );

    {
        std::ofstream os( filename( opt_run_name, "scores" ).c_str() );
        pnt.print_best_scores(os);
    }

    {
        std::ofstream os( filename( opt_run_name, "alignment" ).c_str() );
        std::ofstream os_quality( filename( opt_run_name, "quality" ).c_str() );

        //         pnt.dump_ref_seqs(os);
        //         pnt.align_best_scores(os);
        pnt.write_result_phylip(os, os_quality);
    }

    //ivymike::LN *n = tp.parse();

//     getchar();
    //ivymike::LN::free( n );
//     delete n;


    std::cout << t.elapsed() << std::endl;
    lout << "SUCCESS " << t.elapsed() << std::endl;
    return 0;
//     getchar();
}

