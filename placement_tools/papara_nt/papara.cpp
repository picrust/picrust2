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

#include <stdint.h>
#include <iomanip>
#include <boost/bind.hpp>
#include <boost/dynamic_bitset.hpp>
#include <iterator>

#include "ivymike/fasta.h"
#include "ivymike/demangle.h"
#include "ivymike/time.h"

#include "papara.h"
#include "vec_unit.h"
#include "align_pvec_vec.h"
#include "stepwise_align.h"
#include "align_utils.h"






using namespace ivy_mike;
using namespace ivy_mike::tree_parser_ms;

using namespace papara;

bool papara::g_dump_aux = false;


class log_stream_buffer : public std::streambuf
{

public:
    
    log_stream_buffer() : buffer_(120) {
        
    }

    void post( char overflow, char *start, char *end ) {
        for( std::vector< std::ostream* >::iterator it = log_tees.begin(); it != log_tees.end(); ++it ) {
            std::copy( start, end, std::ostream_iterator<char>( *(*it) ));
            
            if( overflow != 0 ) {
                (*it)->put(overflow);
            }
        }
        
        
        for( std::vector< log_sink* >::iterator it = log_sinks.begin(); it != log_sinks.end(); ++it ) {
            (*it)->post( overflow, start, end );
        }
    }
    
    int_type overflow(int c) {


        post( c, pbase(), epptr() );

        setp(&buffer_.front(), (&buffer_.back()) + 1);

        return 1;
    }

    int sync() {
        //            std::cout << "sync:";
        //            std::copy( pbase(), pptr(), std::ostream_iterator<char>(std::cout));
        //            std::cout << std::endl;
        
        post( 0, pbase(), pptr() );
        setp(&buffer_.front(), (&buffer_.back()) + 1);
        return 0;
    }
    
    void add_tee( std::ostream *os ) {
        log_tees.push_back(os);
    }
    void remove_tee( std::ostream *os ) {
        std::vector<std::ostream *>::iterator it = std::find( log_tees.begin(), log_tees.end(), os );
        if( it != log_tees.end() ) {
            log_tees.erase(it);
        }
    }
    
    void add_sink( log_sink *s ) {
        log_sinks.push_back(s);
    }
    void remove_sink( log_sink *s ) {
        std::vector<log_sink *>::iterator it = std::find( log_sinks.begin(), log_sinks.end(), s );
        if( it != log_sinks.end() ) {
            log_sinks.erase(it);
        }
    }
    
private:
   // copy ctor and assignment not implemented;
    // copying not allowed
    log_stream_buffer(const log_stream_buffer &);
    log_stream_buffer &operator= (const log_stream_buffer &);

    std::vector<char> buffer_;
    std::vector<std::ostream *> log_tees;
    std::vector<log_sink *> log_sinks;
};

// // lifted from http://www.horstmann.com/cpp/iostreams.html
// class log_stream_buffer : public std::filebuf
// {
//     std::vector<std::ostream *> log_tees;
//     
// public:
//     log_stream_buffer() { 
//         //filebuf::open("NUL", ios::out); 
//     }
//     
//     void open(const char *fname);
// //     void close() { 
// // //         _win.close(); filebuf::close();
// //         
// //     } 
//     virtual int sync();
//     
//     void add_tee( std::ostream *os ) {
//         log_tees.push_back(os);
//     }
//     void remove_tee( std::ostream *os ) {
//         std::vector<std::ostream *>::iterator it = std::find( log_tees.begin(), log_tees.end(), os );
//         if( it != log_tees.end() ) {
//             log_tees.erase(it);
//         }
//     }
// private:
//     
// };
// void log_stream_buffer::open(const char *fname)
// {  
//     std::filebuf::close();
//     if( fname != 0 ) {
//         std::cerr << "log file: " << fname << "\n";
//         std::filebuf::open(fname, std::ios::out | std::ios::app | std::ios::trunc );
//         assert(std::filebuf::is_open());
//     }
// 
// } 
// int log_stream_buffer::sync()
// {  
//     //     pbase();
//     //     pptr()
//     //     
//     //     int count = std::filebuf::out_waiting();
//     // _win.append(pbase(), count);
//     std::cerr << "sync: " << std::distance(pbase(), pptr()) << "\n";
//     for( std::vector< std::ostream* >::iterator it = log_tees.begin(); it != log_tees.end(); ++it ) {
//         std::copy( pbase(), pptr(), std::ostream_iterator<char>( *(*it) ));
//     }
//     return std::filebuf::sync(); 
// }
      

static log_stream_buffer ls_buf;
std::ostream papara::lout(&ls_buf);

static ivy_mike::mutex log_buffer_mutex;


// open_log_file::open_log_file( const char *filename ) {
//     
//     ivy_mike::lock_guard<ivy_mike::mutex> lock( log_buffer_mutex );
//     if( log_is_open ) {
//         std::cerr << "papara logfile is already open." << std::endl;
//         abort();
//     }
//     
//     
//     ls_buf.open(filename);
//     log_is_open = true;
// }

// open_log_file::~open_log_file() {
//     ivy_mike::lock_guard<ivy_mike::mutex> lock( log_buffer_mutex );
//     if( !log_is_open ) {
//         std::cerr << "papara logfile is already closed. (This should be impossible) " << std::endl;
//         abort();
//     }
//     ls_buf.close();
//     log_is_open = false;
// }


add_log_tee::add_log_tee( std::ostream &os ) : os_(os) {
    ivy_mike::lock_guard<ivy_mike::mutex> lock( log_buffer_mutex );
    ls_buf.add_tee(&os);
}

add_log_tee::~add_log_tee() {
    ivy_mike::lock_guard<ivy_mike::mutex> lock( log_buffer_mutex );
    ls_buf.remove_tee(&os_);
}


add_log_sink::add_log_sink( log_sink *s ) : s_(s) {
    ivy_mike::lock_guard<ivy_mike::mutex> lock( log_buffer_mutex );
    ls_buf.add_sink(s);
}

add_log_sink::~add_log_sink() {
    ivy_mike::lock_guard<ivy_mike::mutex> lock( log_buffer_mutex );
    ls_buf.remove_sink(s_);
}


std::string papara::get_version_string() {
    return std::string( "2.5" );
}

template<typename seq_tag>
queries<seq_tag>::queries( const std::string &opt_qs_name ) {


//        if( !opt_qs_name.empty() ) {
    //
    // read query sequences
    //
    
    if( !opt_qs_name.empty() ) {
        std::ifstream qsf( opt_qs_name.c_str() );
        
        if( !qsf.good() ) {
            throw std::runtime_error( "cannot open qs file");
        }
        
        // mix them with the qs from the ref alignment <-- WTF? must have been sniffing whiteboard cleaner... the qs are read before the ref seqs...
        read_fasta( qsf, m_qs_names, m_qs_seqs);
    }
    
    //            if( m_qs_names.empty() ) {
        //                throw std::runtime_error( "no qs" );
        //            }
        
        std::for_each( m_qs_names.begin(), m_qs_names.end(), std::ptr_fun( normalize_name ));
        
        //
        // setup qs best-score/best-edge lists
        //
        
        
        m_qs_pvecs.resize( m_qs_names.size() );
        //        }
        //        m_qs_bestscore.resize(m_qs_names.size());
        //        std::fill( m_qs_bestscore.begin(), m_qs_bestscore.end(), 32000);
        //        m_qs_bestedge.resize(m_qs_names.size());
        
        
        
        
}
template<typename seq_tag>
void queries<seq_tag>::preprocess() {
    //
    // preprocess query sequences
    //
    if( m_qs_seqs.empty() ) {
        throw std::runtime_error( "no query sequences" );
    }

    assert( m_qs_seqs.size() == m_qs_names.size() );
    m_qs_pvecs.resize(m_qs_seqs.size());
    m_qs_cseqs.resize(m_qs_seqs.size());

    std::vector<bool> bad_characters( 256, false );
    
    
    
    for( size_t i = 0; i < m_qs_seqs.size(); i++ ) {
//            seq_to_nongappy_pvec( m_qs_seqs[i], m_qs_pvecs[i] );
        //          static void seq_to_nongappy_pvec( std::vector<uint8_t> &seq, std::vector<uint8_t> &pvec ) {

        // the following line means: transform sequence to pvec using seq_model::s2p as mapping
        // function and only append the mapped character into pvec, if it corresponds to a single (=non gap)
        // character.
        
        std::vector<uint8_t> qs_tmp;
        qs_tmp.reserve( m_qs_seqs[i].size() );
        bool bad_char = false;
        for( std::vector<uint8_t>::iterator it = m_qs_seqs[i].begin(), e = m_qs_seqs[i].end(); it != e; ++it ) {
            if( !seq_model::is_known_sstate( *it ) ) {
                bad_characters.at( *it ) = true;
                bad_char = true;
            } else {
                qs_tmp.push_back( *it );
            }
        }
        if( bad_char ) {
            m_qs_seqs[i].swap( qs_tmp );
        }

        std::transform( m_qs_seqs[i].begin(), m_qs_seqs[i].end(),
                        back_insert_ifer(m_qs_cseqs[i], seq_model::cstate_is_single),
                        seq_model::s2c );

        std::transform( m_qs_seqs[i].begin(), m_qs_seqs[i].end(),
                        back_insert_ifer(m_qs_pvecs[i], seq_model::pstate_is_single),
                        seq_model::s2p );

//         std::cout << "preprocess: " << i << " " << m_qs_cseqs[i].size() << " " << m_qs_pvecs[i].size() << "\n";
        
        if( m_qs_cseqs[i].size() != m_qs_pvecs[i].size() ) {
            // check for quirks related to the p-state vs c-state representation
            throw std::runtime_error( "mismatch between lengths of c-state and p-state representations of query sequence." );
        }
        
//            for( unsigned int i = 0; i < seq.size(); i++ ) {
//                seq_model::pars_state_t ps = seq_model::s2p(seq[i]);
//
//                if( seq_model::is_single(ps)) {
//                    pvec.push_back(ps);
//                }
//
//            }



    }

    // print warnings about deleted characters
    bool warn_header = false;
    for( std::vector<bool>::iterator it = bad_characters.begin(), e = bad_characters.end(); it != e; ++it ) {
        if( *it ) {
            if( !warn_header ) {
                std::cout << "WARNING: there were unsupported characters in the query sequences. They will be deleted:\n";
                warn_header = true;
            }
            
            std::cout << "deleted character: '" << uint8_t(std::distance( bad_characters.begin(), it )) << "'\n";
        }
    }
    
//        if( write_testbench ) {
//
//            write_qs_pvecs( "qs.bin" );
//            write_ref_pvecs( "ref.bin" );
//        }

}


// template<typename pvec_t, typename seq_tag>
// void queries<seq_tag>::init_partition_assignments( partassign::part_assignment &part_assign, references<pvec_t,seq_tag> &refs ) {
//     
// }



template<typename seq_tag>
void queries<seq_tag>::add( const std::string& name, std::vector< uint8_t >& qs ) {
    m_qs_names.push_back(name);
    ivy_mike::push_back_swap(m_qs_seqs, qs );
}
template<typename seq_tag>
void queries<seq_tag>::write_pvecs(const char* name) {
    std::ofstream os( name );

    os << m_qs_pvecs.size();
    for( typename std::vector< std::vector< pars_state_t > >::iterator it = m_qs_pvecs.begin(); it != m_qs_pvecs.end(); ++it ) {
        os << " " << it->size() << " ";
        os.write( (char *)it->data(), it->size() );

    }
}
template<typename seq_tag>
size_t queries<seq_tag>::max_name_length() const {
    size_t len = 0;
    for( std::vector <std::string >::const_iterator it = m_qs_names.begin(); it != m_qs_names.end(); ++it ) {
        len = std::max( len, it->size() );
    }

    return len;
}

template<typename seq_tag>
size_t queries<seq_tag>::calc_cups_per_ref(size_t ref_len) const {
    size_t ct = 0;

    typename std::vector<std::vector <pars_state_t> >::const_iterator first = m_qs_pvecs.begin();
    const typename std::vector<std::vector <pars_state_t> >::const_iterator last = m_qs_pvecs.end();

    for(; first != last; ++first ) {
        //ct += (ref_len - first->size()) * first->size();
        ct += ref_len * first->size(); // papara now uses the 'unbanded' aligner
    }



    return ct;
}
template<typename seq_tag>
void queries<seq_tag>::normalize_name(std::string& str) {
    std::string ns;

    // replace runs of one or more whitespaces with an underscores
    bool in_ws_run = false;

    for( std::string::iterator it = str.begin(); it != str.end(); ++it ) {

        if( std::isspace(*it) ) {
            if( !in_ws_run ) {

                if( (*it) != '\r' ) { // just ignore stupid windows line ends
                    ns.push_back( '_' );
                }
                in_ws_run = true;
            }
        } else {
            ns.push_back(*it);
            in_ws_run = false;
        }

    }
    str.swap(ns);

}

//////////////////////////////////////////////////////////////
// references stuff
//////////////////////////////////////////////////////////////

template<typename pvec_t, typename seq_tag>
references<pvec_t,seq_tag>::references(const char* opt_tree_name, const char* opt_alignment_name, queries<seq_tag>* qs) : m_ln_pool(new ln_pool( std::unique_ptr<node_data_factory>(new my_fact<my_adata>) )), spg_(pvec_pgap::pgap_model, &pm_)
{

    //std::cerr << "papara_nt instantiated as: " << typeid(*this).name() << "\n";
    lout << "references container instantiated as: " << ivy_mike::demangle(typeid(*this).name()) << "\n";




//        std::cerr << ivy_mike::isa<papara_nt<pvec_cgap> >(*this) << " " << ivy_mike::isa<papara_nt<pvec_pgap> >(*this) << "\n";
    // load input data: ref-tree, ref-alignment and query sequences

    //
    // parse the reference tree
    //


    ln_pool &pool = *m_ln_pool;
    tree_parser_ms::parser tp( opt_tree_name, pool );
    tree_parser_ms::lnode * n = tp.parse();

    
    
    n = towards_tree( n );
    
    tree_ = std::shared_ptr<im_tree_parser::lnode>(n->get_smart_ptr());
    
    //
    // create map from tip names to tip nodes
    //
    typedef tip_collector<lnode> tc_t;
    tc_t tc;

    visit_lnode( n, tc );

    //boost::dynamic_bitset<> found_tree_taxa( tc.m_nodes.size(), true );



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

        std::vector<my_adata *> tmp_adata;
        boost::dynamic_bitset<> unmasked;

        for( unsigned int i = 0; i < ref_ma.names.size(); i++ ) {

            std::map< std::string, std::shared_ptr<lnode> >::iterator it = name_to_lnode.find(ref_ma.names[i]);

            // process sequences from the ref_ma depending on, if they are contained in the tree.
            // if they are, they are 'swapped' into m_ref_seqs
            // if they are not, into m_qs_seqs. (gaps in the QS are removed later)
            //
            // additionally, all columns that contain only gaps are removed from the reference sequences.

            if( it != name_to_lnode.end() ) {


                std::shared_ptr< lnode > ln = it->second;
                //      adata *ad = ln->m_data.get();

                assert( ivy_mike::isa<my_adata>(ln->m_data.get()) ); //typeid(*ln->m_data.get()) == typeid(my_adata ) );
                my_adata *adata = static_cast<my_adata *> (ln->m_data.get());

                // store the adata ptr corresponding to the current ref sequence for later use.
                // (their indices in m_ref_seqs and tmp_adata correspond.)
                assert( tmp_adata.size() == m_ref_seqs.size() );

                tmp_adata.push_back(adata);
                m_ref_names.push_back(std::string() );
                m_ref_seqs.push_back(std::vector<uint8_t>() );

                m_ref_names.back().swap( ref_ma.names[i] );
                m_ref_seqs.back().swap( ref_ma.data[i] );

                // mark all non-gap positions of the current reference in bit-vector 'unmasked'
                const std::vector<uint8_t> &seq = m_ref_seqs.back();
                if( unmasked.empty() ) {
                    unmasked.resize( seq.size() );
                }
                assert( unmasked.size() == seq.size() );


                for( size_t j = 0, e = seq.size(); j != e; ++j ) {
                    try {
                        unmasked[j] |= !seq_model::pstate_is_gap( seq_model::s2p(seq[j]));

                    } catch( sequence_model::illegal_character x ) {
                        std::cerr << "illegal character in file '" << opt_alignment_name << "'" << std::endl;
                        std::cerr << "row: " << i + 1 << " (name: " << m_ref_names.back() << ")" << std::endl;
                        std::cerr << "col: " << j + 1 << " (char: '" << seq[j] << "')" << std::endl;
                        abort();
                    }
                }

                // erase it from the name to lnode* map, so that it can be used to ideantify tree-taxa without corresponding entries in the alignment
                name_to_lnode.erase(it);

            } else {
                qs->add(ref_ma.names[i], ref_ma.data[i]); // REMARK: the second parameter is 'moved-from' (should be an rvalue-ref)
            }
        }

        if( !name_to_lnode.empty() ) {
            std::cerr << "error: there are " << name_to_lnode.size() << " taxa in the tree with no corresponding sequence in the reference alignment. names:\n";

            for( std::map< std::string, std::shared_ptr< lnode > >::iterator it = name_to_lnode.begin(); it != name_to_lnode.end(); ++it ) {
                std::cout << it->first << "\n";
            }

            throw std::runtime_error( "bailing out due to inconsitent input data\n" );

        }

        {
            // remove all 'pure-gap' columns from the ref sequences

            assert( tmp_adata.size() == m_ref_seqs.size() );

            // retrieve a list of non-gap indices from the bit-vector
            std::vector<size_t> unmasked_idx;
            {
                size_t i = unmasked.find_first();
                while( i != unmasked.npos ) {
//                         std::cout << "um: " << i << "\n";

                    unmasked_idx.push_back(i);
                    i = unmasked.find_next(i);
                }
            }
            for( size_t i = 0, e = m_ref_seqs.size(); i != e; ++i ) {

                std::vector<uint8_t> seq_tmp;
                seq_tmp.reserve(unmasked_idx.size());

                const std::vector<uint8_t> &seq_orig = m_ref_seqs[i];

                // copy all unmasked ref characters to seq_tmp
                for( std::vector<size_t>::iterator it = unmasked_idx.begin(); it != unmasked_idx.end(); ++it ) {
                    assert( *it < seq_orig.size() );


                    seq_tmp.push_back( seq_orig[*it] );
                }
                m_ref_seqs[i].swap( seq_tmp );

                //initialize the corresponding adata object with the cleaned ref seq.
                tmp_adata.at(i)->init_pvec( m_ref_seqs[i] );
            }
        }
    }
    pm_.reset( m_ref_seqs );
    std::cout << "p: " << pm_.setup_pmatrix(0.1) << "\n";

    // initialize empty non-gap map. It is lazily filled as needed when necessary
    ref_ng_map_.resize( m_ref_seqs.size() );
    
    //
    // collect list of edges
    //

    visit_edges( n, m_ec );

    lout << "edges: " << m_ec.m_edges.size() << "\n";

}

template<typename pvec_t, typename seq_tag>
void references<pvec_t,seq_tag>::build_ref_vecs() {
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

        driver<pvec_t,seq_tag>::do_newview( root_pvec, m_ec.m_edges[i].first, m_ec.m_edges[i].second, true );

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

//              std::transform( m_ref_gapp.back().begin(), m_ref_gapp.back().end(), std::ostream_iterator<int>(std::cout), ivy_mike::scaler_clamp<double>(10,0,9) );
//
//              std::cout << "\n";
        }


    }

//     std::cout << "pvecs created: " << t1.elapsed() << "\n";

}
template<typename pvec_t, typename seq_tag>
const std::vector<int> &references<pvec_t,seq_tag>::ng_map_at( size_t i ) {
    std::vector<int> &ng_map = ref_ng_map_.at(i);
    
    if( !ng_map.empty() ) {
        return ng_map;
    }
    
    //std::vector<int> map;
    
    std::vector< uint8_t > &seq = m_ref_seqs.at(i);
    assert( seq.size() < size_t(std::numeric_limits<int>::max()) );
    for( size_t i = 0; i < seq.size(); ++i ) {
        bool is_gap = seq_model::pstate_is_gap( seq_model::s2p(seq[i]));
        
        if( !is_gap ) {
            ng_map.push_back(int(i));
        }
    }
    
    //ng_map.shrink_to_fit();
    std::vector<int>(ng_map).swap(ng_map); // old fashioned shrink_to_fit
    
    return ng_map;
}


template<typename pvec_t, typename seq_tag>
void references<pvec_t,seq_tag>::write_pvecs(const char* name) {
    std::ofstream os( name );

    os << m_ref_pvecs.size();
    for( size_t i = 0; i < m_ref_pvecs.size(); ++i ) {
        os << " " << m_ref_pvecs[i].size() << " ";
        os.write( (char *)m_ref_pvecs[i].data(), m_ref_pvecs[i].size() * sizeof(int));
        os.write( (char *)m_ref_aux[i].data(), m_ref_aux[i].size() * sizeof(unsigned int));
    }
}

template<typename pvec_t, typename seq_tag>
size_t references<pvec_t,seq_tag>::max_name_length() const {
    size_t len = 0;
    for( std::vector <std::string >::const_iterator it = m_ref_names.begin(); it != m_ref_names.end(); ++it ) {
        len = std::max( len, it->size() );
    }

    return len;
}


bool scoring_results::offer(size_t qs, size_t ref, int score) {
    ivy_mike::lock_guard<ivy_mike::mutex> lock(mtx_);

    candss_.at( qs ).offer( score, ref );

    if( best_score_.at(qs) < score || (best_score_.at(qs) == score && ref < best_ref_.at(qs))) {
        best_score_[qs] = score;
        best_ref_.at(qs) = ref;
        return true;
    }

    return false;

}




namespace papara {
void scoring_results::candidates::offer(int score, size_t ref) {

    candidate c( score,ref);

    insert(std::lower_bound( begin(), end(), c), c );

    if( size() > max_num_) {
        pop_back();
    }
}

bool scoring_results::candidate::operator<(const scoring_results::candidate& other) const {
    if( score_ > other.score_ ) {
        return true;
    } else {
        if( score_ == other.score_ ) {
            return ref_ < other.ref_;
        } else {
            return false;
        }
    }
}


////////////////////////////////////////////////////////
// driver stuff
////////////////////////////////////////////////////////

template<typename seq_tag>
class worker {



    const static size_t VW = vu_config<seq_tag>::width;
    typedef typename vu_config<seq_tag>::scalar vu_scalar_t;
    typedef typename block_queue<seq_tag>::block_t block_t;
    typedef model<seq_tag> seq_model;



    block_queue<seq_tag> &block_queue_;
    scoring_results &results_;

    const queries<seq_tag> &qs_;

    const size_t rank_;

    const papara_score_parameters sp_;

    static void copy_to_profile( const block_t &block, aligned_buffer<vu_scalar_t> *prof, aligned_buffer<vu_scalar_t> *aux_prof ) {
        size_t reflen = block.ref_len;



        assert( reflen * VW == prof->size() );
        assert( reflen * VW == aux_prof->size() );

        typename aligned_buffer<vu_scalar_t>::iterator it = prof->begin();
        typename aligned_buffer<vu_scalar_t>::iterator ait = aux_prof->begin();
    //         std::cout << "reflen: " << reflen << " " << size_t(&(*it)) << "\n";

        for( size_t i = 0; i < reflen; ++i ) {
            for( size_t j = 0; j < VW; ++j ) {
    //                 std::cout << "ij: " << i << " " << j << " " << pvecs[j].size() <<  "\n";


                *it = vu_scalar_t(block.seqptrs[j][i]);
                *ait = (block.auxptrs[j][i] == AUX_CGAP) ? vu_scalar_t(-1) : 0;

                ++it;
                ++ait;
            }
        }


        assert( it == prof->end());
        assert( ait == aux_prof->end());

    }

public:
    worker( block_queue<seq_tag> *bq, scoring_results *res, const queries<seq_tag> &qs, size_t rank, const papara_score_parameters &sp ) 
      : block_queue_(*bq), results_(*res), qs_(qs), rank_(rank), sp_(sp) {}
    void operator()() {



        ivy_mike::timer tstatus;
        ivy_mike::timer tprint;

        uint64_t cups_per_ref = -1;


        uint64_t ncup = 0;

        uint64_t inner_iters = 0;
        uint64_t ticks_all = 0;

        uint64_t ncup_short = 0;

        uint64_t inner_iters_short = 0;
        uint64_t ticks_all_short = 0;


        aligned_buffer<vu_scalar_t> pvec_prof;
        aligned_buffer<vu_scalar_t> aux_prof;
        align_vec_arrays<vu_scalar_t> arrays;
        aligned_buffer<vu_scalar_t> out_scores(VW);
        aligned_buffer<vu_scalar_t> out_scores2(VW);
        
        size_t queue_size;
        size_t init_queue_size = -1;
        
        while( true ) {
            block_t block;

            if( !block_queue_.get_block(&block, &queue_size)) {
                break;
            }

            if( init_queue_size == size_t(-1) ) {
                init_queue_size = queue_size;
            }
            
            if( cups_per_ref == uint64_t(-1) ) {
                cups_per_ref = qs_.calc_cups_per_ref(block.ref_len );
            }

#if 1
       //     assert( VW == 8 );

            pvec_prof.resize( VW * block.ref_len );
            aux_prof.resize( VW * block.ref_len );

            copy_to_profile(block, &pvec_prof, &aux_prof );

            pvec_aligner_vec<vu_scalar_t,VW> pav( block.seqptrs, block.auxptrs, block.ref_len, sp_.match, sp_.match_cgap, sp_.gap_open, sp_.gap_extend, seq_model::c2p, seq_model::num_cstates() );

//            const align_pvec_score<vu_scalar_t,VW> aligner( block.seqptrs, block.auxptrs, block.ref_len, score_mismatch, score_match_cgap, score_gap_open, score_gap_extend );
            for( unsigned int i = 0; i < qs_.size(); i++ ) {

                //align_pvec_score_vec<vu_scalar_t, VW, false, typename seq_model::pars_state_t>( pvec_prof, aux_prof, qs_.pvec_at(i), score_match, score_match_cgap, score_gap_open, score_gap_extend, out_scores, arrays );


                std::pair<size_t,size_t> bounds = qs_.get_per_qs_bounds( i );
//		std::cout << "bounds: " << bounds.first << " " << bounds.second << "\n";

                // if no bounds are available, get_per_qs_bounds will return [size_t(-1),size_t(-1)], which align is supposed to interpret as 'full range'
                pav.align( qs_.cseq_at(i).begin(), qs_.cseq_at(i).end(), sp_.match, sp_.match_cgap, sp_.gap_open, sp_.gap_extend, out_scores.begin(), bounds.first, bounds.second );

                //aligner.align(qs_.pvec_at(i).begin(), qs_.pvec_at(i).end());
                //const vu_scalar_t *score_vec = aligner.get_scores();



                //ncup += block.num_valid * block.ref_len * qs_.pvec_at(i).size();
#if 0 // test against old version
                align_pvec_score_vec<vu_scalar_t, VW>( pvec_prof.begin(), pvec_prof.end(), aux_prof.begin(), qs_.pvec_at(i).begin(), qs_.pvec_at(i).end(), score_match, score_match_cgap, score_gap_open, score_gap_extend, out_scores2.begin(), arrays );
                bool eq = std::equal( out_scores.begin(), out_scores.end(), out_scores2.begin() );


                if( !eq ) {
                    std::cout << "meeeeeeep!\n";
                }
                //std::cout << "eq: " << eq << "\n";
#endif

//                 std::cout << "scores: ";
//                 std::copy( out_scores.begin(), out_scores.end(), std::ostream_iterator<int>(std::cout, "\n" ) );
//                 std::cout << "\n";
                results_.offer( i, block.edges, block.edges + block.num_valid, out_scores.begin() );

            }

            ncup += block.num_valid * cups_per_ref;
            ncup_short += block.num_valid * cups_per_ref;

            ticks_all += pav.ticks_all();
            ticks_all_short += pav.ticks_all();

            inner_iters += pav.inner_iters_all();
            inner_iters_short += pav.inner_iters_all();

            if( rank_ == 0 &&  tprint.elapsed() > 10 ) {

                //std::cout << "thread " << rank_ << " " << ncup << " in " << tstatus.elapsed() << " : "
                
                float fdone = (init_queue_size - queue_size) / float(init_queue_size);
                
                lout << fdone * 100 << "% done. ";
                lout << ncup / (tstatus.elapsed() * 1e9) << " gncup/s, " << ticks_all / double(inner_iters) << " tpili (short: " << ncup_short / (tprint.elapsed() * 1e9) << ", " << ticks_all_short / double(inner_iters_short) << ")" << std::endl;

                ncup_short = 0;
                ticks_all_short = 0;
                inner_iters_short = 0;

                tprint = ivy_mike::timer();
            }

        }
#else

//        assert( block.gapp_ptrs[0] != 0 );
//        assert( VW == 4 );
//        const align_pvec_gapp_score<4> aligner( block.seqptrs, block.gapp_ptrs, block.ref_len, score_mismatch, score_match_cgap, score_gap_open, score_gap_extend );
//        for( unsigned int i = 0; i < m_pnt.m_qs_names.size(); i++ ) {
//
//            size_t stride = 1;
//            size_t aux_stride = 1;
//
//            aligner.align(m_pnt.m_qs_pvecs[i]);
//            const float *score_vec = aligner.get_scores();
//
//            ncup += block.num_valid * block.ref_len * m_pnt.m_qs_pvecs[i].size();
//            {
//                ivy_mike::lock_guard<ivy_mike::mutex> lock( m_pnt.m_qmtx );
//
//                for( int k = 0; k < block.num_valid; k++ ) {
//
//
//
//                    if( score_vec[k] < m_pnt.m_qs_bestscore[i] || (score_vec[k] == m_pnt.m_qs_bestscore[i] && block.edges[k] < m_pnt.m_qs_bestedge[i] )) {
//                        const bool validate = false;
//                        if( validate ) {
//                            const int *seqptr = block.seqptrs[k];
//                            const double *gapp_ptr = block.gapp_ptrs[k];
//
////                                std::vector<double> gapp_tmp(gapp_ptr, gapp_ptr + block.ref_len);
//
//
//                            pars_align_gapp_seq pas( seqptr, m_pnt.m_qs_pvecs[i].data(), block.ref_len, m_pnt.m_qs_pvecs[i].size(), stride, gapp_ptr, aux_stride, seq_arrays_gapp, 0, score_gap_open, score_gap_extend, score_mismatch, score_match_cgap );
//                            int res = pas.alignFreeshift(INT_MAX);
//
//                            if( res != score_vec[k] ) {
//
//
//                                std::cout << "meeeeeeep! score: " << score_vec[k] << " " << res << "\n";
//                            }
//                        }
//
//                        m_pnt.m_qs_bestscore[i] = score_vec[k];
//                        m_pnt.m_qs_bestedge[i] = block.edges[k];
//                    }
//                }
//            }
//        }
//    }

#endif
        {
            ivy_mike::lock_guard<ivy_mike::mutex> lock( *block_queue_.hack_mutex() );
            lout << "thread " << rank_ << ": " << ncup / (tstatus.elapsed() * 1e9) << " gncup/s" << std::endl;
        }
    }
};


template <typename pvec_t,typename seq_tag>
void driver<pvec_t,seq_tag>::calc_scores(size_t n_threads, const my_references& refs, const my_queries& qs, scoring_results* res, const papara_score_parameters& sp) {

    //
    // build the alignment blocks
    //


    block_queue<seq_tag> bq;
    build_block_queue(refs, &bq);

    //
    // work
    //
    ivy_mike::timer t1;
    ivy_mike::thread_group tg;
	lout << "papara_core version " << papara::get_version_string() << std::endl;
    lout << "start scoring, using " << n_threads <<  " threads" << std::endl;

    typedef worker<seq_tag> worker_t;

    for( size_t i = 1; i < n_threads; ++i ) {
        tg.create_thread(worker_t(&bq, res, qs, i, sp));
    }

    worker_t w0(&bq, res, qs, 0, sp );

    w0();

    tg.join_all();

    lout << "scoring finished: " << t1.elapsed() << std::endl;

}

template <typename pvec_t,typename seq_tag>
void driver<pvec_t,seq_tag>::do_newview(pvec_t& root_pvec, lnode* n1, lnode* n2, bool incremental) {
    typedef my_adata_gen<pvec_t, seq_tag > my_adata;

    std::deque<rooted_bifurcation<im_tree_parser::lnode> > trav_order;

    //std::cout << "traversal for branch: " << *(n1->m_data) << " " << *(n2->m_data) << "\n";

    rooted_traversal_order( n1, n2, trav_order, incremental );
    //     std::cout << "traversal: " << trav_order.size() << "\n";

    for( std::deque< rooted_bifurcation< ivy_mike::tree_parser_ms::lnode > >::iterator it = trav_order.begin(); it != trav_order.end(); ++it ) {
        //         std::cout << *it << "\n";

        my_adata *p = it->parent->m_data->get_as<my_adata>();
        my_adata *c1 = it->child1->m_data->get_as<my_adata>();
        my_adata *c2 = it->child2->m_data->get_as<my_adata>();
        //         rooted_bifurcation<ivy_mike::tree_parser_ms::lnode>::tip_case tc = it->tc;



        //         std::cout << "tip case: " << (*it) << "\n";


        pvec_t::newview(p->get_pvec(), c1->get_pvec(), c2->get_pvec(), it->child1->backLen, it->child2->backLen, it->tc);

    }





    {
        my_adata *c1 = n1->m_data->get_as<my_adata>();
        my_adata *c2 = n2->m_data->get_as<my_adata>();

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


template <typename pvec_t,typename seq_tag>
void driver<pvec_t,seq_tag>::build_block_queue(const my_references& refs, my_block_queue* bq) {
    // creates the list of ref-block to be consumed by the worker threads.  A ref-block onsists of N ancestral state sequences, where N='width of the vector unit'.
    // The vectorized alignment implementation will align a QS against a whole ref-block at a time, rather than a single ancestral state sequence as in the
    // sequencial algorithm.

    const static size_t VW = vu_config<seq_tag>::width;

    typedef typename block_queue<seq_tag>::block_t block_t;


    size_t n_groups = (refs.num_pvecs() / VW);
    if( (refs.num_pvecs() % VW) != 0 ) {
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
            if( edge < refs.num_pvecs()) {
                block.edges[i] = edge;
                block.num_valid++;

                block.seqptrs[i] = refs.pvec_at(edge).data();
                block.auxptrs[i] = refs.aux_at(edge).data();

                //                if( !m_ref_gapp[edge].empty() ) {
                //                    block.gapp_ptrs[i] = m_ref_gapp[edge].data();
                //                } else {
                block.gapp_ptrs[i] = 0;
                //                }

                block.ref_len = refs.pvec_size();
                //                     do_newview( root_pvec, m_ec.m_edges[edge].first, m_ec.m_edges[edge].second, true );
                //                     root_pvec.to_int_vec(seqlist[i]);
                //                     root_pvec.to_aux_vec(auxlist[i]);
                //
                //                     seqptrs[i] = seqlist[i].data();
                //                     auxptrs[i] = auxlist[i].data();

                num_valid++;
            } else {
                if( i < 1 ) {
                    std::cout << "edge: " << edge << " " << refs.num_pvecs() << std::endl;

                    throw std::runtime_error( "bad integer mathematics" );
                }
                block.edges[i] = block.edges[i-1];

                block.seqptrs[i] = block.seqptrs[i-1];
                block.auxptrs[i] = block.auxptrs[i-1];
                block.gapp_ptrs[i] = block.gapp_ptrs[i-1];
            }

        }
        bq->push_back(block);
    }
}

template <typename pvec_t,typename seq_tag>
void driver<pvec_t,seq_tag>::seq_to_position_map(const std::vector< uint8_t >& seq, std::vector< int >& map) {
    typedef model<seq_tag> seq_model;

    for( size_t i = 0; i < seq.size(); ++i ) {
        if( seq_model::pstate_is_single(seq_model::s2p(seq[i]))) {
            map.push_back(int(i));
        }
    }
}

template <typename pvec_t,typename seq_tag>
void driver<pvec_t,seq_tag>::print_best_scores(std::ostream& os, const my_queries& qs, const scoring_results& res) {
    boost::io::ios_all_saver ioss(os);
    os << std::setfill ('0');
    for( unsigned int i = 0; i < qs.size(); i++ ) {
        os << qs.name_at(i) << " "  << std::setw (4) << res.bestedge_at(i) << " " << std::setw(5) << res.bestscore_at(i) << "\n";
    }
}

template <typename pvec_t,typename seq_tag>
std::vector< std::vector< uint8_t > > driver<pvec_t,seq_tag>::generate_traces(std::ostream& os_quality, std::ostream& os_cands, const my_queries& qs, const my_references& refs, const scoring_results& res, const papara_score_parameters& sp) {


    typedef typename queries<seq_tag>::pars_state_t pars_state_t;
    typedef model<seq_tag> seq_model;


    lout << "generating best scoring alignments\n";
    ivy_mike::timer t1;



    std::vector<pars_state_t> out_qs_ps;
    align_arrays_traceback<int> arrays;

    std::vector<std::vector<uint8_t> > qs_traces( qs.size() );

    std::vector<uint8_t> cand_trace;

    std::deque<size_t> bounded_bad_scores;
    
    for( size_t i = 0; i < qs.size(); i++ ) {
        size_t best_edge = res.bestedge_at(i);

        assert( size_t(best_edge) < refs.num_pvecs() );

        int score = -1;



        const std::vector<pars_state_t> &qp = qs.pvec_at(i);

        score = align_freeshift_pvec<int>(
                    refs.pvec_at(best_edge).begin(), refs.pvec_at(best_edge).end(),
                    refs.aux_at(best_edge).begin(),
                    qp.begin(), qp.end(),
                    sp.match, sp.match_cgap, sp.gap_open, sp.gap_extend, qs_traces.at(i), arrays
                );


//         std::cout << "scores: " << score << " " << res.bestscore_at(i) << "\n";
        
        std::pair<size_t,size_t> bounds = qs.get_per_qs_bounds( i );
        
        
        if( bounds.first == size_t(-1) ) {
            if( score != res.bestscore_at(i) ) {
                std::cout << "meeeeeeep! score: " << res.bestscore_at(i) << " " << score << "\n";
                throw std::runtime_error( "alignment scores differ between the vectorized and sequential alignment kernels.");
            }
        } else {
            if( score != res.bestscore_at(i) ) {
                bounded_bad_scores.push_back(i);
            }
        }





        if( os_cands.good() ) {
            const scoring_results::candidates &cands = res.candidates_at(i);

            std::vector<std::vector<uint8_t> > unique_traces;

            for( size_t j = 0; j < cands.size(); ++j ) {
                const scoring_results::candidate &cand = cands[j];

                cand_trace.clear();

                score = align_freeshift_pvec<int>(
                            refs.pvec_at(cand.ref()).begin(), refs.pvec_at(cand.ref()).end(),
                            refs.aux_at(cand.ref()).begin(),
                            qp.begin(), qp.end(),
                            sp.match, sp.match_cgap, sp.gap_open, sp.gap_extend, cand_trace, arrays
                        );
                out_qs_ps.clear();

                std::vector<std::vector<uint8_t> >::iterator it;// = unique_traces.end(); //

                if( !true ) {
                    it = std::lower_bound(unique_traces.begin(), unique_traces.end(), cand_trace );
                } else {
                    it = unique_traces.end();
                }

                if( it == unique_traces.end() || *it != cand_trace ) {

                    gapstream_to_alignment_no_ref_gaps(cand_trace, qp, &out_qs_ps, seq_model::gap_pstate() );

                    unique_traces.insert(it, cand_trace);
                    os_cands << i << " " << cand.ref() << " " << cand.score() << "\t";

                    //os << std::setw(pad) << std::left << qs.name_at(i);
                    std::transform( out_qs_ps.begin(), out_qs_ps.end(), std::ostream_iterator<char>(os_cands), seq_model::p2s );
                    os_cands << std::endl;
                }

            }
        }

    }
    
    if( !bounded_bad_scores.empty() ) {
        std::cout << "There were internal problems handling per-gene QS. This is most likely due to overhangs into another partition. The overhangs will be chopped off, but the alignment may be wrong.\n";
    
        std::cout << "QS names";
            
            if( bounded_bad_scores.size() > 20 ) {
                std::cout << " (showing only first 20 of " << bounded_bad_scores.size() << " QS names):\n";
            } else {
                std::cout << " :\n";
            }
            
            size_t m = std::min( bounded_bad_scores.size(), size_t(20) );
                        
            for( size_t i = 0; i < m; ++i ) {
                std::cout << qs.name_at( bounded_bad_scores[i] ) << "\n";
            }
        
    }


    return qs_traces;
}

template <typename pvec_t,typename seq_tag>
void driver<pvec_t,seq_tag>::align_best_scores2(std::ostream& os, std::ostream& os_quality, std::ostream& os_cands, const my_queries& qs, const my_references& refs, const scoring_results& res, size_t pad, const bool ref_gaps, const papara_score_parameters& sp) {

    typedef typename queries<seq_tag>::pars_state_t pars_state_t;
    typedef model<seq_tag> seq_model;


    std::vector<std::vector<uint8_t> > qs_traces = generate_traces(os_quality, os_cands, qs, refs, res, sp );
    std::vector<pars_state_t> out_qs_ps;
    for( size_t i = 0; i < qs.size(); ++i ) {
        const std::vector<pars_state_t> &qp = qs.pvec_at(i);

        out_qs_ps.clear();


        gapstream_to_alignment_no_ref_gaps(qs_traces.at(i), qp, &out_qs_ps, seq_model::gap_pstate() );
        //boost::dynamic_bitset<> bs( out_qs_ps.size() );
        std::vector<bool> bs( out_qs_ps.size() );
        std::transform( out_qs_ps.begin(), out_qs_ps.end(), bs.begin(), seq_model::pstate_is_gap );
    }


}

template <typename pvec_t,typename seq_tag>
void driver<pvec_t,seq_tag>::align_best_scores(std::ostream& os, std::ostream& os_quality, std::ostream& os_cands, const my_queries& qs, const my_references& refs, const scoring_results& res, size_t pad, const bool ref_gaps, const papara_score_parameters& sp) {
    // create the actual alignments for the best scoring insertion position (=do the traceback)

    typedef typename queries<seq_tag>::pars_state_t pars_state_t;
    typedef model<seq_tag> seq_model;


    //     lout << "generating best scoring alignments\n";
    //     ivy_mike::timer t1;



    double mean_quality = 0.0;
    double n_quality = 0.0;
    //


    std::vector<pars_state_t> out_qs_ps;

    // create the best alignment traces per qs


    std::vector<std::vector<uint8_t> > qs_traces = generate_traces(os_quality, os_cands, qs, refs, res, sp );


    // collect ref gaps introduiced by qs
    ref_gap_collector rgc( refs.pvec_size() );
    for( std::vector<std::vector<uint8_t> >::iterator it = qs_traces.begin(); it != qs_traces.end(); ++it ) {
        rgc.add_trace(*it);
    }




    if( ref_gaps ) {
        os << refs.num_seqs() + qs.size() << " " << rgc.transformed_ref_len() << "\n";
    } else {
        os << refs.num_seqs() + qs.size() << " " << refs.pvec_size() << "\n";
    }
    // write refs (and apply the ref gaps)

    for( size_t i = 0; i < refs.num_seqs(); i++ ) {
        os << std::setw(pad) << std::left << refs.name_at(i);

        if( ref_gaps ) {
            rgc.transform( refs.seq_at(i).begin(), refs.seq_at(i).end(), std::ostream_iterator<char>(os), '-' );
        } else {
            std::transform( refs.seq_at(i).begin(), refs.seq_at(i).end(), std::ostream_iterator<char>(os), seq_model::normalize);
        }

        //std::transform( m_ref_seqs[i].begin(), m_ref_seqs[i].end(), std::ostream_iterator<char>(os), seq_model::normalize );
        os << "\n";
    }



    for( size_t i = 0; i < qs.size(); i++ ) {

        const std::vector<pars_state_t> &qp = qs.pvec_at(i);

        out_qs_ps.clear();

        if( ref_gaps ) {
            gapstream_to_alignment(qs_traces.at(i), qp, &out_qs_ps, seq_model::gap_pstate(), rgc);
        } else {
            gapstream_to_alignment_no_ref_gaps(qs_traces.at(i), qp, &out_qs_ps, seq_model::gap_pstate() );
        }

        os << std::setw(pad) << std::left << qs.name_at(i);
        std::transform( out_qs_ps.begin(), out_qs_ps.end(), std::ostream_iterator<char>(os), seq_model::p2s );
        os << std::endl;





        if( os_quality.good() && qs.seq_at(i).size() == refs.pvec_size()) {


            std::vector<int> map_ref;
            std::vector<int> map_aligned;
            seq_to_position_map( qs.seq_at(i), map_ref );
            align_utils::trace_to_position_map( qs_traces[i], &map_aligned );


            if( map_ref.size() != map_aligned.size() ) {
                throw std::runtime_error( "alignment quirk: map_ref.size() != map_aligned.size()" );
            }

            size_t num_equal = ivy_mike::count_equal( map_ref.begin(), map_ref.end(), map_aligned.begin() );

            //std::cout << "size: " << map_ref.size() << " " << map_aligned.size() << " " << m_qs_seqs[i].size() << "\n";
            //std::cout << num_equal << " equal of " << map_ref.size() << "\n";

            double score = num_equal / double(map_ref.size());
            //double score = alignment_quality( out_qs, m_qs_seqs[i], debug );

            os_quality << qs.name_at(i) << " " << score << "\n";

            mean_quality += score;
            n_quality += 1;
        }


    }

    lout << "mean quality: " << mean_quality / n_quality << "\n";

}

template <typename pvec_t,typename seq_tag>
void driver<pvec_t,seq_tag>::align_best_scores_oa( output_alignment *oa, const my_queries &qs, const my_references &refs, const scoring_results &res, size_t pad, const bool ref_gaps, const papara_score_parameters &sp ) {
    typedef typename queries<seq_tag>::pars_state_t pars_state_t;
    typedef model<seq_tag> seq_model;


    //     lout << "generating best scoring alignments\n";
    //     ivy_mike::timer t1;



    double mean_quality = 0.0;
    double n_quality = 0.0;
    //


    std::vector<pars_state_t> out_qs_ps;

    

    // supply non-open ofstreams to keep it quiet. This is actually quite a bad interface...
    std::ofstream os_quality;
    std::ofstream os_cands;
    
    
    // create the best alignment traces per qs
    std::vector<std::vector<uint8_t> > qs_traces = generate_traces(os_quality, os_cands, qs, refs, res, sp );


    // collect ref gaps introduiced by qs
    ref_gap_collector rgc( refs.pvec_size() );
    for( std::vector<std::vector<uint8_t> >::iterator it = qs_traces.begin(); it != qs_traces.end(); ++it ) {
        rgc.add_trace(*it);
    }




    if( ref_gaps ) {
        oa->set_size(refs.num_seqs() + qs.size(), rgc.transformed_ref_len());
    } else {
        oa->set_size(refs.num_seqs() + qs.size(), refs.pvec_size());
    }
    oa->set_max_name_length( pad );
    
    // write refs (and apply the ref gaps)

    std::vector<char> tmp;
    for( size_t i = 0; i < refs.num_seqs(); i++ ) {
        tmp.clear();
        
        

        if( ref_gaps ) {
            rgc.transform( refs.seq_at(i).begin(), refs.seq_at(i).end(), std::back_inserter(tmp), '-' );
        } else {
            std::transform( refs.seq_at(i).begin(), refs.seq_at(i).end(), std::back_inserter(tmp), seq_model::normalize);
        }

        oa->push_back( refs.name_at(i), tmp, output_alignment::type_ref );
        //std::transform( m_ref_seqs[i].begin(), m_ref_seqs[i].end(), std::ostream_iterator<char>(os), seq_model::normalize );
        
    }

    std::deque<size_t> overhang_qs;

    for( size_t i = 0; i < qs.size(); i++ ) {
        tmp.clear();
        const std::vector<pars_state_t> &qp = qs.pvec_at(i);

        out_qs_ps.clear();

        if( ref_gaps ) {
            gapstream_to_alignment(qs_traces.at(i), qp, &out_qs_ps, seq_model::gap_pstate(), rgc);
        } else {
            gapstream_to_alignment_no_ref_gaps(qs_traces.at(i), qp, &out_qs_ps, seq_model::gap_pstate() );
            
            
            
            // chop off QS parts that hang over into other partition
            std::pair<size_t,size_t> bounds = qs.get_per_qs_bounds(i);
            
            if( bounds.first != size_t(-1) ) {
                assert( bounds.second != size_t(-1) );
                
                assert( bounds.first < out_qs_ps.size() );
                assert( bounds.second <= out_qs_ps.size() );
                assert( bounds.first < bounds.second );
                bool pre_overhang = false;
                bool post_overhang = false;
                
                const typename seq_model::pars_state_t gap_state = seq_model::gap_pstate();
                for( size_t j = 0; j < bounds.first; ++j ) {
                    if( out_qs_ps[j] != gap_state ) {
                        out_qs_ps[j] = gap_state;
                        pre_overhang = true;
                    }
                }
                
                for( size_t j = bounds.second + 1; j < out_qs_ps.size(); ++j ) {
                    if( out_qs_ps[j] != gap_state ) {
                        out_qs_ps[j] = gap_state;
                        post_overhang = true;
                    }
                }
                
                if( pre_overhang || post_overhang ) {
                    overhang_qs.push_back(i);
                }
                
            }
            
        }

        if( !overhang_qs.empty() ) {
            std::cout << "WARNING: per-gene alignment, with overhangs into other partitons. chopped off.\nQS names";
            
            if( overhang_qs.size() > 20 ) {
                std::cout << " (showing only first 20 of " << overhang_qs.size() << " QS names):\n";
            } else {
                std::cout << " :\n";
            }
            
            size_t m = std::min( overhang_qs.size(), size_t(20) );
            for( size_t i = 0; i < m; ++i ) {
                std::cout << qs.name_at( overhang_qs[i] ) << "\n";
            }
        }
        
        //os << std::setw(pad) << std::left << qs.name_at(i);
        std::transform( out_qs_ps.begin(), out_qs_ps.end(), std::back_inserter(tmp), seq_model::p2s );
        
        oa->push_back( qs.name_at(i), tmp, output_alignment::type_qs );





        if( os_quality.good() && qs.seq_at(i).size() == refs.pvec_size()) {


            std::vector<int> map_ref;
            std::vector<int> map_aligned;
            seq_to_position_map( qs.seq_at(i), map_ref );
            align_utils::trace_to_position_map( qs_traces[i], &map_aligned );


            if( map_ref.size() != map_aligned.size() ) {
                throw std::runtime_error( "alignment quirk: map_ref.size() != map_aligned.size()" );
            }

            size_t num_equal = ivy_mike::count_equal( map_ref.begin(), map_ref.end(), map_aligned.begin() );

            //std::cout << "size: " << map_ref.size() << " " << map_aligned.size() << " " << m_qs_seqs[i].size() << "\n";
            //std::cout << num_equal << " equal of " << map_ref.size() << "\n";

            double score = num_equal / double(map_ref.size());
            //double score = alignment_quality( out_qs, m_qs_seqs[i], debug );

            os_quality << qs.name_at(i) << " " << score << "\n";

            mean_quality += score;
            n_quality += 1;
        }


    }
}
void output_alignment_phylip::write_seq_phylip(const std::string& name, const out_seq& seq) {
    size_t pad = std::max( max_name_len_, name.size() + 1 );

    os_ << std::setw(pad) << std::left << name;

    std::copy( seq.begin(), seq.end(), std::ostream_iterator<char>(os_) );
    os_ << "\n";
}
void output_alignment_phylip::push_back(const std::string& name, const out_seq& seq, output_alignment::seq_type t) {
    if( !header_flushed_ ) {
        os_ << num_rows_ << " " << num_cols_ << "\n";
        header_flushed_ = true;
    }
    
    write_seq_phylip( name, seq );
}

void output_alignment_fasta::push_back(const std::string& name, const out_seq& seq, output_alignment::seq_type t) {
    os_ << ">" << name << "\n";
    
    std::copy( seq.begin(), seq.end(), std::ostream_iterator<char>(os_) );
    os_ << "\n";
}
output_alignment::~output_alignment() {}

double alignment_quality_very_strict(const std::vector<uint8_t> &s1, const std::vector<uint8_t> &s2, bool debug) {
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

double alignment_quality(const std::vector<uint8_t> &s1, const std::vector<uint8_t> &s2, bool debug) {
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

std::string filename(const std::string &run_name, const char *type) {
    std::stringstream ss;

    ss << "papara_" << type << "." << run_name;

    return ss.str();
}

bool file_exists(const char *filename)
{
    std::ifstream is(filename);
    return is.good();
}

}

namespace papara
{
    // explicit template instantiations of the queries/references classes for the different supported data types and gap models,
    // i.e., controlled combinatorial detonation happens here...
    template class queries<tag_dna>;
    template class queries<tag_aa>;

    template class references<pvec_pgap,tag_dna>;
    template class references<pvec_cgap,tag_dna>;

    template class references<pvec_cgap,tag_aa>;
    template class references<pvec_pgap,tag_aa>;


    template class driver<pvec_pgap,tag_dna>;
    template class driver<pvec_cgap,tag_dna>;

    template class driver<pvec_cgap,tag_aa>;
    template class driver<pvec_pgap,tag_aa>;
}
