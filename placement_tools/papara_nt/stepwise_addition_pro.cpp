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

// #include "pch.h"

// someone down the line tries to sabotage us by including a windows header, 
// which will clutter the whole namespace with macros. I don't know who it is, 
// but disable_shit.h takes care of restricting the damage...
#include "ivymike/disable_shit.h"
#include <cctype>

#include <algorithm>
#include <functional>
#include <vector>
#include <iostream>
#include <fstream>
#include <iterator>
#include <deque>
#include <thread>
#include <iomanip>
#define BOOST_UBLAS_NDEBUG


#include <boost/bind.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "pairwise_seq_distance.h"
#include "stepwise_align.h"
#include "raxml_interface.h"
#include "sequence_model.h"
#include "pvec.h"
#include "ivymike/time.h"
#include "ivymike/getopt.h"
#include "ivymike/smart_ptr.h"
#include "ivymike/tree_parser.h"
#include "ivymike/tdmatrix.h"
#include "ivymike/algorithm.h"
#include "ivymike/tree_traversal_utils.h"
#include "ivymike/tree_split_utils.h"
#include "ivymike/flat_map.h"
#include "ivymike/fasta.h"
#include "ivymike/thread.h"

namespace tree_parser = ivy_mike::tree_parser_ms;

using ivy_mike::tree_parser_ms::ln_pool;
using ivy_mike::tree_parser_ms::lnode;

using sequence_model::tag_aa;
using sequence_model::tag_dna;

using ivy_mike::apply_lnode;
using ivy_mike::apply_lnode;
using ivy_mike::back_insert_ifer;
using ivy_mike::iterate_lnode;
using ivy_mike::rooted_bifurcation;
using ivy_mike::back_insert_ifer;
using ivy_mike::iterate_lnode;
using ivy_mike::rooted_bifurcation;
using ivy_mike::scoring_matrix;

namespace ublas = boost::numeric::ublas;


typedef std::vector<unsigned char> sequence;

namespace {
 
//     double g_delta = log(0.1);
//     double g_epsilon = log(0.5);
}

class log_odds {
public:

    log_odds( double bg_prob ) : bg_prob_(bg_prob) {}

    inline double operator()( double p ) {
        return std::max( -100.0, log( p / bg_prob_ ));
    }

private:
    const double bg_prob_;
};


static std::vector<uint8_t> gapstream_to_alignment( const std::vector<uint8_t> &gaps, const std::vector<uint8_t> &raw, uint8_t gap_char, bool upper ) {

    std::vector<uint8_t> out;
    
    std::vector<uint8_t>::const_reverse_iterator rit = raw.rbegin();


    // 'gap indicator': if upper is set, insert gap into reference (=if gaps[i] == 2).
    const uint8_t gap_ind = upper ? 2 : 1;



    for ( std::vector<uint8_t>::const_iterator git = gaps.begin(); git != gaps.end(); ++git ) {



        if ( *git == gap_ind ) {
            out.push_back(gap_char);
        } else {
            
            if( rit != raw.rend() ) {
                out.push_back(*rit);
                ++rit;
            } else {
                out.push_back( 'X' );
            }
        }
    }
    if( rit != raw.rend() ) {

        std::cerr << "too short tb: " << raw.rend() - rit << " upper: " << upper << "\n";
    }
    std::reverse(out.begin(), out.end());
    
    return out;
}

class log_odds_viterbi {
    typedef ublas::matrix<double> dmat;
    typedef std::vector<double> dsvec;


    // lof_t: log-odds-float = float type good enough to hold/calculate log-odds scores.
    // 32bit float should be enough
    typedef float lof_t;
    typedef ublas::matrix<lof_t> lomat;

    typedef std::vector<lof_t> losvec;

public:
    log_odds_viterbi( const dmat &state, const dmat &gap, boost::array<double,4> state_freq )
    : 
      ref_state_prob_(state), ref_gap_prob_(gap), ref_len_(state.size2()),
      state_freq_(state_freq),
      neg_inf_( -std::numeric_limits<lof_t>::infinity() ),
      m_(state.size2() + 1),
      d_(state.size2() + 1),
      i_(state.size2() + 1),
      max_matrix_height_(0),
      delta_(lof_t(log(0.01))),
      epsilon_(lof_t(log(0.1)))
    {
        ref_gap_prob_log_ = ref_gap_prob_;
        {
            auto &d = ref_gap_prob_log_.data();
            //std::transform( d.begin(), d.end(), d.begin(), log );
			std::transform( d.begin(), d.end(), d.begin(), [](double v){ return log(v); } );
            
        }
        
//         for( auto it1 = ref_gap_prob_.begin1(); it1 != ref_gap_prob_.end1(); ++it1 ) {
//             std::copy( it1.begin(), it1.end(), std::ostream_iterator<float>( std::cout, " " ));
//             std::cout << "\n";
//             
//         }
        
        precalc_log_odds();
    }

    void setup( size_t qlen ) {
        assert( ref_gap_prob_.size2() == ref_len_ );



        // init first rows
        std::fill( m_.begin(), m_.end(), 0.0 );
        
        
        //std::fill( d_.begin(), d_.end(), 0.0 );
        {
            int i = 0; // TODO: still wearing my STL hat. Not really sure if std::generate + lambda is better than a for loop...
            std::generate( d_.begin(), d_.end(), [&]() { return delta_ + (i++) * epsilon_; });
        
            m_.assign( d_.begin(), d_.end() );
            
        }
        std::fill( i_.begin(), i_.end(), neg_inf_ );

        std::fill( traceback_.begin1().begin(), traceback_.begin1().end(), tb_d_to_d | tb_m_to_d );
        traceback_(0, 0) = tb_m_to_m;
        
        // init first columns
        m_[0] = 0.0;
        d_[0] = neg_inf_;
        i_[0] = delta_;



    }


    void precalc_log_odds() {
        ref_state_lo_.resize( ref_state_prob_.size1(), ref_state_prob_.size2() );

        for( size_t i = 0; i < 4; ++i ) {
            const ublas::matrix_row<dmat> prow( ref_state_prob_, i );
            ublas::matrix_row<lomat> lorow( ref_state_lo_, i );
            std::transform( prow.begin(), prow.end(), lorow.begin(), log_odds(state_freq_[i]));
        }



//         {
//             odds odds_ngap( 1 - g_gap_freq );
//             odds odds_gap( g_gap_freq );
// 
//             ref_ngap_odds_.resize(ref_gap_prob_.size());
//             ref_gap_odds_.resize(ref_gap_prob_.size());
//             for( size_t i = 0; i < ref_gap_prob_.size(); ++i ) {
//                 ref_ngap_odds_[i] = odds_ngap(1 - ref_gap_prob_[i]);
//                 ref_gap_odds_[i] = odds_gap( ref_gap_prob_[i] );
//             }
//         }
    }

    template<typename T, typename T2>
    static inline std::pair<T,T2> max2( const T &a, const T &b, const T2 &a2, const T2 &b2 ) {
        if( a > b ) {
            return std::make_pair( a, a2 );
        } else {
            return std::make_pair( b, b2 );
        }
    }
    
    
    
    template<typename T, typename T2>
    static inline std::pair<T,T2> max3( const T &a, const T &b, const T &c, const T2 &a2, const T2 &b2, const T2 &c2, bool dump = false ) {
        //return std::max( a, std::max( b, c ));
        if( dump ) {
            std::cout << "(" << a << " " << b << " " << c << ")";
        }
        
        if( a > b ) {
            return max2( a, c, a2, c2 );
        } else {
            return max2( b, c, b2, c2 );
        }
    }
    
            
    
    

    double align( const std::vector<uint8_t> &qs ) {
        const size_t qlen = qs.size();

        traceback_.resize(qlen + 1, m_.size(), false);
        
        setup( qlen );

        //dmat ref_state_trans = trans(ref_state_prob_);


        assert( m_.size() == ref_len_ + 1 );

        best_score_ = neg_inf_;
        best_i_ = -1;
        best_j_ = -1;
        short last_state = 0;

        for( size_t i = 1; i < qlen + 1; ++i ) {
            const int b = qs[i-1];
            //          std::cout << "b: " << b << "\n";

            //const double b_freq = state_freq_.at(b);
            //const ublas::matrix_column<dmat> b_state( ref_state_prob_, b );
            const ublas::matrix_row<lomat> b_state_lo( ref_state_lo_, b );

            //          const ublas::matrix_column<dmat> ngap_prob( ref_gap_prob_, 0 );
            //          const ublas::matrix_column<dmat> gap_prob( ref_gap_prob_, 1 );

            i_[0] = delta_ + (i-1) * epsilon_;
            m_[0] = i_[0];
            
            lof_t diag_m = m_[0];
            lof_t diag_d = d_[0];
            lof_t diag_i = i_[0];

            losvec::iterator m0 = m_.begin() + 1;
            losvec::iterator d0 = d_.begin() + 1;
            losvec::iterator i0 = i_.begin() + 1;

            losvec::iterator m1 = m_.begin();
            losvec::iterator d1 = d_.begin();
            losvec::iterator i1 = i_.begin();
            ublas::matrix_row<lomat>::const_iterator bsl = b_state_lo.begin();
            
            ublas::matrix_row<ublas::matrix<short>> traceback_row( traceback_, i );
            auto traceback_iter = traceback_row.begin();
            *traceback_iter = tb_i_to_i | tb_m_to_i;
            ++traceback_iter;
            
//             losvec::iterator rg = ref_gap_odds_.begin();
//             losvec::iterator rng = ref_ngap_odds_.begin();

            auto gap_prob_log = ref_gap_prob_log_.begin2();
            
            const losvec::iterator m_end = m_.end();

            for( size_t j = 1; m0 != m_end; m1 = m0++, d1 = d0++, i1 = i0++, ++bsl, ++gap_prob_log, ++traceback_iter, ++j ) {
                //ublas::matrix_row<dmat> a_state(ref_state_prob_, j-1 );
                //ublas::matrix_row<dmat> a_gap(ref_gap_prob_, j-1 );

                //double match_log_odds = log( b_state[j-1] / b_freq );
                //lof_t match_log_odds = b_state_lo[j-1];
                const lof_t match_log_odds = *bsl;



                //lof_t gap_log_odds = ref_gap_lo_[j-1];
                //lof_t ngap_log_odds = ref_ngap_lo_[j-1];
//                 const lof_t gap_odds = *rg;
//                 const lof_t ngap_odds = *rng;


                lof_t p_ngap_log = *gap_prob_log;
                lof_t p_gap_log = *(gap_prob_log.begin() + 1);
                
//                 p_gap = std::max( p_gap, 0.01f );
                
//                 lof_t m_log_sum = log(
//                            exp(diag_m) * p_ngap
//                          + exp(diag_d) * p_gap
//                          + exp(diag_i)
//                 );
                
             //   std::cout << "logp: " << log(p_gap) << " ";
//                 std::cout << "\t";
                auto m_log_max = max3<float>(
                           diag_m + p_ngap_log,
                           diag_d,
                           diag_i,
                           tb_m_to_m,
                           tb_m_to_d,
                           tb_m_to_i,
                           !true
                );
                
                diag_m = *m0;
                *m0 = m_log_max.first + match_log_odds;
                *traceback_iter = m_log_max.second; // do not or, because it's the first value tritten to the tb matrix
//                 std::cout << *m0 << " ";
                
#if 0
                std::cout << i << " " << j << " " << m_(i,j) << " : " << m_(i-1, j-1) + ngap_log_odds
                        << " " << d_(i-1, j-1) + gap_log_odds << " " << i_(i-1, j-1) + gap_log_odds << " " << match_log_odds << " " << gap_log_odds << " " << ngap_log_odds << " max: " << m_max << "\n";
#endif

                diag_i = *i0;

                // the two 'diags' have already been updated, so they're both actually containing the current 'aboves',
                // which is exactly what we need to calculate the new i

//                 lof_t i_log_sum = log(
//                           exp(diag_m) /** delta_*/
//                         + exp(diag_i) /** epsilon_*/
//                 );
                auto i_log_max = max3<float>(
                          diag_m + delta_,
                          diag_i + epsilon_,
                          diag_d + delta_,
                          tb_i_to_m,
                          tb_i_to_i,
                          tb_i_to_d
                );
                
                
                *i0 = i_log_max.first;
                *traceback_iter |= i_log_max.second;
                
#if 1
//                 lof_t d_log_sum = log(
//                           exp(*m1) /** delta_*/
//                         + exp(*d1) /** epsilon_*/
//                 );
                auto d_log_max = max3<float>(
                          *m1 + delta_,
                          *d1 + epsilon_,
                          *i1 + delta_,
                          tb_d_to_m,
                          tb_d_to_d,
                          tb_d_to_i
                );
#else
                lof_t d_log_sum = *m1 + math_approx::log(
                        delta_
                      + math_approx::exp(*d1 - *m1) * epsilon_
                );
#endif
                diag_d = *d0;
                *d0 = d_log_max.first + p_gap_log;
                *traceback_iter |= d_log_max.second;
               
                last_state = m_log_max.second;
                
//                 if( m0 == m_end - 1 || i == qlen ) {
//                     if( *m0 > best_score_ ) {
//                         best_score_ = *m0;
//                         best_i_ = i;
//                         best_j_ = j;
//                         best_state_ = m_log_max.second;
//                     }
//                     
//                 }
                //end_state_ = max3( *m0, *i0, *d0, tb_m_to_m, tb_m_to_i, tb_m_to_d );
                //std::cout << end_state_.first << " ";
//                 std::cout << match_log_odds << " ";
                //lof_t old_m = m_[j];

            }
//             std::cout << "\n";
          
        }

        //return m_.back();
        
        
        best_score_ = m_.back();
        best_i_ = qlen;
        best_j_ = m_.size() - 1;
        best_state_ = last_state;
        
        return best_score_;
    }
    
    
    std::vector<uint8_t> traceback() {
        std::cout << "end state: " << std::hex << int(end_state_.second) << std::dec << "\n";
        
        int ia = traceback_.size1() - 1;
        int ib = traceback_.size2() - 1;
        
        int max_a = best_i_;
        int max_b = best_j_;
        
        assert( ia == max_a || ib == max_b );
        
        
//         for( auto it1 = traceback_.begin1(); it1 != traceback_.end1(); ++it1 ) {
//             std::cout << std::hex;
//             std::copy( it1.begin(), it1.end(), std::ostream_iterator<int>( std::cout, " " ));
//             std::cout << std::dec << "\n";
//             
//         }
        
        
        bool in_l = false;
        bool in_u = false;
     
        std::vector<uint8_t> tb_out;
        tb_out.reserve( ia + ib );
        
        while( ia > max_a ) {
            tb_out.push_back(2);
            --ia;
        }
        
        while( ib > max_b ) {
            tb_out.push_back(1);
            --ib;
        }
        
        {
            char tb_end = best_state_;
        
            in_u = (tb_end & tb_m_to_i) != 0;
            in_l = (tb_end & tb_m_to_d) != 0;    
            
            assert( !in_u || !in_l );
        }
        
        while( ia > 0 && ib > 0 ) {
            auto tb = traceback_( ia, ib );
//             std::cout << std::hex;
//             std::cout << "tb: " << int(tb) << "\n";
//             std::cout << std::dec;
            if( !in_l && !in_u ) {
                in_l = (tb & tb_m_to_d) != 0;
                in_u = (tb & tb_m_to_i) != 0;
                
//                 if( !in_l && !in_u ) {
//                     
//                     
//                 }
                tb_out.push_back(0);
                --ia;
                --ib;
            } else if( in_u ) {
                tb_out.push_back(2);
                --ia;
                
                in_u = (tb & tb_i_to_i) != 0;
                in_l = (tb & tb_i_to_d) != 0;
            } else if( in_l ) {
                tb_out.push_back(1);
                --ib;
                
                in_l = (tb & tb_d_to_d) != 0;
                in_u = (tb & tb_d_to_i) != 0;
            }
            
            
        }
        
        while( ia > 0 ) {
            tb_out.push_back(2);
            --ia;
        }
        
        while( ib > 0 ) {
            tb_out.push_back(1);
            --ib;
        }
        
        
        //return std::vector<char>( tb_out.rbegin(), tb_out.rend() );
        
        return tb_out;
    }
private:
    ublas::matrix<short> traceback_;
    static const short tb_i_to_i;// = 0x1;
    static const short tb_i_to_m;// = 0x2;
    static const short tb_d_to_d;// = 0x4;
    static const short tb_d_to_m;// = 0x8;
    static const short tb_m_to_m;// = 0x10;
    static const short tb_m_to_i;// = 0x20;
    static const short tb_m_to_d;// = 0x40;
    static const short tb_i_to_d;// = 0x80;
    static const short tb_d_to_i;// = 0x100;
    
    dmat ref_state_prob_;
    dmat ref_gap_prob_;
    dmat ref_gap_prob_log_;
    
    
    lomat ref_state_lo_;
//     losvec ref_gap_odds_;
//     losvec ref_ngap_odds_;

    const size_t ref_len_;
    const boost::array<double,4> state_freq_;

    const float neg_inf_;

    losvec m_;
    losvec d_;
    losvec i_;

    std::pair<lof_t, char> end_state_;
    
    size_t max_matrix_height_;

    const lof_t delta_;// = log(0.1);
    const lof_t epsilon_;// = log(0.5);

//     size_t max_col_;
//     size_t max_row_;
    double best_score_;
    size_t best_i_;
    size_t best_j_;
    short best_state_;



};

const short log_odds_viterbi::tb_i_to_i = 0x1;
const short log_odds_viterbi::tb_i_to_m = 0x2;
const short log_odds_viterbi::tb_d_to_d = 0x4;
const short log_odds_viterbi::tb_d_to_m = 0x8;
const short log_odds_viterbi::tb_m_to_m = 0x10;
const short log_odds_viterbi::tb_m_to_i = 0x20;
const short log_odds_viterbi::tb_m_to_d = 0x40;
const short log_odds_viterbi::tb_i_to_d = 0x80;
const short log_odds_viterbi::tb_d_to_i = 0x100;

namespace {

template<class pvec_t,typename seq_tag>
class my_adata_gen : public ivy_mike::tree_parser_ms::adata {

    
public:
//     int m_ct;
    my_adata_gen() {

//         std::cout << "my_adata\n";

    }

    virtual ~my_adata_gen() {

//         std::cout << "~my_adata\n";

    }

    
    
    void init_sequence( const sequence &seq ) {
        assert( seq_.empty() );
        
        seq_.assign( seq.begin(), seq.end() );
    }
    
    void update_pvec() {
        pvec_.init2( seq_, sequence_model::model<seq_tag>() );
    }
    
    pvec_t &pvec() {
        return pvec_;
    }
    
    const ublas::matrix<float> &calculate_anc_gap_probs() {
        
        const auto &pgap = pvec_.get_pgap();
        
        anc_probs_.resize(pgap.size1(), pgap.size2());
        
        
        auto oit = anc_probs_.begin2();
        for( auto iit = pgap.begin2(); iit != pgap.end2(); ++iit, ++oit ) {
            double v1 = *iit * (1-pvec_pgap::pgap_model->gap_freq());
            double v2 = *(iit.begin() + 1) * pvec_pgap::pgap_model->gap_freq();
            
            *oit = v1 / (v1 + v2);
            *(oit.begin() + 1) = v2 / (v1 + v2);
            
        }
        
        return anc_probs_;
        
//         for( it = 
        
    }
    
private:
    
    

    
    pvec_t pvec_;
    ublas::matrix<float> anc_probs_;
    
    sequence seq_;
    

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
   
typedef sequence_model::model<tag_dna> seq_model;
typedef my_adata_gen<pvec_pgap,tag_dna> my_adata;
}

template<typename iiter> 
bool is_sorted( iiter start, iiter end ) {
    if( std::distance(start,end) < 2 ) {
        return true;
    }
    
    while( start != end - 1) {
        if( start > start++ ) {
            return false;
        }
    }
    
    return true;
}

class tree_deconstructor {
public:
    tree_deconstructor( lnode *t, const std::vector<std::string> &addition_order ) {
        remove_order_.assign( addition_order.rbegin(), addition_order.rend() );
        
        
        ivy_mike::flat_map<std::string,size_t> tip_name_map;

        std::vector<lnode*> edges;
        std::vector<boost::dynamic_bitset<> > splits;
        std::vector<lnode *> sorted_tips;
        ivy_mike::get_all_splits_by_node( t,  edges, splits, sorted_tips );
        
        for( size_t i = 0; i < sorted_tips.size(); ++i ) {
            
            
            tip_name_map.put_fast( sorted_tips[i]->m_data->tipName, i );
        }
        tip_name_map.sort();
        
        assert( remove_order_.size() > 2 );
        
        for( size_t i = 0; i < remove_order_.size(); ++i ) {
            auto name = remove_order_.at(i);
            
            const size_t * idx_ptr = tip_name_map.get(name);
            
            assert( idx_ptr != nullptr );
            lnode *tip = sorted_tips.at(*idx_ptr);
            
            
            if( !tip->m_data->isTip ) {
                tip = tip->back;
            }
            
            assert( tip->m_data->isTip );
            
            lnode *remove_node = tip->back;
            
            
//             ivy_mike::tree_parser_ms::prune_with_rollback pwr(remove_node);
//             pwr_stack_.push_back( std::move(pwr));
             
            pwr_stack_.emplace_back(remove_node);
            
        }
    }
    ~tree_deconstructor() {
        while( !pwr_stack_.empty() ) {
            pwr_stack_.pop_back();
        }
    }
    
    
    lnode *get_save_node() const {
        // return node from the current tree
        return pwr_stack_.back().get_save_node();
    }
    
    void pop() {
        
        assert( !pwr_stack_.empty() );
        pwr_stack_.pop_back();
    }
    bool empty() const {
        return pwr_stack_.empty();
    }
    
    size_t size() const {
        return pwr_stack_.size();
    }
        
  
  
private:
    
    std::vector<ivy_mike::tree_parser_ms::prune_with_rollback> pwr_stack_;
    std::vector<std::string> remove_order_;
};


class addition_order {
public:



    addition_order( const scoring_matrix &sm, const std::vector<sequence> &mapped_seqs )
     : used_seqs_( mapped_seqs.size() )
    {
        const size_t num_seqs = mapped_seqs.size();

        std::cout << "size: " << num_seqs << "\n";
        pw_dist_.init_size(num_seqs, num_seqs);


        ivy_mike::tdmatrix<int> out_scores( mapped_seqs.size(), mapped_seqs.size() );

        const size_t num_ali_threads = 4;
        pairwise_seq_distance(mapped_seqs, out_scores, sm, -5, -2, num_ali_threads, 64);
        init_pw_dist_from_msa_score_matrix(out_scores);

    }


    size_t find_next_candidate() {

        if( pw_dist_.size() == 0 ) {
            throw std::runtime_error( "find_next_candidate called with empty pw-dist matrix");
        }

        if( dist_acc_.empty() ) {
            //
            // create initial distance accumulator if it does not exist.
            //
            size_t f = used_seqs_.find_first();

            std::vector<float> dist_sum;

            while( f != used_seqs_.npos ) {
                ivy_mike::odmatrix<float> slice = pw_dist_[f];
                if( dist_sum.empty() ) {
                    dist_sum.assign( slice.begin(), slice.end() );
                } else {
                    std::transform( dist_sum.begin(), dist_sum.end(), slice.begin(), dist_sum.begin(), std::plus<float>() );
                }

                f = used_seqs_.find_next(f);
            }
            dist_acc_.swap( dist_sum );
        }

        float min_dist = 1e8;
        size_t min_element = size_t(-1);
        for( size_t i = 0; i < dist_acc_.size(); i++ ) {
            if( !used_seqs_[i] && dist_acc_[i] < min_dist ) {
                min_dist = dist_acc_[i];
                min_element = i;
            }
        }

        assert( min_element != size_t(-1) || used_seqs_.count() == used_seqs_.size() );

        if( min_element != size_t(-1) ) {

            // update accumulator
            assert( min_element != size_t(-1));
            ivy_mike::odmatrix<float> slice = pw_dist_[min_element];

            assert( slice.size() == dist_acc_.size() );

            // element-wise calculate dist_acc_ = dist_acc_ + slice;
            std::transform( dist_acc_.begin(), dist_acc_.end(), slice.begin(), dist_acc_.begin(), std::plus<float>() );
            used_seqs_[min_element] = true;
        }


        return min_element;

    }


    std::pair<size_t,size_t> first_pair() const {
        return first_pair_;
    }

private:
    void init_pw_dist_from_msa_score_matrix( ivy_mike::tdmatrix<int> &out_scores ) {
        size_t li = -1, lj = -1;
        float lowest_dist = 1e8;
        int min = *(std::min_element( out_scores.begin(), out_scores.end() ));
        int max = *(std::max_element( out_scores.begin(), out_scores.end() ));

        for( size_t i = 0; i < out_scores.size(); i++ ) {

            for( size_t j = 0; j < out_scores[i].size(); j++ ) {

                // three modes for normalizing: min, max and mean
                //const float norm = min( ma[i][i], ma[j][j] );
                //             const float norm = max( ma[i][i], ma[j][j] );
                const float norm = (out_scores[i][j] - min) / float(max-min);


                const float dist = 1.0 - norm;
                pw_dist_[i][j] = dist;

                if( i != j && dist < lowest_dist ) {
                    lowest_dist = dist;
                    li = i;
                    lj = j;
                }

            }

        }

        used_seqs_[li] = used_seqs_[lj] = true;

        first_pair_ = std::make_pair( li, lj );

    }
    ivy_mike::tdmatrix<float> pw_dist_;
    std::vector<float> dist_acc_;
    boost::dynamic_bitset<> used_seqs_;

    std::pair<size_t,size_t> first_pair_;
   // scoring_matrix scoring_matrix_;
};

template<typename K, typename V>
class flat_map {
public:
    flat_map() : sorted_(true) {}
    
    void sort() {
        if( !sorted_ ) {
            std::sort( pairs_.begin(), pairs_.end() );
            sorted_ = true;
        }
    }
    
    void put_fast( const K &key, const V &value ) {
        pairs_.emplace_back( key, value );
        
        sorted_ = false;
    }
    
    void put_fast( const K &key, V &&value ) {
        pairs_.emplace_back( key, value );
        
        sorted_ = false;
    }
    template<typename... Args>
    void emplace_fast( const K &key, Args&&... args) {
        pairs_.emplace_back( key, args... );
        
        sorted_ = false;
    
    }
    
    
    void put( const K &key, const V &value ) {
        if( !sorted_ ) {
            throw std::runtime_error( "flat_map::put on unsorted map" );
        }
        
        ipair p{key, value};
        auto lb = std::lower_bound( pairs_.begin(), pairs_.end(), p );
        
        pairs_.insert(lb, std::move(p));
    }
    
    const V * get( const K &key ) const {
        if( !sorted_ ) {
            throw std::runtime_error( "flat_map::get on unsorted map" );
        }
        auto lb = std::lower_bound( pairs_.begin(), pairs_.end(), ipair(key, V()) );
        
        if( lb == pairs_.end() || lb->key_ != key ) {
            return nullptr;
        } else {
            return &lb->value_;
        }
    }
    
    void reserve( size_t s ) {
        pairs_.reserve(s);
    }
        
    
private:
    struct ipair {
        K key_;
        V value_;
        
        ipair( const K &key, const V &value ) : key_(key), value_(value) {}
        
        
        inline bool operator<( const ipair &other ) const {
            return key_ < other.key_;
        }
        
        
    };
    
    
    std::vector<ipair> pairs_;
    bool sorted_;
};

class sequences {

public:
    //sequences() = delete;
    
    sequences( std::istream &is ) : pw_scoring_matrix_( 3, 0 ) {
        assert( is.good() );

//         std::cout << "here\n";
        ivy_mike::read_fasta( is, names_, seqs_ );

        std::for_each( seqs_.begin(), seqs_.end(), boost::bind( &sequences::normalize_seq, this, _1) );

        std::vector<std::vector<uint8_t> > qs_mapped;
        mapped_seqs_.reserve(names_.size() );

        // pre-map raw qs seqs to 'state numbers' (=scoring matrix rows/columns)
        for( auto it = seqs_.begin(); it != seqs_.end(); ++it)
        {
            mapped_seqs_.push_back(std::vector< uint8_t >());//(it->size()));
            mapped_seqs_.back().reserve(it->size());
           
            std::for_each( it->begin(), it->end(), scoring_matrix::valid_state_appender<std::vector< uint8_t > >(pw_scoring_matrix_, mapped_seqs_.back() ));

//             std::copy( mapped_seqs_.back().begin(), mapped_seqs_.back().end(), std::ostream_iterator<int>( std::cout, " " ));

            // the raw sequences stored in 'seqs_' are filtered for invalid states (e.g., gaps and characters not
            // present in the scoring matrix already. So the members of 'mapped_seqs_' must have the same length.

            assert( mapped_seqs_.back().size() == it->size() );


        }

        max_name_len_ = 0;
        name_to_index_.reserve(names_.size());
        for( size_t i = 0; i < names_.size(); ++i ) {
            name_to_index_.put_fast( names_[i], i );
            max_name_len_ = std::max( max_name_len_, names_[i].size() );
        }
        name_to_index_.sort();
        //std::sort( name_to_index_.begin(), name_to_index_.end() );



//        calc_dist_matrix( false );

    }

    const std::vector<sequence> &mapped_seqs() const {
        return mapped_seqs_;
    }
    const std::vector<sequence> &seqs() const {
        return seqs_;
    }


    const sequence &seq_at( size_t i ) const {
        return seqs_.at(i);
    }

    const sequence &mapped_seq_at( size_t i ) const {
        return mapped_seqs_.at(i);
    }

    
    const scoring_matrix &pw_scoring_matrix() const {
        return pw_scoring_matrix_;
    }

    const std::string &name_at( size_t i ) const {
        return names_.at(i);
    }

    
    size_t name_to_index( const std::string &name ) const {
        auto rp = name_to_index_.get(name);
        
        if( rp == nullptr ) {
            std::cerr << "name: " << name << "\n";
            throw std::runtime_error( "name not found" );
        } else {
            return *rp;
        }
        
//         auto it = std::lower_bound( name_to_index_.begin(), name_to_index_.end(), name_to_index_pair( name, size_t(-1) ));
//         
//         if( it == name_to_index_.end() || it->name_ != name ) {
//             std::cerr << "name: " << name << "\n";
//             throw std::runtime_error( "name not found" );
//         } else {
//             return it->index();
//         }
    }

    size_t clone_seq( size_t i, const std::string &name ) {
        std::vector<uint8_t> seq = seqs_.at(i);
        std::vector<uint8_t> mapped_seq = mapped_seqs_.at(i);

        ivy_mike::push_back_swap(seqs_, seq);
        ivy_mike::push_back_swap(mapped_seqs_, mapped_seq);

//        seqs_.push_back( seq );
//        mapped_seqs_.push_back( mapped_seq );
        names_.push_back(name);

//         name_to_index_pair nni( name, names_.size() - 1 );
//         name_to_index_.insert( std::lower_bound( name_to_index_.begin(), name_to_index_.end(), nni ), nni ); 
        
        name_to_index_.put( name, names_.size() - 1 );
        
        return seqs_.size() - 1;
    }

    size_t max_name_len() const {
        return max_name_len_;
    }
    
private:

    void normalize_seq( std::vector<uint8_t> &seq ) {
        std::vector<uint8_t> nseq;
      //  nseq.reserve(seq.size());

        for( std::vector<uint8_t>::iterator it = seq.begin(); it != seq.end(); ++it ) {
            uint8_t c = std::toupper(*it);

            if( pw_scoring_matrix_.state_valid(c)) {
                nseq.push_back(c);
            }
        }

        // shrink-to-fit into original vector 'seq'
       // std::vector<uint8_t> tmp( nseq.begin(), nseq.end() );
        seq.swap(nseq);

    }

    class name_to_index_pair {
    public:
        
        
        name_to_index_pair( const std::string &name, size_t idx ) : name_(name), index_(idx) {}
        
        bool operator<( const name_to_index_pair &other ) const {
            return name_ < other.name_;
        }
        
        inline size_t index() const {
            return index_;
        }
        
    private:
        std::string name_;
        size_t index_;
    };
    
    std::vector<std::vector<uint8_t> > seqs_;
    std::vector<std::vector<uint8_t> > mapped_seqs_;
    std::vector<std::string> names_;
    size_t max_name_len_;
    //std::vector<name_to_index_pair> name_to_index_;
    flat_map<std::string,size_t> name_to_index_;

    scoring_matrix pw_scoring_matrix_;


};


static void make_tip( lnode *n, const std::string &name ) {
    assert( n != 0 );
    assert( n->m_data != 0 );


    // check if the node cn be a valid tip: at least two back pointers must be null.
    size_t null_back = 0;
    if( n->back == 0 ) {
        ++null_back;
    }
    if( n->next->back == 0 ) {
        ++null_back;
    }
    if( n->next->next->back == 0 ) {
        ++null_back;
    }

    assert( null_back >= 2 );

    n->m_data->isTip = true;
    n->m_data->setTipName( name );
}

static bool has_node_label( lnode *n ) {
    assert( n != 0 );
    assert( n->m_data != 0 );
    return !n->m_data->nodeLabel.empty();
}

class tree_builder {
public:

    tree_builder( sequences * const seqs, addition_order * const order, ln_pool * const pool, lnode *destiny_tree )
     : seqs_(*seqs),
       order_(order),
       pool_(pool),
       destiny_tree_(destiny_tree),
       destiny_tree_pin_( destiny_tree_, *pool_ ),
       clones_pruned_(false)
    {
        // build the initial tree. This is mostly based on black magic.

        std::pair<size_t,size_t> first = order->first_pair();

        size_t seqa = first.first; // confusing, isn't it?
        size_t seqb = first.second;

    //    std::copy( seqs.seq_at( seqa ).begin(), seqs.seq_at( seqa ).end(), std::ostream_iterator<char>(std::cout));
    //    std::cout << "\n";
    //    std::copy( seqs.seq_at( seqb ).begin(), seqs.seq_at( seqb ).end(), std::ostream_iterator<char>(std::cout));
    //    std::cout << "\n";


        sequence aligned_a = seqs->seq_at(seqa);
        sequence aligned_b = seqs->seq_at(seqb);

        std::string name_clonea = seqs->name_at( seqa ) + "_clone";
        std::string name_cloneb = seqs->name_at( seqb ) + "_clone";

        cloned_names_.push_back(name_clonea);
        cloned_names_.push_back(name_cloneb);
        
        size_t seqa_clone = seqs->clone_seq( seqa, name_clonea );
        size_t seqb_clone = seqs->clone_seq( seqb, name_cloneb );

        used_seqs_.resize( seqs->seqs().size(), false );
        used_seqs_[seqa] = true;
        used_seqs_[seqa_clone] = true;
        used_seqs_[seqb] = true;
        used_seqs_[seqb_clone] = true;

        aligned_seqs_.resize( seqs->seqs().size() );

    //    std::copy( aligned_a.begin(), aligned_a.end(), std::ostream_iterator<char>(std::cout));
    //    std::cout << "\n";
    //    std::copy( aligned_b.begin(), aligned_b.end(), std::ostream_iterator<char>(std::cout));
    //    std::cout << "\n";


        lnode *nx = lnode::create( *pool );
        lnode *ny = lnode::create( *pool );
        tree_parser::twiddle_nodes(nx, ny, 1.0, "MOAL", 0 );


        lnode *na1 = lnode::create( *pool );
        lnode *na2 = lnode::create( *pool );

        lnode *nb1 = lnode::create( *pool );
        lnode *nb2 = lnode::create( *pool );

        make_tip( na1, seqs->name_at( seqa ));
        make_tip( na2, name_clonea);
        make_tip( nb1, seqs->name_at( seqb ));
        make_tip( nb2, name_cloneb);

        
        

        tree_parser::twiddle_nodes(na1, nx->next, 1.0, "I1", 0 );
        tree_parser::twiddle_nodes(na2, nx->next->next, 1.0, "I2", 0 );
        tree_parser::twiddle_nodes(nb1, ny->next, 1.0, "I3", 0 );
        tree_parser::twiddle_nodes(nb2, ny->next->next, 1.0, "I4", 0 );

        align_freeshift( seqs->pw_scoring_matrix(), aligned_a, aligned_b, -5, -3 );
        assert( aligned_a.size() == aligned_b.size() );

        
        aligned_seqs_[seqa_clone] = aligned_a;
        aligned_seqs_[seqa].swap( aligned_a );
        aligned_seqs_[seqb_clone] = aligned_b;
        aligned_seqs_[seqb].swap( aligned_b );

        
        
        tree_ = nx;
        
        {
            std::vector<lnode *> dt;
            
            ivy_mike::iterate_lnode( destiny_tree, ivy_mike::back_insert_ifer(dt, ivy_mike::is_tip));
            
            for( lnode *n : dt ) {
                //    destiny_tree_tips_.push_back
                
                assert( n->m_data->isTip );
                destiny_tree_tips_.put_fast( n->m_data->tipName, n );
            
            }
            
            destiny_tree_tips_.sort();
        }
        

    }
    double calc_gap_freq () {
        size_t ngaps = 0;
        size_t nres = 0;

        size_t idx = used_seqs_.find_first();
        
        //for( std::vector< std::vector< uint8_t > >::const_iterator it = seqs.begin(); it != seqs.end(); ++it ) {
        while( idx != used_seqs_.npos ) {
            const sequence &seq = aligned_seqs_.at(idx);
            assert( !seq.empty() );
            
            nres += seq.size();
            ngaps += std::count_if( seq.begin(), seq.end(), []( unsigned char c ) {return c == '-'; } ); // TODO: this should depend on seq_model
            
            idx = used_seqs_.find_next(idx);
        }

        double rgap = double(ngaps) / nres;
        std::cout << "gap rate: " << ngaps << " " << nres << "\n";
        std::cout << "gap rate: " << rgap << "\n";
        return rgap;
    }
    
    void align_ref_seqs( std::ostream &os, const std::vector<uint8_t> &tb ) {
        auto idx = used_seqs_.find_first();
        
        while( idx != used_seqs_.npos ) {
            const sequence &seq = aligned_seqs_.at(idx);
            auto ali = gapstream_to_alignment( tb, seq, '-', true );
            
            std::copy( ali.begin(), ali.end(), std::ostream_iterator<char>(std::cout) );
            std::cout << "\n";
            
            aligned_seqs_.at(idx) = std::move(ali);
            
            idx = used_seqs_.find_next(idx);
        }
    }
    void init_tree_sequences() {
        // initialize the tip nodes with the aligned sequence data
        
        apply_lnode( tree_, [&]( lnode *n ) {
            if( n->m_data->isTip ) {
                size_t idx = seqs_.name_to_index( n->m_data->tipName );
                const sequence &seq = aligned_seqs_.at(idx);
                
                assert( !seq.empty() );
                std::cout << "init: " << n->m_data->tipName << "\n";
                n->m_data->get_as<my_adata>()->init_sequence( seq );
                n->m_data->get_as<my_adata>()->update_pvec();
            }
        } );
//                 na1->m_data->get_as<my_adata>()->init_sequence(aligned_a);
//         na2->m_data->get_as<my_adata>()->init_sequence(aligned_a);
//         nb1->m_data->get_as<my_adata>()->init_sequence(aligned_b);
//         na2->m_data->get_as<my_adata>()->init_sequence(aligned_b);



    }


    void print_matrix( const ublas::matrix<double> &m ) {
       // for( ; first != last)
    }


    
    bool insertion_step_destiny( const size_t cand_id, const std::vector<std::string> &insertion_pos ) {
        
                
        
        std::vector<ublas::matrix<double> > pvecs;
        write_ali_and_tree_for_raxml();

        
        // generate anc state pvecs using external raxml, and replace the current 'main-tree' with the one 
        // from raxml.
        tree_ = generate_marginal_ancestral_state_pvecs( *pool_, "sa_tree", "sa_ali", &pvecs );
        init_tree_sequences();
        
        double gap_freq = calc_gap_freq();
        
//         size_t cand_id = order_->find_next_candidate();
        std::cout << "cand_id: " << cand_id << "\n";
        
        if( cand_id == size_t(-1) ) { 
            return false;
            
        }

        
        const auto &cand_seq = seqs_.seq_at(cand_id);
        const auto &cand_mapped_seq = seqs_.mapped_seq_at(cand_id);
        
        lnode *insertion_edge = nullptr;
        
        {        
            // find insertion position of cand_seq, according to destiny tree
            
            
            std::vector<lnode*> edges;
            std::vector<boost::dynamic_bitset<> > splits;
            std::vector<lnode *> sorted_tips;
            ivy_mike::get_all_splits_by_node( tree_,  edges, splits, sorted_tips );
            
            std::vector<std::string> tip_names;
            
            // transform the tips in sorted_tips into (sorted) list of taxon names
            std::transform( sorted_tips.begin(), sorted_tips.end(), std::back_inserter( tip_names ), [](lnode *n){ return n->m_data->tipName; } );
            std::sort( tip_names.begin(), tip_names.end() );
            
            
            if( !clones_pruned_ ) {
                // delete the 'cloned' nodes of the two initial taxa.
                delete_cloned_nodes( splits, sorted_tips );
            }
            
            const auto &cand_name = seqs_.name_at(cand_id);
#if 0                   
            // look up the tip in the destiny tree, which corresponds to cand_name
            lnode * const * dnode_ptr = destiny_tree_tips_.get( cand_name );
            assert( dnode_ptr != 0 );
            lnode *dnode = *dnode_ptr;
            
            // the node to be (temporarily) pruned from the destiny tree
            lnode *pnode = dnode->back;
            
            assert( pnode != 0 );
            assert( !pnode->m_data->isTip );
            
            ivy_mike::tree_parser_ms::prune_with_rollback pwr( pnode );
            
            
            // generate split set from the destiny tree
            std::vector<std::string> split_set = ivy_mike::get_split_set_by_edge(pwr.get_save_node());
            
           
     
#if 1
            // remove taxon names from the (destiny tree) split set which are not (yet) 
            // contained in the constructed tree
            
            // functor that returns true if a tip name is not contained in sorted_tips
            auto is_not_in_tree = [&](const std::string &s){ return !std::binary_search( tip_names.begin(), tip_names.end(), s); };
            
            // remove all tip names from split_set (filter by is_not_in_tree)
   
            split_set.erase( std::remove_if( split_set.begin(), split_set.end(), is_not_in_tree ), split_set.end() );
            
//             std::sort( split_set.begin(), spit_set.end() );
#endif

            std::cout << "cand name: " << cand_name << "\n";
            
            std::cout << "destiny tree split: ";
            std::copy( split_set.begin(), split_set.end(), std::ostream_iterator<std::string>(std::cout, " " ));
            std::cout << "\n";
            
            // find the split in the constructed three that corresponds to the split from the destiny tree
            std::sort( split_set.begin(), split_set.end() );
#else
            auto split_set = insertion_pos;
            std::sort( split_set.begin(), split_set.end() );
#endif
            std::cout << "cand name: " << cand_name << "\n";
            
            std::cout << "destiny tree split: ";
            std::copy( split_set.begin(), split_set.end(), std::ostream_iterator<std::string>(std::cout, " " ));
            std::cout << "\n";
            
            
            auto bs_size = splits.front().size();
            boost::dynamic_bitset<> bsplit_set;
            
            for( size_t i = 0; i < bs_size; ++i ) {
                std::string tip_name = sorted_tips.at(i)->m_data->tipName;
                std::cout << "tip name: " << tip_name << "\n";
                bool do_set = std::binary_search( split_set.begin(), split_set.end(), tip_name );
                bsplit_set.push_back(do_set);
            }
            auto bsplit_set_comp = bsplit_set;
            bsplit_set_comp.flip();
            
            for( size_t i = 0; i < bsplit_set.size(); ++i ) {
                std::cout << bsplit_set[i] << " ";
            }
            
            
            std::cout << "\n";
            
            // now (hopefully) one of the splits in 'splits' should be equal to bsplit_set or its complement
            
            auto sit = splits.begin();
            for( auto e = splits.end(); sit != e; ++sit ) {
//                 for( size_t i = 0; i < bsplit_set.size(); ++i ) {
//                     std::cout << (*sit)[i] << " ";
//                 }
//                 
//                 
//                 std::cout << "\n";
                
                if( bsplit_set == *sit || bsplit_set_comp == *sit ) {
                    break;
                }
            }
            
            assert( sit != splits.end() );
            assert( splits.size() == edges.size() );
            insertion_edge = edges.at( std::distance( splits.begin(), sit ));
        }
        
        
        probgap_model gpm(gap_freq);
        ivy_mike::stupid_ptr_guard<probgap_model> spg( pvec_pgap::pgap_model, &gpm );
        
        
        std::cout << "pvecs: " << pvecs.size() << "\n";

//         std::cout << n->backLabel << " " << n->next->backLabel << " " << n->next->next->backLabel << "\n";


//        std::deque<rooted_bifurcation<lnode> > to;
//        rooted_traveral_order_rec(n->next, to, false );
//
//        for( std::deque<rooted_bifurcation<lnode> >::iterator it = to.begin(); it != to.end(); ++it ) {
//            std::cout << *it << "\n";
//        }


//         for( auto it = pvecs.begin(); it != pvecs.end(); ++it ) {
//             std::cout << "vec:\n";
//             
//             
//             
//             for( auto it1 = it->begin2(); it1 != it->end2(); ++it1 ) {
//                 std::transform( it1.begin(), it1.end(), std::ostream_iterator<double>( std::cout, "\t" ), [](double x) {return x; /*std::max(1.0,-log(x));*/});
//                 std::cout << "\n";
//             }
//             
//         }
        
        
        std::vector<lnode *> labelled_nodes;
        iterate_lnode(tree_, back_insert_ifer( labelled_nodes, has_node_label ));

        
        
        std::cout << "num labelled: " << labelled_nodes.size() << "\n";

        std::sort( labelled_nodes.begin(), labelled_nodes.end(), [](lnode *n1, lnode *n2){return n1->m_data->nodeLabel < n2->m_data->nodeLabel;} ); 
        
        
        lnode *virtual_root = lnode::create( *pool_ );
        
        
//         bool incremental = false;
        
        double best_score = -std::numeric_limits<double>::infinity();
        std::vector<uint8_t> best_tb;
        lnode *best_np = nullptr;
        
        
        //for( lnode *np : labelled_nodes ) {
        {
            lnode *np = nullptr;
            
            if( has_node_label(insertion_edge)) {
                np = insertion_edge;
            } else if( has_node_label( insertion_edge->back ) ) {
                np = insertion_edge->back;
            }
            
            assert( np != nullptr );
            
             
            // splice virtual root into current insertion edge.
            // NOTE: splice_with_rollback will automatically undo the insertion on scope-exit
            ivy_mike::tree_parser_ms::splice_with_rollback swr( np, virtual_root );
            
            
            std::deque<rooted_bifurcation<lnode> > rto;
            rooted_traveral_order_rec( virtual_root, rto, false );
//             incremental = true;
            for( auto it = rto.begin(); it != rto.end(); ++it ) {
                my_adata *p = it->parent->m_data->get_as<my_adata>();
                my_adata *c1 = it->child1->m_data->get_as<my_adata>();
                my_adata *c2 = it->child2->m_data->get_as<my_adata>();
//                  std::cout << "newview: " << *it << " " << p << " " << it->parent << "\n";
                //         std::cout << "tip case: " << (*it) << "\n";
                 
                auto z1 = it->child1->backLen;
                auto z2 = it->child2->backLen;
                 
//                 z1 = 0.0001;
//                 z2 = 0.0001;
                
                pvec_pgap::newview(p->pvec(), c1->pvec(), c2->pvec(), z1, z2, it->tc);
                
            }
            
            const bool verbose_scoring = false;
            if( verbose_scoring ) {
                std::cout << "vr: " << *virtual_root->m_data << " " << virtual_root << "\n";
            }
            const pvec_pgap &rpp = virtual_root->m_data->get_as<my_adata>()->pvec();
            const boost::numeric::ublas::matrix< double > &pm = rpp.get_gap_prob();
            const auto &anc_gap = virtual_root->m_data->get_as<my_adata>()->calculate_anc_gap_probs();
            
            if( verbose_scoring ) {
                std::cerr << pm.size1() << " " << pm.size2() << "\n";
            }
//             for( auto it1 = pm.begin2(); it1 != pm.end2(); ++it1 ) {
//                 std::transform( it1.begin(), it1.end(), std::ostream_iterator<double>( std::cerr, "\t" ), [](double x) {return x; /*std::max(1.0,-log(x));*/});
//                 std::cerr << "\n";
//             }
//             std::cerr << "===============\n";
//             for( auto it1 = anc_gap.begin2(); it1 != anc_gap.end2(); ++it1 ) {
//                 std::transform( it1.begin(), it1.end(), std::ostream_iterator<double>( std::cerr, "\t" ), [](double x) {return x; /*std::max(1.0,-log(x));*/});
//                 std::cerr << "\n";
//             }
//             std::cerr << "node label: " << np->m_data->nodeLabel << "\n";
            size_t node_label = size_t(-1);
            {
                std::stringstream ss( np->m_data->nodeLabel );
                ss >> node_label;
            }
            assert( node_label != size_t(-1) );
            
            
            if( verbose_scoring ) {
                std::cerr << "node label: " << node_label << "\n";
            }
            auto const & anc_state = pvecs.at( node_label );
            
            
            boost::array<double,4> bg_state{0.25, 0.25, 0.25, 0.25};
            log_odds_viterbi lov(anc_state, anc_gap, bg_state );
            auto score = lov.align(cand_mapped_seq);
            
            std::cerr << "cand:\n";
            std::copy( cand_seq.begin(), cand_seq.end(), std::ostream_iterator<char>(std::cerr));
            std::cerr << "\n";
            
            if( verbose_scoring ) {
                std::cerr << "score: " << score << "\n";
            }
            
            
            auto tb = lov.traceback();
            auto qs_ali = gapstream_to_alignment(tb, cand_seq, '-', false );
            std::copy( tb.begin(), tb.end(), std::ostream_iterator<int>( std::cerr, " " ) );
            std::cerr << std::endl;
//            
//             dump_ref_seqs( std::cout, tb );
//             
//             std::copy( qs_ali.begin(), qs_ali.end(), std::ostream_iterator<char>( std::cout ) );
//             std::cout << "\n";
            
            if( score > best_score ) {
                best_score = score;
                best_tb = lov.traceback();
                best_np = np;
            }
            
//             std::cout << np->m_data->nodeLabel << " " << np->m_data->isTip << "\n";
//             
//             
//             tree_parser::print_newick( np, std::cout, false );
//             std::cout << "\n";
        }
        
        {

            assert( best_np != nullptr );
            
            ivy_mike::tree_parser_ms::splice_with_rollback swr( best_np, virtual_root );
            swr.commit();
            virtual_root->back = lnode::create(*pool_);
            virtual_root->back->m_data->setTipName( seqs_.name_at(cand_id));
            virtual_root->back->m_data->isTip = true;
            
            auto qs_ali = gapstream_to_alignment(best_tb, cand_seq, '-', false );
//             std::copy( tb.begin(), tb.end(), std::ostream_iterator<int>( std::cout, " " ) );
//             std::cout << "\n";
           
            align_ref_seqs( std::cout, best_tb );
            
//             std::copy( qs_ali.begin(), qs_ali.end(), std::ostream_iterator<char>( std::cout ) );
//             std::cout << "\n";
            
            aligned_seqs_.at(cand_id) = std::move(qs_ali);
            
            
            
        }
        
        // officially add cand_id to the set of completed sequences.
        used_seqs_[cand_id] = true;
        
        // write aligned seqeuences
        {
            if( !inc_log_.good() ) {
                inc_log_.open("inc_log" );
                
            }
            
            inc_log_ << cand_id << "\n\n";
            size_t idx = used_seqs_.find_first();
            
            while( idx != used_seqs_.npos ) {
                
                const sequence &seq = aligned_seqs_.at(idx);
                
                inc_log_ << seqs_.name_at(idx) << " ";
                std::copy( seq.begin(), seq.end(), std::ostream_iterator<char>(inc_log_) );                
                inc_log_ << "\n";
                
                idx = used_seqs_.find_next(idx);
            }
            inc_log_.flush();
            
        }
        

        //pool_->clear();
        //pool_->mark(tree_);
        //pool_->sweep();
        
        return true;
        
    }
    
    bool insertion_step() {

        std::vector<ublas::matrix<double> > pvecs;
        write_ali_and_tree_for_raxml();

        
        // generate anc state pvecs using external raxml, and replace the current 'main-tree' with the one 
        // from raxml.
        tree_ = generate_marginal_ancestral_state_pvecs( *pool_, "sa_tree", "sa_ali", &pvecs );
        init_tree_sequences();
        
        double gap_freq = calc_gap_freq();
        
        size_t cand_id = order_->find_next_candidate();
        std::cout << "cand_id: " << cand_id << "\n";
        
        if( cand_id == size_t(-1) ) { 
            return false;
            
        }

        
        const auto &cand_seq = seqs_.seq_at(cand_id);
        const auto &cand_mapped_seq = seqs_.mapped_seq_at(cand_id);
        
        
        probgap_model gpm(gap_freq);
        ivy_mike::stupid_ptr_guard<probgap_model> spg( pvec_pgap::pgap_model, &gpm );
        
        
        std::cout << "pvecs: " << pvecs.size() << "\n";

//         std::cout << n->backLabel << " " << n->next->backLabel << " " << n->next->next->backLabel << "\n";


//        std::deque<rooted_bifurcation<lnode> > to;
//        rooted_traveral_order_rec(n->next, to, false );
//
//        for( std::deque<rooted_bifurcation<lnode> >::iterator it = to.begin(); it != to.end(); ++it ) {
//            std::cout << *it << "\n";
//        }


//         for( auto it = pvecs.begin(); it != pvecs.end(); ++it ) {
//             std::cout << "vec:\n";
//             
//             
//             
//             for( auto it1 = it->begin2(); it1 != it->end2(); ++it1 ) {
//                 std::transform( it1.begin(), it1.end(), std::ostream_iterator<double>( std::cout, "\t" ), [](double x) {return x; /*std::max(1.0,-log(x));*/});
//                 std::cout << "\n";
//             }
//             
//         }
        
        
        std::vector<lnode *> labelled_nodes;
        iterate_lnode(tree_, back_insert_ifer( labelled_nodes, has_node_label ));

        
        
        std::cout << "num labelled: " << labelled_nodes.size() << "\n";

        std::sort( labelled_nodes.begin(), labelled_nodes.end(), [](lnode *n1, lnode *n2){return n1->m_data->nodeLabel < n2->m_data->nodeLabel;} ); 
        
        
        lnode *virtual_root = lnode::create( *pool_ );
        
        
        bool incremental = false;
        
        double best_score = -std::numeric_limits<double>::infinity();
        std::vector<uint8_t> best_tb;
        lnode *best_np = nullptr;
        
        
        for( lnode *np : labelled_nodes ) {

            // splice virtual root into current insertion edge.
            // NOTE: splice_with_rollback will automatically undo the insertion on scope-exit
            ivy_mike::tree_parser_ms::splice_with_rollback swr( np, virtual_root );
            
            
            std::deque<rooted_bifurcation<lnode>> rto;
            rooted_traveral_order_rec( virtual_root, rto, incremental );
            incremental = true;
            for( auto it = rto.begin(); it != rto.end(); ++it ) {
                my_adata *p = it->parent->m_data->get_as<my_adata>();
                my_adata *c1 = it->child1->m_data->get_as<my_adata>();
                my_adata *c2 = it->child2->m_data->get_as<my_adata>();
//                  std::cout << "newview: " << *it << " " << p << " " << it->parent << "\n";
                //         std::cout << "tip case: " << (*it) << "\n";
                 
                auto z1 = it->child1->backLen;
                auto z2 = it->child2->backLen;
                 
//                 z1 = 0.0001;
//                 z2 = 0.0001;
                
                pvec_pgap::newview(p->pvec(), c1->pvec(), c2->pvec(), z1, z2, it->tc);
                
            }
            
            const bool verbose_scoring = false;
            if( verbose_scoring ) {
                std::cout << "vr: " << *virtual_root->m_data << " " << virtual_root << "\n";
            }
            const pvec_pgap &rpp = virtual_root->m_data->get_as<my_adata>()->pvec();
            const boost::numeric::ublas::matrix< double > &pm = rpp.get_gap_prob();
            const auto &anc_gap = virtual_root->m_data->get_as<my_adata>()->calculate_anc_gap_probs();
            
            if( verbose_scoring ) {
                std::cout << pm.size1() << " " << pm.size2() << "\n";
            }
//             for( auto it1 = pm.begin2(); it1 != pm.end2(); ++it1 ) {
//                 std::transform( it1.begin(), it1.end(), std::ostream_iterator<double>( std::cout, "\t" ), [](double x) {return x; /*std::max(1.0,-log(x));*/});
//                 std::cout << "\n";
//             }
//             std::cout << "===============\n";
//             for( auto it1 = anc_gap.begin2(); it1 != anc_gap.end2(); ++it1 ) {
//                 std::transform( it1.begin(), it1.end(), std::ostream_iterator<double>( std::cout, "\t" ), [](double x) {return x; /*std::max(1.0,-log(x));*/});
//                 std::cout << "\n";
//             }
            
            size_t node_label = size_t(-1);
            {
                std::stringstream ss( np->m_data->nodeLabel );
                ss >> node_label;
            }
            assert( node_label != size_t(-1) );
            
            
            if( verbose_scoring ) {
                std::cout << "node label: " << node_label << "\n";
            }
            auto const & anc_state = pvecs.at( node_label );
            
            
            boost::array<double,4> bg_state{0.25, 0.25, 0.25, 0.25};
            log_odds_viterbi lov(anc_state, anc_gap, bg_state );
            auto score = lov.align(cand_mapped_seq);
            
//             std::cout << "cand:\n";
//             std::copy( cand_seq.begin(), cand_seq.end(), std::ostream_iterator<char>(std::cout));
//             std::cout << "\n";
            
            if( verbose_scoring ) {
                std::cout << "score: " << score << "\n";
            }
            
            
//             auto tb = lov.traceback();
//             auto qs_ali = gapstream_to_alignment(tb, cand_seq, '-', false );
//             std::copy( tb.begin(), tb.end(), std::ostream_iterator<int>( std::cout, " " ) );
//             std::cout << "\n";
//            
//             dump_ref_seqs( std::cout, tb );
//             
//             std::copy( qs_ali.begin(), qs_ali.end(), std::ostream_iterator<char>( std::cout ) );
//             std::cout << "\n";
            
            if( score > best_score ) {
                best_score = score;
                best_tb = lov.traceback();
                best_np = np;
            }
            
//             std::cout << np->m_data->nodeLabel << " " << np->m_data->isTip << "\n";
//             
//             
//             tree_parser::print_newick( np, std::cout, false );
//             std::cout << "\n";
        }
        
        {

            assert( best_np != nullptr );
            
            ivy_mike::tree_parser_ms::splice_with_rollback swr( best_np, virtual_root );
            swr.commit();
            virtual_root->back = lnode::create(*pool_);
            virtual_root->back->m_data->setTipName( seqs_.name_at(cand_id));
            virtual_root->back->m_data->isTip = true;
            
            auto qs_ali = gapstream_to_alignment(best_tb, cand_seq, '-', false );
//             std::copy( tb.begin(), tb.end(), std::ostream_iterator<int>( std::cout, " " ) );
//             std::cout << "\n";
           
            align_ref_seqs( std::cout, best_tb );
            
//             std::copy( qs_ali.begin(), qs_ali.end(), std::ostream_iterator<char>( std::cout ) );
//             std::cout << "\n";
            
            aligned_seqs_.at(cand_id) = std::move(qs_ali);
            
            
            
        }
        
        // officially add cand_id to the set of completed sequences.
        used_seqs_[cand_id] = true;
        
        // write aligned seqeuences
        {
            if( !inc_log_.good() ) {
                inc_log_.open("inc_log" );
                
            }
            
            inc_log_ << cand_id << "\n\n";
            size_t idx = used_seqs_.find_first();
            
            while( idx != used_seqs_.npos ) {
                
                const sequence &seq = aligned_seqs_.at(idx);
                
            
                std::copy( seq.begin(), seq.end(), std::ostream_iterator<char>(inc_log_) );                
                inc_log_ << "\n";
                
                idx = used_seqs_.find_next(idx);
            }
            
            
        }
        

        pool_->clear();
        pool_->mark(tree_);
        pool_->sweep();
        
        return true;
    }

    void write_ali_and_tree( const char *tree_name, const char *ali_name, bool pad = true ) {
        
        std::cout << "write tree: " << tree_name << "\n";
        {
            std::ofstream os( tree_name );
            tree_parser::print_newick( tree_, os );
        }



        std::ofstream os( ali_name );


        size_t pos = used_seqs_.find_first();

        const bool write_clones = !pad && !clones_pruned_;
        if( !write_clones ) {
            os << (used_seqs_.count() - cloned_names_.size()) << " " << aligned_seqs_.at(pos).size() << "\n";
        } else {
            os << used_seqs_.count() << " " << aligned_seqs_.at(pos).size() << "\n";
        }

        
        std::sort( cloned_names_.begin(), cloned_names_.end() ); // most likely they are already sorted from a previous call of prune_cloned_nodes
        
        std::copy( cloned_names_.begin(), cloned_names_.end(), std::ostream_iterator<std::string>(std::cout, "\n"));
        while( pos != used_seqs_.npos ) {
            const auto &name = seqs_.name_at(pos);
            
            
//             std::cerr << "name: " << name << " " << std::binary_search( cloned_names_.begin(), cloned_names_.begin(), name ) << "\n";
            // FIXME: kind of hack: don't print out cloned seqs when pad is set (=on final printout)
            
            const bool is_clone = std::binary_search( cloned_names_.begin(), cloned_names_.end(), name );
            
            
            
           
               // smart-ass: material implication: is_clone -> write_clones
               if( !is_clone || write_clones ) {
                assert( pos < aligned_seqs_.size() );
                
                os << name;
                
                // pad column
                if( pad ) {
                    const size_t max_len = seqs_.max_name_len();
                    //                 std::cout << "name: " << name << "\n";
                    assert( name.size() <= max_len );
                    
                    
                    for( size_t i = 0; i < max_len - name.size() + 1; ++i ) {
                        os << " ";
                    }
                } else {
                    os << " ";
                }
                std::copy( aligned_seqs_[pos].begin(), aligned_seqs_[pos].end(), std::ostream_iterator<char>(os));
                os << "\n";
            }
            pos = used_seqs_.find_next(pos);
        }

    }
    
    void prune_cloned_nodes() {
        std::vector<lnode *>rn;
        std::sort( cloned_names_.begin(), cloned_names_.end() );
        apply_lnode( tree_, [&](lnode *n) {
        
            if( n->m_data->isTip ) {
                // TODO: lol, I guess linear search would be better for the 2 elements in cloned_names_...
                if( std::binary_search( cloned_names_.begin(), cloned_names_.end(), n->m_data->tipName ) ) {
                    rn.push_back( n );
                }
            }
        });
        
        std::cout << "rn: " << rn.size() << "\n";
        
        for( lnode *n : rn ) {
            tree_parser::prune_with_rollback pwr(n->back);
            pwr.commit();
          
            tree_ = pwr.get_save_node();
        }
        
        clones_pruned_ = true;
        
        
    }
    void delete_cloned_nodes(std::vector< boost::dynamic_bitset<> > &splits, std::vector< lnode* > &sorted_tips) {
        ivy_mike::flat_map<std::string,size_t> tip_name_map;
        
        for( size_t i = 0; i < sorted_tips.size(); ++i ) {
            
            
            tip_name_map.put_fast( sorted_tips[i]->m_data->tipName, i );
        }
        tip_name_map.sort();
        
        
        std::vector<size_t> delete_idx;
        
        for( const auto &cn : cloned_names_ ) {
            const size_t *idx = tip_name_map.get(cn);
            
            assert( idx != 0 );
            
            delete_idx.push_back(*idx);
        }
        
        std::sort( delete_idx.begin(), delete_idx.end() );
        
        
        for( auto &split : splits ) {
            boost::dynamic_bitset<> out_split;
            
            size_t l = split.size();
            
            for( size_t i = 0; i < l; ++i ) {
                if( !binary_search( delete_idx.begin(), delete_idx.end(), i )) {
                    out_split.push_back( split[i] );
                }
            }
            split.swap( out_split );
            
        }
        
        for( size_t i = 0; i < delete_idx.size(); ++i ) {
            size_t di = delete_idx[i];
            std::cerr << "di: " << di << std::endl;
            auto dit = sorted_tips.begin() + di - i;
            
            sorted_tips.erase(dit); 
        }
        
        
        
    }
    
private:

    void write_ali_and_tree_for_raxml() {
        write_ali_and_tree( "sa_tree", "sa_ali", false );
        
    }
//         {
//             std::ofstream os( "sa_tree" );
//             tree_parser::print_newick( tree_, os );
//         }
// 
// 
// 
//         std::ofstream os( "sa_ali" );
// 
// 
//         size_t pos = used_seqs_.find_first();
// 
//         os << used_seqs_.count() << " " << aligned_seqs_.at(pos).size() << "\n";
// 
//         while( pos != used_seqs_.npos ){
//             assert( pos < aligned_seqs_.size() );
// 
//             os << seqs_.name_at(pos) << " ";
//             std::copy( aligned_seqs_[pos].begin(), aligned_seqs_[pos].end(), std::ostream_iterator<char>(os));
//             os << "\n";
// 
//             pos = used_seqs_.find_next(pos);
//         }
// 
//     }



    const sequences &seqs_;
    addition_order * const order_;
    ln_pool * const pool_;
    boost::dynamic_bitset<> used_seqs_;
    lnode * tree_;
    lnode * destiny_tree_;
    
    ivy_mike::tree_parser_ms::ln_pool_pin destiny_tree_pin_;
    ivy_mike::flat_map<std::string, lnode *> destiny_tree_tips_;
    
    std::vector<sequence> aligned_seqs_;

    std::ofstream inc_log_;
    std::vector<std::string> cloned_names_;
    bool clones_pruned_;
};

//void insertion_loop( sequences *seqs, addition_order *order, ln_pool * const pool ) {
//
//    std::pair<size_t,size_t> first = order->first_pair();
//
//    size_t seqa = first.first; // confusing, isn't it?
//    size_t seqb = first.second;
//
////    std::copy( seqs.seq_at( seqa ).begin(), seqs.seq_at( seqa ).end(), std::ostream_iterator<char>(std::cout));
////    std::cout << "\n";
////    std::copy( seqs.seq_at( seqb ).begin(), seqs.seq_at( seqb ).end(), std::ostream_iterator<char>(std::cout));
////    std::cout << "\n";
//
//
//    sequence aligned_a = seqs->seq_at(seqa);
//    sequence aligned_b = seqs->seq_at(seqb);
//
//    std::string name_clonea = seqs->name_at( seqa ) + "_clone";
//    std::string name_cloneb = seqs->name_at( seqb ) + "_clone";
//
//    size_t seqa_clone = seqs->clone_seq( seqa, name_clonea );
//    size_t seqb_clone = seqs->clone_seq( seqb, name_cloneb );
//
//
//
//
////    std::copy( aligned_a.begin(), aligned_a.end(), std::ostream_iterator<char>(std::cout));
////    std::cout << "\n";
////    std::copy( aligned_b.begin(), aligned_b.end(), std::ostream_iterator<char>(std::cout));
////    std::cout << "\n";
//
//
//    lnode *nx = lnode::create( *pool );
//    lnode *ny = lnode::create( *pool );
//    tree_parser::twiddle_nodes(nx, ny, 1.0, "MOAL", 0 );
//
//
//    lnode *na1 = lnode::create( *pool );
//    lnode *na2 = lnode::create( *pool );
//
//    lnode *nb1 = lnode::create( *pool );
//    lnode *nb2 = lnode::create( *pool );
//
//    make_tip( na1, seqs->name_at( seqa ));
//    make_tip( na2, name_clonea);
//    make_tip( nb1, seqs->name_at( seqb ));
//    make_tip( nb2, name_cloneb);
//
//
//    tree_parser::twiddle_nodes(na1, nx->next, 1.0, "I1", 0 );
//    tree_parser::twiddle_nodes(na2, nx->next->next, 1.0, "I2", 0 );
//    tree_parser::twiddle_nodes(nb1, ny->next, 1.0, "I3", 0 );
//    tree_parser::twiddle_nodes(nb2, ny->next->next, 1.0, "I4", 0 );
//
//    {
//        std::ofstream os( "sa_tree" );
//        tree_parser::print_newick( nx, os );
//    }
//
//    align_freeshift( seqs->pw_scoring_matrix(), aligned_a, aligned_b, -5, -3 );
//    assert( aligned_a.size() == aligned_b.size() );
//    {
//        std::ofstream os( "sa_ali" );
//        os << "4 " << aligned_a.size() << "\n";
//
//
//        os << seqs->name_at( seqa ) << " ";
//        std::copy( aligned_a.begin(), aligned_a.end(), std::ostream_iterator<char>(os));
//        os << "\n";
//        os << seqs->name_at( seqa_clone ) << " ";
//        std::copy( aligned_a.begin(), aligned_a.end(), std::ostream_iterator<char>(os));
//        os << "\n";
//
//
//
//
//        os << seqs->name_at( seqb ) << " ";
//        std::copy( aligned_b.begin(), aligned_b.end(), std::ostream_iterator<char>(os));
//        os << "\n";
//        os << seqs->name_at( seqb_clone ) << " ";
//        std::copy( aligned_b.begin(), aligned_b.end(), std::ostream_iterator<char>(os));
//        os << "\n";
//
//
//    }
//
//
////    size_t next;
////    while( (next = order->find_next_candidate()) != size_t(-1)) {
////        std::cout << "next: " << next << "\n";
////        std::copy( seqs->seq_at( next ).begin(), seqs->seq_at( next ).end(), std::ostream_iterator<char>(std::cout));
////        std::cout << "\n";
////
////    }
//}

int main( int argc, char *argv[] ) {
//    {
//        size_t n = 1024 * 1024 * 1024;
//        boost::dynamic_bitset<> bs( n, false );
//
//        for( size_t i = 0; i < n; ++i ) {
//            if( std::rand() < RAND_MAX / 2 ) {
//                bs.set( i, true );
//            }
//        }
//
//
//        ivy_mike::timer t1;
//        size_t x = 0;
//        size_t pos = bs.find_first();
//        while( pos != bs.npos) {
//            x += pos;
//            pos = bs.find_next(pos);
//        }
//
//
//        std::cout << "elapsed: " << t1.elapsed() << "\n";
//
//        return 0;
//    }



    ivy_mike::getopt::parser igp;
    std::string opt_seq_file;

	
    int num_cores = std::thread::hardware_concurrency();

    int opt_num_ali_threads;
    int opt_num_nv_threads;
    bool opt_load_scores;

    igp.add_opt('h', false );
    igp.add_opt('f', ivy_mike::getopt::value<std::string>(opt_seq_file) );
    igp.add_opt('j', ivy_mike::getopt::value<int>(opt_num_ali_threads).set_default(num_cores) );
    igp.add_opt('k', ivy_mike::getopt::value<int>(opt_num_nv_threads).set_default(1) );
    igp.add_opt('l', ivy_mike::getopt::value<bool>(opt_load_scores, true).set_default(false) );
    bool ret = igp.parse(argc, argv);

    if( igp.opt_count('h') != 0 || !ret ) {
        std::cout <<
        "  -h        print help message\n";
        return 0;

    }


    if( igp.opt_count('f') != 1 ) {
        std::cerr << "missing option -f\n";
#ifndef WIN32 // hack. make it easier to start inside visual studio
        return 0;
#endif
        opt_seq_file = "test_218/218.fa";
    }

    const char *filename = opt_seq_file.c_str();



    std::map<std::string, std::vector<uint8_t> >out_msa1;
//    std::shared_ptr<ln_pool> pool;//(new ln_pool(std::auto_ptr<node_data_factory>(new my_fact()) ));

    ln_pool pool( ln_pool::fact_ptr_type(new my_fact<my_adata>));



    std::ifstream sis( filename );
    sequences seqs( sis );

    addition_order order( seqs.pw_scoring_matrix(), seqs.mapped_seqs() );

    
        
    lnode *destiny_tree;
    
    {
        const char * destiny_name = "RAxML_randomTree.rnd_150_ss";
        
        ivy_mike::tree_parser_ms::parser p( destiny_name, pool );
        
        destiny_tree = p.parse();
        
    }
    
    assert( destiny_tree != 0 );
    
    std::vector<size_t> insert_id_order;
    std::vector<std::string> insert_order;
    
    while( true ) {
        size_t cand = order.find_next_candidate();
        
        
        
        if( cand == size_t(-1) ) {
            break;
        }
        
        insert_id_order.push_back(cand);
        insert_order.push_back( seqs.name_at(cand) );
    }
    
    
    tree_builder builder( &seqs, &order, &pool, destiny_tree );
    tree_deconstructor td( destiny_tree, insert_order );
    
    auto cand_it = insert_id_order.begin();
    size_t step = 0;
    
    while( !td.empty() ) {
        lnode *t = td.get_save_node();
        
        assert( cand_it != insert_id_order.end() );
        
        
        // t = ivy_mike::tree_parser_ms::next_non_tip(t);
        auto ss = ivy_mike::get_split_set_by_edge(t);
        auto cand_id = *cand_it++;
        
        builder.insertion_step_destiny(cand_id, ss);
        
        std::cout << td.size() << " ";
        std::copy( ss.begin(), ss.end(), std::ostream_iterator<std::string>(std::cout, " " ) );
        std::cout << "\n";
            
        {
            std::stringstream ss;
            ss << std::setfill('0') << std::setw(3) << std::right << step;
            
            std::string tree_name( "sa_tree_inc_" );
            tree_name += ss.str();
            
            std::string ali_name( "sa_ali_inc_" );
            ali_name += ss.str();
            
            builder.write_ali_and_tree( tree_name.c_str(), ali_name.c_str() );
            ++step;
        }
        
        
        td.pop();
    }
    
    
//     return 0;

    
#if 0
//    insertion_loop( &seqs, &order, &pool );
    tree_builder builder( &seqs, &order, &pool, destiny_tree );


    
    size_t step = 0;
    while(true) {
        bool valid = builder.insertion_step_destiny();
        
        {
            std::stringstream ss;
            ss << std::setfill('0') << std::setw(3) << std::right << step;
            
            std::string tree_name( "sa_tree_inc_" );
            tree_name += ss.str();
            
            std::string ali_name( "sa_ali_inc_" );
            ali_name += ss.str();
            
            builder.write_ali_and_tree( tree_name.c_str(), ali_name.c_str() );
            ++step;
            
            if( step == 2 ) {
                builder.prune_cloned_nodes();
                
            }
        }
        
        if( !valid ) {
            break;
        }
    }
#endif
    builder.prune_cloned_nodes();
    builder.write_ali_and_tree( "sa_tree_end", "sa_ali_end" );
    
    
//    builder.write_ali_and_tree_for_raxml();


}
