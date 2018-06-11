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


#define BOOST_UBLAS_NDEBUG 1


//#include <sys/time.h>
//#include <sys/resource.h>
//#include <sys/syscall.h>
#include <sys/types.h>
//#include <unistd.h>

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <vector>
#include <deque>
#include <map>
#include <functional>
#include <cstring>
#include <limits>
#include <cmath>

#include <boost/io/ios_state.hpp>
#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/bind.hpp>
#include <boost/array.hpp>
//#include <boost/thread.hpp>
#include "ivymike/thread.h"


#include "parsimony.h"
#include "pvec.h"
#include "align_utils.h"
#include "ivymike/fasta.h"
#include "raxml_interface.h"

#include "math_approx.h"
#include "ivymike/tree_traversal_utils.h"
#include "ivymike/multiple_alignment.h"

#include "ivymike/tree_parser.h"
#include "ivymike/time.h"
#include "ivymike/getopt.h"
#include "ivymike/demangle.h"
#include "ivymike/stupid_ptr.h"
#include "ivymike/algorithm.h"
#include "ivymike/smart_ptr.h"

#include "papara.h"
#include "sequence_model.h"



using ivy_mike::tree_parser_ms::lnode;
using ivy_mike::tree_parser_ms::ln_pool;
using ivy_mike::tree_parser_ms::node_data_factory;

using ivy_mike::back_insert_ifer;
using ivy_mike::rooted_bifurcation;
using ivy_mike::apply_lnode;
using ivy_mike::iterate_lnode;


namespace ublas = boost::numeric::ublas;





//class my_adata : public ivy_mike::tree_parser_ms::adata {
////     static int ct;
//    //std::vector<parsimony_state> m_pvec;
//    pvec_t m_pvec;
//
//public:
////     int m_ct;
//    my_adata_gen() {
//
////         std::cout << "my_adata\n";
//
//    }
//
//    virtual ~my_adata_gen() {
//
////         std::cout << "~my_adata\n";
//
//    }
//
//    virtual void visit() {
////         std::cout << "tr: " << m_ct << "\n";
//    }
//    void init_pvec(const std::vector< uint8_t >& seq) {
//
//
//        m_pvec.init( seq );
////         std::cout << "init_pvec: " << m_pvec.size() << "\n";
////                 m_pvec.reserve(seq.size());
////         for( std::vector< uint8_t >::const_iterator it = seq.begin(); it != seq.end(); ++it ) {
////             m_pvec.push_back(dna_to_parsimony_state(*it));
////
////         }
//    }
//    pvec_t &get_pvec() {
//        return m_pvec;
//    }
//
//
//};

namespace {

const double g_delta = 0.0001;
const double g_epsilon = 0.001;
const double g_gap_freq = 0.715;


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

#if 0


class log_odds_aligner_score_only {
    typedef ublas::matrix<double> dmat;
    typedef std::vector<double> dsvec;


    // lof_t: log-odds-float = float type good enough to hold/calculate log-odds scores.
    // 32bit float should be enough
    typedef float lof_t;
    typedef ublas::matrix<lof_t> lomat;

    typedef std::vector<lof_t> losvec;

public:
    log_odds_aligner_score_only( const dmat &state, const dsvec &gap, boost::array<double,4> state_freq )
    : ref_state_prob_(state), ref_gap_prob_(gap), ref_len_(state.size1()),
      state_freq_(state_freq),
      neg_inf_( -std::numeric_limits<lof_t>::infinity() ),
      m_(state.size1() + 1),
      d_(state.size1() + 1),
      i_(state.size1() + 1),
      max_matrix_height_(0),
      delta_log_(log(g_delta)),
      epsilon_log_(log(g_epsilon))
    {
        precalc_log_odds();
    }

    void setup( size_t qlen ) {
        assert( ref_gap_prob_.size() == ref_len_ );



        // init first rows
        std::fill( m_.begin(), m_.end(), 0.0 );
        std::fill( d_.begin(), d_.end(), 0.0 );
        std::fill( i_.begin(), i_.end(), neg_inf_ );

        // init first columns
        m_[0] = 0.0;
        d_[0] = neg_inf_;
        i_[0] = 0.0;



    }


    void precalc_log_odds() {
        ref_state_lo_.resize( ref_state_prob_.size2(), ref_state_prob_.size1() );

        for( size_t i = 0; i < 4; ++i ) {
            const ublas::matrix_column<dmat> pcol( ref_state_prob_, i );
            ublas::matrix_row<lomat> lorow( ref_state_lo_, i );
            std::transform( pcol.begin(), pcol.end(), lorow.begin(), log_odds(state_freq_[i]));
        }



        {
            log_odds lo_ngap( 1 - g_gap_freq );
            log_odds lo_gap( g_gap_freq );

            ref_ngap_lo_.resize(ref_gap_prob_.size());
            ref_gap_lo_.resize(ref_gap_prob_.size());
            for( size_t i = 0; i < ref_gap_prob_.size(); ++i ) {
                ref_ngap_lo_[i] = lo_ngap(1 - ref_gap_prob_[i]);
                ref_gap_lo_[i] = lo_gap( ref_gap_prob_[i] );
            }
        }
    }

    template<typename T>
    static inline T max3( const T &a, const T &b, const T &c ) {
        return std::max( a, std::max( b, c ));
    }


    double align( const std::vector<uint8_t> &qs ) {
        const size_t qlen = qs.size();

        setup( qlen );

        //dmat ref_state_trans = trans(ref_state_prob_);


        assert( m_.size() == ref_len_ + 1 );




        for( size_t i = 1; i < qlen + 1; ++i ) {
            const int b = qs[i-1];
            //			std::cout << "b: " << b << "\n";

            //const double b_freq = state_freq_.at(b);
            //const ublas::matrix_column<dmat> b_state( ref_state_prob_, b );
            const ublas::matrix_row<lomat> b_state_lo( ref_state_lo_, b );

            //			const ublas::matrix_column<dmat> ngap_prob( ref_gap_prob_, 0 );
            //			const ublas::matrix_column<dmat> gap_prob( ref_gap_prob_, 1 );

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
            losvec::iterator rg = ref_gap_lo_.begin();
            losvec::iterator rng = ref_ngap_lo_.begin();

            const losvec::iterator m_end = m_.end();

            for( ; m0 != m_end; m1 = m0++, d1 = d0++, i1 = i0++, ++rg, ++rng, ++bsl ) {
                //ublas::matrix_row<dmat> a_state(ref_state_prob_, j-1 );
                //ublas::matrix_row<dmat> a_gap(ref_gap_prob_, j-1 );

                //double match_log_odds = log( b_state[j-1] / b_freq );
                //lof_t match_log_odds = b_state_lo[j-1];
                const lof_t match_log_odds = *bsl;



                //lof_t gap_log_odds = ref_gap_lo_[j-1];
                //lof_t ngap_log_odds = ref_ngap_lo_[j-1];
                const lof_t gap_log_odds = *rg;
                const lof_t ngap_log_odds = *rng;

                lof_t m_max = max3(
                        diag_m + ngap_log_odds,
                        diag_d + gap_log_odds,
                        diag_i + gap_log_odds
                );

                diag_m = *m0;
                *m0 = m_max + match_log_odds;

#if 0
                std::cout << i << " " << j << " " << m_(i,j) << " : " << m_(i-1, j-1) + ngap_log_odds
                        << " " << d_(i-1, j-1) + gap_log_odds << " " << i_(i-1, j-1) + gap_log_odds << " " << match_log_odds << " " << gap_log_odds << " " << ngap_log_odds << " max: " << m_max << "\n";
#endif

                diag_i = *i0;

                // the two 'diags' have already been updated, so they're both actually containing the current 'aboves',
                // which is exactly what we need to calculate the new i
                lof_t i_max = std::max(
                        diag_m + delta_log_,
                        diag_i + epsilon_log_
                );

                *i0 = i_max;

                lof_t d_max = std::max(
                        *m1 + delta_log_,
                        *d1 + epsilon_log_
                );

                diag_d = *d0;
                *d0 = d_max;


                //lof_t old_m = m_[j];

            }
        }
        {

            losvec::iterator max_it;

            max_it = std::max_element( m_.begin() + qlen, m_.end() );
            max_col_ = std::distance(m_.begin(), max_it);
            max_row_ = qlen;
            max_score_ = *max_it;



            return max_score_;
        }
    }

    dmat ref_state_prob_;
    dsvec ref_gap_prob_;

    lomat ref_state_lo_;
    losvec ref_gap_lo_;
    losvec ref_ngap_lo_;

    const size_t ref_len_;
    const boost::array<double,4> state_freq_;

    const float neg_inf_;

    losvec m_;
    losvec d_;
    losvec i_;

    size_t max_matrix_height_;

    const lof_t delta_log_;// = log(0.1);
    const lof_t epsilon_log_;// = log(0.5);

    size_t max_col_;
    size_t max_row_;
    double max_score_;

};
#endif

class odds {
public:

    odds( double bg_prob ) : bg_prob_(bg_prob) {}

    inline double operator()( double p ) {
        return p / bg_prob_;
    }

private:
    const double bg_prob_;
};


class log_odds_viterbi {
    typedef ublas::matrix<double> dmat;
    typedef std::vector<double> dsvec;


    // lof_t: log-odds-float = float type good enough to hold/calculate log-odds scores.
    // 32bit float should be enough
    typedef float lof_t;
    typedef ublas::matrix<lof_t> lomat;

    typedef std::vector<lof_t> losvec;

public:
    log_odds_viterbi( const dmat &state, const dsvec &gap, boost::array<double,4> state_freq )
    : ref_state_prob_(state), ref_gap_prob_(gap), ref_len_(state.size2()),
      state_freq_(state_freq),
      neg_inf_( -std::numeric_limits<lof_t>::infinity() ),
      m_(ref_len_ + 1),
      d_(ref_len_ + 1),
      i_(ref_len_ + 1),
      max_matrix_height_(0),
      delta_(g_delta),
      epsilon_(g_epsilon)
    {
        precalc_log_odds();
    }

    void setup( size_t qlen ) {
        
//         std::cerr << "len: " << ref_gap_prob_.size() << " " << ref_len_ << "\n";
        
        assert( ref_gap_prob_.size() == ref_len_ );



        // init first rows
        std::fill( m_.begin(), m_.end(), lof_t(0.0) );
        std::fill( d_.begin(), d_.end(), lof_t(0.0) );
        std::fill( i_.begin(), i_.end(), lof_t(0.0) /*neg_inf_*/ );

        // init first columns
        m_[0] = 0.0;
        d_[0] = 0.0; /*neg_inf_*/
        i_[0] = 0.0;



    }


    void precalc_log_odds() {
        ref_state_lo_.resize( ref_state_prob_.size1(), ref_state_prob_.size2() );

        for( size_t i = 0; i < 4; ++i ) {
            const ublas::matrix_row<dmat> pcol( ref_state_prob_, i );
            ublas::matrix_row<lomat> lorow( ref_state_lo_, i );
            std::transform( pcol.begin(), pcol.end(), lorow.begin(), log_odds(state_freq_[i]));
        }



        {
            odds odds_ngap( 1 - g_gap_freq );
            odds odds_gap( g_gap_freq );

            ref_ngap_odds_.resize(ref_gap_prob_.size());
            ref_gap_odds_.resize(ref_gap_prob_.size());
            for( size_t i = 0; i < ref_gap_prob_.size(); ++i ) {
                ref_ngap_odds_[i] = odds_ngap(1 - ref_gap_prob_[i]);
                ref_gap_odds_[i] = odds_gap( ref_gap_prob_[i] );
            }
        }
    }

    template<typename T>
    static inline T max3( const T &a, const T &b, const T &c ) {
        return std::max( a, std::max( b, c ));
    }


    double align( const std::vector<uint8_t> &qs ) {
        const size_t qlen = qs.size();

        setup( qlen );

        //dmat ref_state_trans = trans(ref_state_prob_);


        assert( m_.size() == ref_len_ + 1 );


        const bool verbose = false;
        if( verbose ) {
            std::cout << "viterbi:\n";
        }
        
        for( size_t i = 1; i < qlen + 1; ++i ) {
            const int b = qs[i-1];
            //          std::cout << "b: " << b << "\n";

            //const double b_freq = state_freq_.at(b);
            //const ublas::matrix_column<dmat> b_state( ref_state_prob_, b );
            const ublas::matrix_row<lomat> b_state_lo( ref_state_lo_, b );

            //          const ublas::matrix_column<dmat> ngap_prob( ref_gap_prob_, 0 );
            //          const ublas::matrix_column<dmat> gap_prob( ref_gap_prob_, 1 );

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
            losvec::iterator rg = ref_gap_odds_.begin();
            losvec::iterator rng = ref_ngap_odds_.begin();

            const losvec::iterator m_end = m_.end();

            for( ; m0 != m_end; m1 = m0++, d1 = d0++, i1 = i0++, ++rg, ++rng, ++bsl ) {
                //ublas::matrix_row<dmat> a_state(ref_state_prob_, j-1 );
                //ublas::matrix_row<dmat> a_gap(ref_gap_prob_, j-1 );

                //double match_log_odds = log( b_state[j-1] / b_freq );
                //lof_t match_log_odds = b_state_lo[j-1];
                const lof_t match_log_odds = *bsl;



                //lof_t gap_log_odds = ref_gap_lo_[j-1];
                //lof_t ngap_log_odds = ref_ngap_lo_[j-1];
                const lof_t gap_odds = *rg;
                const lof_t ngap_odds = *rng;

#if 0
                lof_t m_log_sum = log(
                           exp(diag_m) * ngap_odds
                         + exp(diag_d) * gap_odds
                         + exp(diag_i) * gap_odds
                );
                              
//                 std::cout << m_log_sum << " " << m_log_sum2 << "\n";
                
#else
                lof_t m_log_sum = diag_m + math_approx::log(
                           ngap_odds
                         + math_approx::exp(diag_d-diag_m) * gap_odds
                         + math_approx::exp(diag_i-diag_m) * gap_odds
                );
#endif

                diag_m = *m0;
                *m0 = m_log_sum + match_log_odds;
                if( verbose ) {
                    std::cout << std::setw(10) << *m0 << std::setw(10) << diag_d << std::setw(10) << diag_i << ";";
                    //                     std::cout << std::setw(10) << match_log_odds;
                }

               // assert( std::isfinite(*m0) );
				assert( *m0 == *m0 ); // FIXME: is this a valid replacement for std::isfinite?
#if 0
                std::cout << i << " " << j << " " << m_(i,j) << " : " << m_(i-1, j-1) + ngap_log_odds
                        << " " << d_(i-1, j-1) + gap_log_odds << " " << i_(i-1, j-1) + gap_log_odds << " " << match_log_odds << " " << gap_log_odds << " " << ngap_log_odds << " max: " << m_max << "\n";
#endif

                diag_i = *i0;

                // the two 'diags' have already been updated, so they're both actually containing the current 'aboves',
                // which is exactly what we need to calculate the new i
#if 0
                lof_t i_log_sum = log(
                          exp(diag_m) * delta_
                        + exp(diag_i) * epsilon_
                );
#else
                lof_t i_log_sum = diag_m + math_approx::log(
                        delta_
                      + math_approx::exp(diag_i-diag_m) * epsilon_
                );
#endif
                *i0 = i_log_sum;

#if 0
                lof_t d_log_sum = log(
                          exp(*m1) * delta_
                        + exp(*d1) * epsilon_
                );
#else
                lof_t d_log_sum = *m1 + math_approx::log(
                        delta_
                      + math_approx::exp(*d1 - *m1) * epsilon_
                );
#endif
                diag_d = *d0;
                *d0 = d_log_sum;


                //lof_t old_m = m_[j];

            }
            
            if( verbose ) {
                std::cout << "\n";
            }
        }

        
        float max_score = 0;
        {
            size_t lr_start = std::min( qlen, m_.size() - 1 ); 
            //auto max_it = std::max_element( m_.begin() + lr_start, m_.end() );
            
            for( size_t i = lr_start; i != m_.size(); ++i ) {
                max_score = max_score + log( 1 + exp( m_[i] - max_score ));
            }
            
            if( max_score == 0 ) {
                std::cout << "meeeep\n";
                
                std::copy( m_.begin() + lr_start, m_.end(), std::ostream_iterator<float>( std::cout, " " ));
                std::cout << "\n";
                getchar();
                
            }
            
//             assert( max_it != m_.end() );
//             max_score = *max_it;
            
        }
        return max_score;
//         return m_.back();
    }

    dmat ref_state_prob_;
    dsvec ref_gap_prob_;

    lomat ref_state_lo_;
    losvec ref_gap_odds_;
    losvec ref_ngap_odds_;

    const size_t ref_len_;
    const boost::array<double,4> state_freq_;

    const float neg_inf_;

    losvec m_;
    losvec d_;
    losvec i_;

    size_t max_matrix_height_;

    const lof_t delta_;// = log(0.1);
    const lof_t epsilon_;// = log(0.5);

    size_t max_col_;
    size_t max_row_;
    double max_score_;





};


template<bool NP>
class bin_log_odds {
public:
    inline double operator()( double v1, double v2 ) {
        double lo;

        if( !NP ) {
            lo = log( v1 / v2 );
        } else {
            lo = log( (1-v1) / (1-v2));
        }

        return std::max(-100.0, lo);
    }
};

class log_odds_aligner {
    typedef ublas::matrix<double> dmat;
    typedef std::vector<double> dsvec;


    // lof_t: log-odds-float = float type good enough to hold/calculate log-odds scores.
    // 32bit float should be enough
    typedef float lof_t;

    typedef ublas::matrix<lof_t> lomat;
    typedef std::vector<lof_t> losvec;

public:
    log_odds_aligner( const dmat &state, const dsvec &gap, boost::array<double,4> state_freq, const dsvec &pc_gap_freq )
    : ref_state_prob_(state), ref_gap_prob_(gap), ref_len_(state.size2()),
      state_freq_(state_freq),
      pc_gap_freq_( pc_gap_freq ),
      neg_inf_( -std::numeric_limits<lof_t>::infinity() ),
      max_matrix_height_(0),
      delta_log_(log(g_delta)),
      epsilon_log_(log(g_epsilon))
    {
        precalc_log_odds();
    }

    void setup( size_t qlen ) {
        assert( ref_gap_prob_.size() == ref_len_ );

        if( qlen + 1 > max_matrix_height_ ) {

            m_.resize( qlen + 1, ref_len_ + 1, false );
            d_.resize( qlen + 1, ref_len_ + 1, false );
            i_.resize( qlen + 1, ref_len_ + 1, false );
            max_matrix_height_ = qlen + 1;
        }
        // matrix organization: ref->cols, q->rows

        // init first rows
        std::fill( m_.begin1().begin(), m_.begin1().end(), lof_t(0.0) );
        std::fill( d_.begin1().begin(), d_.begin1().end(), lof_t(0.0) );
        std::fill( i_.begin1().begin(), i_.begin1().end(), neg_inf_ );

        // init first columns
        std::fill( m_.begin2().begin(), m_.begin2().end(), lof_t(0.0) );
        std::fill( d_.begin2().begin(), d_.begin2().end(), neg_inf_ );
        std::fill( i_.begin2().begin(), i_.begin2().end(), lof_t(0.0) );


    }

    void precalc_log_odds() {
        ref_state_lo_.resize( ref_state_prob_.size1(), ref_state_prob_.size2() );

        const bool print = false;
        
        for( size_t i = 0; i < 4; ++i ) {
            const ublas::matrix_row<dmat> pcol( ref_state_prob_, i );
            ublas::matrix_row<lomat> lorow( ref_state_lo_, i );
            std::transform( pcol.begin(), pcol.end(), lorow.begin(), log_odds(state_freq_[i]));
            
            if( print ) {
                std::cout << "lo row: " << i << "\n";
                if( i == 0 ) {
                    size_t j = 0;
                    std::for_each( lorow.begin(), lorow.end(), [&](double v) { std::cout << std::setw(10) << j++; });
                    std::cout << "\n";
                    
                }
                
                std::for_each( lorow.begin(), lorow.end(), [&](double v) { std::cout << std::setw(10) << v; });
                std::cout << "\n";
                
//                 getchar();

            }
        }
        
        if( print ) {
            throw std::runtime_error( "exit" );
        }
    //const double gap_freq = 0.83;


#if 0
        {
            log_odds lo_ngap( 1 - g_gap_freq );
            log_odds lo_gap( g_gap_freq );

            ref_ngap_lo_.resize(ref_gap_prob_.size());
            ref_gap_lo_.resize(ref_gap_prob_.size());
            for( size_t i = 0; i < ref_gap_prob_.size(); ++i ) {
                ref_ngap_lo_[i] = lo_ngap(1 - ref_gap_prob_[i]);
                ref_gap_lo_[i] = lo_gap( ref_gap_prob_[i] );
            }
        }
#endif
        ref_ngap_lo_.resize(ref_gap_prob_.size());
        ref_gap_lo_.resize(ref_gap_prob_.size());

        std::transform( ref_gap_prob_.begin(), ref_gap_prob_.end(), pc_gap_freq_.begin(), ref_ngap_lo_.begin(), bin_log_odds<true>());
        std::transform( ref_gap_prob_.begin(), ref_gap_prob_.end(), pc_gap_freq_.begin(), ref_gap_lo_.begin(), bin_log_odds<false>());

#if 0
        //		std::copy( ref_gap_lo_.begin(), ref_gap_lo_.end(), std::ostream_iterator<double>(std::cout, "\t"));
        //		std::cout << "\n";
        //		std::copy( ref_ngap_lo_.begin(), ref_ngap_lo_.end(), std::ostream_iterator<double>(std::cout, "\t"));
        //		std::cout << "\n";
        //		std::copy( ref_gap_prob_.begin(), ref_gap_prob_.end(), std::ostream_iterator<double>(std::cout, "\t"));
        //		std::cout << "\n";
        //		std::copy( pc_gap_freq_.begin(), pc_gap_freq_.end(), std::ostream_iterator<double>(std::cout, "\t"));
        //		std::cout << "\n";

        for( size_t i = 0; i < ref_gap_lo_.size(); ++i ) {
            std::cout << ref_gap_lo_[i] << "\t" << ref_ngap_lo_[i] << "\t" << ref_gap_prob_[i] << "\t" << pc_gap_freq_[i] << "\n";
        }
        std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
        //		throw std::runtime_error("");
#endif
    }

    template<typename T>
    static T max3( const T &a, const T &b, const T &c ) {
        return std::max( a, std::max( b, c ));
    }

    double align( const std::vector<uint8_t> &qs ) {
        const size_t qlen = qs.size();

        setup( qlen );

        //dmat ref_state_trans = trans(ref_state_prob_);

        const bool verbose = false;

        if( verbose ) {
            std::cout << "align\n";
        }
        
        //auto m_row = m_.begin1();
        
        for( size_t i = 1; i < qlen + 1; ++i ) {
            const int b = qs[i-1];
            //			std::cout << "b: " << b << "\n";

            //const double b_freq = state_freq_.at(b);
            //const ublas::matrix_column<dmat> b_state( ref_state_prob_, b );
            const ublas::matrix_row<lomat> b_state_lo( ref_state_lo_, b );

            //			const ublas::matrix_column<dmat> ngap_prob( ref_gap_prob_, 0 );
            //			const ublas::matrix_column<dmat> gap_prob( ref_gap_prob_, 1 );

            auto m_11 = m_.find2(0,i-1,0);
            auto d_11 = d_.find2(0,i-1,0);
            auto i_11 = i_.find2(0,i-1,0);
            
            auto m_00 = m_.find2(0,i,1);
            auto i_00 = i_.find2(0,i,1);
            
            auto m_01 = m_.find2(0,i,0);
            auto d_01 = d_.find2(0,i,0);
            
//             auto d_11 = d_.find2(i-1,0);
//             auto i_11 = i_.find2(i-1,0);
            
            for( size_t j = 1; j < ref_len_ + 1; ++j ) {
                //ublas::matrix_row<dmat> a_state(ref_state_prob_, j-1 );
                //ublas::matrix_row<dmat> a_gap(ref_gap_prob_, j-1 );

                //double match_log_odds = log( b_state[j-1] / b_freq );
                lof_t match_log_odds = b_state_lo[j-1];

                if( match_log_odds < -100 ) {
                    //	std::cout << "odd: " << match_log_odds << b_state[j-1] << " " << b_freq << " " << b << "\n";
                    match_log_odds = -100;
                }

                
#if 1
                lof_t gap_log_odds = ref_gap_lo_[j-1];
                lof_t ngap_log_odds = ref_ngap_lo_[j-1];
//                 std::cout << gap_log_odds << " " << ngap_log_odds << "\n";
//                 std::cout << "match: " << match_log_odds << "\n";
#else
                lof_t gap_log_odds =  0;
                lof_t ngap_log_odds = 0;
#endif
                
                lof_t m_max = max3(
                        *m_11 + ngap_log_odds,
                        *d_11 + gap_log_odds,
                        *i_11 + gap_log_odds
                );

                *m_00 = m_max + match_log_odds;

#if 0
                std::cout << i << " " << j << " " << m_(i,j) << " : " << m_(i-1, j-1) + ngap_log_odds
                        << " " << d_(i-1, j-1) + gap_log_odds << " " << i_(i-1, j-1) + gap_log_odds << " " << match_log_odds << " " << gap_log_odds << " " << ngap_log_odds << " max: " << m_max << "\n";
#endif
                        
                auto m_10 = ++m_11;
                auto i_10 = ++i_11;
                if( verbose ) {
                    std::cout << std::setw(10) << m_(i,j);
//                     std::cout << std::setw(10) << match_log_odds;
                }
                lof_t i_max = std::max(
                        *m_10 + delta_log_,
                        *i_10 + epsilon_log_
                );
                *i_00 = i_max;

                lof_t d_max = std::max(
                        *m_01 + delta_log_,
                        *d_01 + epsilon_log_
                );

                auto d_00 = ++d_01;
                *d_00 = d_max;
                
                ++d_11;
                ++m_00;
                ++i_00;
                ++m_01;
                

            }
            if( verbose ) {
                std::cout << "\n";
            }
        }
        {

            ublas::matrix_row<lomat> m_last(m_, qlen);
            ublas::matrix_row<lomat>::iterator max_it;
            
            
            // start point for searching the maximum in the last row
            // this is either qlen (=intersection of last row and diagonal) or the last
            // column (if qlen > reflen).
            size_t lr_start = std::min( qlen, m_last.size() - 1 ); 
            
            max_it = std::max_element( m_last.begin() + lr_start, m_last.end() );
            max_col_ = std::distance(m_last.begin(), max_it);
            max_row_ = qlen;
            max_score_ = *max_it;



            return max_score_;
        }
    }


    double align2( const std::vector<uint8_t> &qs ) {
        const size_t qlen = qs.size();

        setup( qlen );

        //dmat ref_state_trans = trans(ref_state_prob_);

        const bool verbose = false;

        if( verbose ) {
            std::cout << "align\n";
        }
        
        for( size_t i = 1; i < qlen + 1; ++i ) {
            const int b = qs[i-1];
            //          std::cout << "b: " << b << "\n";

            //const double b_freq = state_freq_.at(b);
            //const ublas::matrix_column<dmat> b_state( ref_state_prob_, b );
            const ublas::matrix_row<lomat> b_state_lo( ref_state_lo_, b );

            //          const ublas::matrix_column<dmat> ngap_prob( ref_gap_prob_, 0 );
            //          const ublas::matrix_column<dmat> gap_prob( ref_gap_prob_, 1 );

            for( size_t j = 1; j < ref_len_ + 1; ++j ) {
                //ublas::matrix_row<dmat> a_state(ref_state_prob_, j-1 );
                //ublas::matrix_row<dmat> a_gap(ref_gap_prob_, j-1 );

                //double match_log_odds = log( b_state[j-1] / b_freq );
                lof_t match_log_odds = b_state_lo[j-1];

                if( match_log_odds < -100 ) {
                    //  std::cout << "odd: " << match_log_odds << b_state[j-1] << " " << b_freq << " " << b << "\n";
                    match_log_odds = -100;
                }

                lof_t gap_log_odds = ref_gap_lo_[j-1];
                lof_t ngap_log_odds = ref_ngap_lo_[j-1];
                lof_t m_max = max3(
                        m_(i-1,j-1) + ngap_log_odds,
                        d_(i-1,j-1) + gap_log_odds,
                        i_(i-1,j-1) + gap_log_odds
                );

                m_(i,j) = m_max + match_log_odds;

#if 0
                std::cout << i << " " << j << " " << m_(i,j) << " : " << m_(i-1, j-1) + ngap_log_odds
                        << " " << d_(i-1, j-1) + gap_log_odds << " " << i_(i-1, j-1) + gap_log_odds << " " << match_log_odds << " " << gap_log_odds << " " << ngap_log_odds << " max: " << m_max << "\n";
#endif
                        
                        
                if( verbose ) {
                    std::cout << std::setw(10) << m_(i,j);
//                     std::cout << std::setw(10) << match_log_odds;
                }
                lof_t i_max = std::max(
                        m_(i-1,j) + delta_log_,
                        i_(i-1,j) + epsilon_log_
                );
                i_(i,j) = i_max;

                lof_t d_max = std::max(
                        m_(i,j-1) + delta_log_,
                        d_(i,j-1) + epsilon_log_
                );

                d_(i,j) = d_max;

            }
            if( verbose ) {
                std::cout << "\n";
            }
        }
        {

            ublas::matrix_row<lomat> m_last(m_, qlen);
            ublas::matrix_row<lomat>::iterator max_it;
            
            
            // start point for searching the maximum in the last row
            // this is either qlen (=intersection of last row and diagonal) or the last
            // column (if qlen > reflen).
            size_t lr_start = std::min( qlen, m_last.size() - 1 ); 
            
            max_it = std::max_element( m_last.begin() + lr_start, m_last.end() );
            max_col_ = std::distance(m_last.begin(), max_it);
            max_row_ = qlen;
            max_score_ = *max_it;



            return max_score_;
        }
    }
    
    void traceback( std::vector<uint8_t> *tb ) {
        assert( tb != 0 );

        //const ublas::matrix_column<dmat> ngap_prob( ref_gap_prob_, 0 );
        //const ublas::matrix_column<dmat> gap_prob( ref_gap_prob_, 1 );

        //		std::cout << "tb: max_col: " << max_col_ << "\n";

        size_t i = max_row_;
        size_t j = m_.size2() - 1;

        while( j > max_col_ ) {
            tb->push_back(1);
            --j;
        }

        bool in_d = false;
        bool in_i = false;
        while( j > 0 && i > 0 ) {

            if( in_d ) {
                assert( !in_i );
                tb->push_back(1);
                --j;
                in_d = m_(i,j) + delta_log_ < d_(i,j) + epsilon_log_;
            } else if( in_i ) {
                assert( !in_d );

                tb->push_back(2);

                --i;
                in_i = m_(i,j) + delta_log_ < i_(i,j) + epsilon_log_;
            } else {

                // the j - 1 corrects for the +1 offset, it does not mean 'last column'!
                const double gap_freq = 0.83;
                double ngap_log_odds = log( (1-ref_gap_prob_[j-1]) / (1-gap_freq) );
                double gap_log_odds = log( ref_gap_prob_[j-1] / gap_freq );
                gap_log_odds = std::max(-100.0, gap_log_odds);
                ngap_log_odds = std::max(-100.0, ngap_log_odds);


                tb->push_back(0);

                --i;
                --j;

                double vm = m_(i,j) + ngap_log_odds;
                double vi = i_(i,j) + gap_log_odds;
                in_i = vm < vi;

                if( in_i ) {
                    in_d = vi < d_(i,j) + gap_log_odds;
                    in_i = !in_d; // choose between d and i
                } else {
                    in_d = vm < d_(i,j) + gap_log_odds;
                }
            }

        }

        while( j > 0 ) {
            tb->push_back(1);
            --j;
        }
        while( i > 0 ) {
            tb->push_back(2);
            --i;
        }

    }

private:
    dmat ref_state_prob_;
    dsvec ref_gap_prob_;

    lomat ref_state_lo_;
    losvec ref_gap_lo_;
    losvec ref_ngap_lo_;

    const size_t ref_len_;
    const boost::array<double,4> state_freq_;
    const dsvec pc_gap_freq_;
    const float neg_inf_;

    lomat m_;
    lomat d_;
    lomat i_;

    size_t max_matrix_height_;

    const lof_t delta_log_;// = log(0.1);
    const lof_t epsilon_log_;// = log(0.5);

    size_t max_col_;
    size_t max_row_;
    double max_score_;


};

class my_adata : public ivy_mike::tree_parser_ms::adata {

public:

    typedef boost::numeric::ublas::matrix<double> apvecs;


    my_adata() : anc_gap_probs_valid_(false), max_depth_(0)
    {
    }

    void init_gap_vec(const std::vector< uint8_t >& seq) {
        gap_vec_.init(seq);

        update_ancestral_gap_prob();
    }

    void init_anc_prob_vecs( const apvecs &apv ) {
        anc_state_probs_ = apv; // if this uses too much memory, change to move semantics.
    }

    const pvec_pgap &get_gap_vec() const {
        return gap_vec_;
    }

    pvec_pgap *get_gap_vec_ptr() {
        anc_gap_probs_valid_ = false; // we are exposing internal state, so this _might_ invalidate the cache
        return &gap_vec_;
    }

    void update_ancestral_gap_prob() { // aahrg, this is stupid

        anc_gap_prob_.resize(gap_vec_.size() );
        gap_vec_.to_ancestral_gap_prob(anc_gap_prob_.begin());


        //    	anc_gap_probs_ = trans( gap_vec_.get_pgap() );
        anc_gap_probs_valid_ = true;
    }

    void print_vecs( std::ostream &os ) {


        const size_t len = anc_state_probs_.size1();

        assert( len == anc_gap_prob_.size() );


        for( size_t i = 0; i < len; ++i ) {
            os << "( ";
            for( size_t j = 0; j < 4; ++j ) {
                os << anc_state_probs_(i,j) << " ";
            }
            os << ") " << anc_gap_prob_[i] << "\n";
        }
    }

    const apvecs &state_probs() const {
        return anc_state_probs_;
    }

    const std::vector<double> &gap_probs() const {
        assert( anc_gap_probs_valid_ && "extra anal check because of the stupid cached transposed gap state vector." );
        return anc_gap_prob_;
    }

    size_t max_depth() const {
        return max_depth_;
    }

    void max_depth( size_t max_depth ) {
        max_depth_ = max_depth;
    }

private:
    apvecs anc_state_probs_;
    //apvecs anc_gap_probs_; // this is actually a transposed version of gap_vec_.get_pgap()
    pvec_pgap gap_vec_;
    std::vector<double> anc_gap_prob_;

    bool anc_gap_probs_valid_;
    size_t max_depth_;
};


class my_fact : public ivy_mike::tree_parser_ms::node_data_factory {

    virtual my_adata *alloc_adata() {
        return new my_adata;
    }

};


class queries {
public:
    typedef sequence_model::model<sequence_model::tag_dna4> seq_model;
    
    void load_fasta( const char *name ) {
        std::ifstream is( name );
        assert( is.good() );
        assert( names_.size() == raw_seqs_.size() );
        ivy_mike::read_fasta( is, names_, raw_seqs_ );
    }


    
    void add_move( const std::string &name, std::vector<uint8_t> &&seq ) {
        names_.push_back(name);
        //ivy_mike::push_back_swap(raw_seqs_, seq );
        raw_seqs_.emplace_back( seq );
    }



    void preprocess() {
        seqs_.resize(raw_seqs_.size());

        for( size_t i = 0; i < raw_seqs_.size(); ++i ) {
            seqs_.at(i).clear();

            // recode (ascii character to canonical state space) non-gap characters in raw_seq into seq
            
            
            std::for_each( raw_seqs_[i].begin(), raw_seqs_[i].end(),
                           [&]( uint8_t rc ) {
                               uint8_t cc = seq_model::s2c(rc);
                               if( !seq_model::cstate_is_gap(cc) ) {
                                   seqs_[i].push_back(cc);
                               }
                               
                           });
                            
//             std::cout << "qs: " << i << " " << names_.at(i) << "\n";
           
//             std::transform( m_qs_seqs[i].begin(), m_qs_seqs[i].end(),
//                             back_insert_ifer(m_qs_cseqs[i], std::not1( std::ptr_fun(seq_model::cstate_is_gap) )),
//                             seq_model::s2c );
//             std::transform( raw_seqs_[i].begin(), raw_seqs_[i].end(),
//                     back_insert_ifer(seqs_[i], std::bind2nd(std::less_equal<uint8_t>(), 3 ) ),
//                     encode_dna );

        }
    }

    size_t size() const {
        return names_.size();
    }
    const std::vector<uint8_t> &get_raw( size_t i ) {
        return raw_seqs_.at(i);
    }
    const std::vector<uint8_t> &get_recoded( size_t i ) const {
        return seqs_.at(i);
    }

    const std::string &get_name( size_t i ) {
        return names_.at(i);
    }

    size_t max_name_length() {
        size_t m = 0;
        for( std::vector<std::string>::iterator it = names_.begin(); it != names_.end(); ++it ) {
            m = std::max( m, it->size() );
        }
        return m;
    }

    #if 0
    static void seq_to_position_map( const std::vector<uint8_t> &seq, std::vector<int> *map ) {
        for( size_t i = 0; i < seq.size(); ++i ) {
            if( is_dna(seq[i]) ) {
                map->push_back(int(i));
            }
        }
    }


    static bool is_dna( uint8_t c ) {
        switch( c ) {
        case 'a':
        case 'c':
        case 'g':
        case 't':
        case 'A':
        case 'C':
        case 'G':
        case 'T':
        case 'U':
        case 'u':
            return true;
        default:
        {}
        }
        return false;
    }

    static uint8_t encode_dna( uint8_t c ) {
        switch( c ) {
        case 'a':
        case 'A':
            return 0;

        case 'c':
        case 'C':
            return 1;

        case 'g':
        case 'G':
            return 2;

        case 't':
        case 'T':
        case 'u':
        case 'U':
            return 3;

        default:
        {}
        }
        return std::numeric_limits<uint8_t>::max();
    }
#endif
private:
    std::vector<std::string> names_;
    std::vector<std::vector<uint8_t> > raw_seqs_; // seqs in the source alphabet (e.g. ACGT)

    //
    // stuff initialized by preprocess
    //
    std::vector<std::vector<uint8_t> > seqs_; // seqs in recoded alphabet (= 0 1 2 3...)
};

class references {
public:
    typedef sequence_model::model<sequence_model::tag_dna4> seq_model;
    
    references( std::shared_ptr<ln_pool> pool, const std::string &tree_name, const std::string &ali_name )
    : pool_(pool), tree_name_(tree_name), ali_name_(ali_name)
    {}


    void preprocess( queries &qs ) {

        std::vector<my_adata::apvecs> all_pvecs;

        ivy_mike::timer t1;

        lnode *t = generate_marginal_ancestral_state_pvecs(*pool_, tree_name_, ali_name_, &all_pvecs );

        std::cout << "generate: " << t1.elapsed() << "\n";

        std::vector<lnode *> nodes;
        iterate_lnode( t, std::back_inserter(nodes) );

        {
            // load reference seqs and filter out the seqs that are not in the tree
            // (they are added to the query seqs.)

            ivy_mike::multiple_alignment ref_ma;
            ref_ma.load_phylip( ali_name_.c_str() );


            std::vector<std::string> tip_names;

            for( std::vector<lnode *>::iterator it = nodes.begin(); it != nodes.end(); ++it ) {
                if( (*it)->m_data->isTip ) {
                    tip_names.push_back( (*it)->m_data->tipName );
                }
            }

            std::sort( tip_names.begin(), tip_names.end() );
            for( size_t i = 0; i < ref_ma.names.size(); ++i ) {
                if( std::binary_search(tip_names.begin(), tip_names.end(), ref_ma.names[i] )) {
                    ivy_mike::push_back_swap( ref_names_, ref_ma.names[i]);
                    ivy_mike::push_back_swap( ref_seqs_, ref_ma.data[i]);
                } else {
                    qs.add_move( ref_ma.names[i], std::move(ref_ma.data[i]) );
                }
            }
        }


        //
        // inject the probgap model ...
        //

        probgap_model pm( ref_seqs_ );
        std::cout << "p: " << pm.setup_pmatrix(0.1) << "\n";
        ivy_mike::stupid_ptr_guard<probgap_model> spg( pvec_pgap::pgap_model, &pm );


        //
        // initialize node/tip data
        //
        lnode *root = 0;
        for( std::vector<lnode *>::iterator it = nodes.begin(); it != nodes.end(); ++it ) {
            const std::string &nl = (*it)->m_data->nodeLabel;

            my_adata *mad = (*it)->m_data->get_as<my_adata>();

            if( !nl.empty() ) {
                assert( std::count_if(nl.begin(), nl.end(), isdigit)  == ptrdiff_t(nl.size()) );

                size_t nodenum = atoi( nl.c_str() );


                mad->init_anc_prob_vecs( all_pvecs.at(nodenum) );
            } else {
                assert( root == 0 );
                root = *it;
            }

            if( mad->isTip ) {
                size_t idx = std::find(ref_names_.begin(), ref_names_.end(), mad->tipName ) - ref_names_.begin();
                //				std::cout << "tip name: " << mad->tipName << " " << idx << "\n";
                //size_t idx = mp->second;
                mad->init_gap_vec( ref_seqs_.at(idx));
            }
        }



        assert( root != 0 );
        std::cout << "pvecs: " << all_pvecs.size() << "\n";

        base_freqs_ = count_base_freqs( ref_seqs_.begin(), ref_seqs_.end());


        //
        // generate rooted traversal order
        //

        trav_order_.clear();
        rooted_preorder_traversal( root->back, std::front_inserter(trav_order_), false );
        rooted_preorder_traversal( root->next->back, std::front_inserter(trav_order_), false );
        rooted_preorder_traversal( root->next->next->back, std::front_inserter(trav_order_), false );


        if( false ) { // unit test...
            std::deque<rooted_bifurcation<lnode> > trav_order_old;
            rooted_traversal_order( root->back, root->next->back, root->next->next->back, trav_order_old, false );

            assert( trav_order_.size() == trav_order_old.size() );
            assert( std::equal( trav_order_.begin(), trav_order_.end(), trav_order_old.begin()) );
        }


        //
        // do a whole-tree newview to calculate the gap probabilities
        //
        for( std::deque< rooted_bifurcation< ivy_mike::tree_parser_ms::lnode > >::iterator it = trav_order_.begin(); it != trav_order_.end(); ++it ) {
            //         std::cout << *it << "\n";

            my_adata *p = dynamic_cast<my_adata *>( it->parent->m_data.get());
            const my_adata *c1 = dynamic_cast<my_adata *>( it->child1->m_data.get());
            const my_adata *c2 = dynamic_cast<my_adata *>( it->child2->m_data.get());
            //         rooted_bifurcation<ivy_mike::tree_parser_ms::lnode>::tip_case tc = it->tc;

            //			std::cout << "tip case: " << (*it) << "\n";
            p->get_gap_vec_ptr()->newview( c1->get_gap_vec(), c2->get_gap_vec(), it->child1->backLen, it->child2->backLen, it->tc);

            p->update_ancestral_gap_prob();
            p->max_depth( 1 + std::max( c1->max_depth(), c2->max_depth() ));
            //p->cache_transposed_gap_probs();
            //p->print_vecs(std::cout);

        }



        calc_per_column_gap_freq();

        std::ofstream node_depth("node_depth.txt");
        for( size_t i = 0; i < node_size(); ++i ) {
            node_depth << i << " " << (get_node(i)->m_data->get_as<my_adata>()->max_depth()) << "\n";
        }

    }

    size_t node_size() const {
        return trav_order_.size();
    }
    const lnode *get_node( size_t i ) const {
        return trav_order_.at(i).parent;
    }

    const boost::array<double,4> &base_freqs() const {
        return base_freqs_;
    }

    const std::vector<std::string> &ref_names() const {
        return ref_names_;
    }

    const std::vector<std::vector<uint8_t> > &ref_seqs() const {
        return ref_seqs_;
    }
    
    size_t ref_len() const {
        return ref_seqs_.front().size();
    }

    const std::vector<double> &per_column_gap_freq() const {
        return per_column_gap_freq_;
    }

private:
    template<typename iiter>
    static boost::array<double,4> count_base_freqs( iiter start, iiter end ) {
        boost::array<double,4> counts;
        counts.fill(0);

        while( start != end ) {
            const std::vector<uint8_t> &seq = *start++;

            for( std::vector<uint8_t>::const_iterator it = seq.begin(); it != seq.end(); ++it ) {
                switch( *it ) {
                case 'a':
                case 'A':
                    ++counts[0];
                    break;

                case 'c':
                case 'C':
                    ++counts[1];
                    break;

                case 'g':
                case 'G':
                    ++counts[2];
                    break;

                case 't':
                case 'T':
                case 'u':
                case 'U':
                    ++counts[3];
                    break;

                default:
                {}
                }
            }
        }

        // this STL stuff is becoming kind of a bad habit ...
        // the following code divides the elements of 'counts' by the sum of all elements
        double all_count = std::accumulate( counts.begin(), counts.end(), 0 );
        std::transform( counts.begin(), counts.end(), counts.begin(), std::bind2nd(std::divides<double>(), all_count ));

        std::copy( counts.begin(), counts.end(), std::ostream_iterator<double>(std::cout, " "));
        std::cout << "\n";
        return counts;
    }


    void calc_per_column_gap_freq() {
        assert( !ref_seqs_.empty() );

        per_column_gap_freq_.assign(ref_seqs_.front().size(), g_gap_freq );
        return;

        std::vector<int> gap_counts( ref_seqs_.front().size(), 1 );

        for( std::vector<std::vector<uint8_t> >::iterator it = ref_seqs_.begin(); it != ref_seqs_.end(); ++it ) {
            assert( it->size() == gap_counts.size() );

            for( size_t i = 0; i < gap_counts.size(); ++i ) {
                if( !seq_model::sstate_is_character((*it)[i]) ) {
                    gap_counts[i]++;
                }
            }
        }

        per_column_gap_freq_.resize(gap_counts.size() );
        const double nseq = ref_seqs_.size();

        std::transform( gap_counts.begin(), gap_counts.end(), per_column_gap_freq_.begin(), std::bind2nd(std::divides<double>(),nseq + 2 ));

        //		std::copy( per_column_gap_freq_.begin(), per_column_gap_freq_.end(), std::ostream_iterator<double>( std::cout, " " ));
        //		std::cout << "\n";

    }
//     static bool is_dna( uint8_t c ) {
//         switch( c ) {
//             case 'a':
//             case 'c':
//             case 'g':
//             case 't':
//             case 'A':
//             case 'C':
//             case 'G':
//             case 'T':
//             case 'U':
//             case 'u':
//                 return true;
//             default:
//             {}
//         }
//         return false;
//     }
    
    std::shared_ptr<ln_pool> pool_;

    const std::string tree_name_;
    const std::string ali_name_;


    //
    // stuff initialized by preprocess
    //
    std::vector <std::string > ref_names_;
    std::vector <std::vector<uint8_t> > ref_seqs_;

    boost::array<double,4> base_freqs_;

    std::deque<rooted_bifurcation<lnode> > trav_order_;

    std::vector<double> per_column_gap_freq_;
};

std::string filename( const std::string &run_name, const char *type ) {
    std::stringstream ss;

    ss << "propara_" << type << "." << run_name;

    return ss.str();
}

bool file_exists(const char *filename)
{
    std::ifstream is(filename);
    return is.good();
}





struct scoring_results {


    scoring_results( size_t num_qs )
    : best_score_(num_qs, -std::numeric_limits<double>::infinity()),
      best_ref_(num_qs, size_t(-1))
//       os_("scores.txt")
    {}


    bool delta_equal( double v1, double v2 ) {
        return fabs( v1 - v2 ) < 0.0001; // TODO: what's a good delta for log-odds scores (and does anyone care...)
    }

    // the offer and set_trace methods are separate, so the lock can be release during traceback creation

    // if offer return true, this means that the qs/ref pair is currently best-scoring
    // the thread can crate the traceback and set it later with 'set_trace'
    bool offer( size_t qs, size_t ref, double score_ali, double score_vit ) {
//         boost::lock_guard<boost::mutex> lock(mtx_);

        if( qs == 0 ) {
            os_ << qs << " " << ref << " " << score_ali << " " << score_vit << "\n";
        }
//         const double score = score_vit;
        const double score = score_ali;
        if( best_score_.at(qs) < score || (delta_equal(best_score_.at(qs), score) && ref < best_ref_.at(qs))) {
            best_score_[qs] = score;
            best_ref_.at(qs) = ref;
            return true;
        }

        return false;
    }

//     // offer a traceback for a qs/ref pair. It is declined if it is no longer the best scoring pair
//     void set_trace( size_t qs, size_t ref, std::vector<uint8_t> *trace ) {
//         boost::lock_guard<boost::mutex> lock(mtx_);
// 
//         if( best_ref_.at(qs) == ref ) {
//             best_tb_.at(qs).swap(*trace);
//         } else {
//             std::cerr << "declined stale traceback\n";
//         }
//     }

    void merge( const scoring_results &other ) {
        ivy_mike::lock_guard<ivy_mike::mutex> lock(mtx_);
        
        const size_t size = best_score_.size();
        assert( size == best_ref_.size() );
        
        assert( other.best_score_.size() == other.best_ref_.size() ); // make sure that other is sane itself.
        
        if( size != other.best_score_.size() ) {
            throw std::runtime_error( "inconsistent sizes in scoring_result::merge" );
            
        }
        
        for( size_t i = 0; i < best_score_.size(); ++i ) {
            double other_score = other.best_score_[i];
            size_t other_ref = other.best_ref_[i];
            
            if( best_score_[i] < other_score || (delta_equal(best_score_[i], other_score) && other_ref < best_ref_[i])) {
                best_score_[i] = other_score;
                best_ref_[i] = other_ref;
            }
        }
        
    }
    
    std::vector<double> best_score_;
    std::vector<size_t> best_ref_;
//     std::vector<std::vector<uint8_t> > best_tb_;

    std::ofstream os_;

    ivy_mike::mutex mtx_;

};





class ScoringWorker {
public:

    ScoringWorker( const queries &qs, const references &refs, scoring_results *res, size_t rank, size_t num_workers )
    : qs_(qs),
      refs_(refs),
      res_(res),
      rank_(rank),
      num_workers_(num_workers)
    {}


    void operator()() {
        std::vector<uint8_t> tb_tmp;
        std::vector<uint8_t> rtb_tmp;
        
        uint64_t cups = 0;
        ivy_mike::timer t1;



        scoring_results local_res( res_->best_score_.size() );
        
        for( size_t i = rank_; i < refs_.node_size(); i+=num_workers_ ) {
            const lnode *a = refs_.get_node(i);
            assert( a->towards_root );
            assert( a->m_data != 0 );

            const my_adata *ma = dynamic_cast<const my_adata *>(a->m_data.get());

//             log_odds_viterbi vit( ma->state_probs(), ma->gap_probs(), refs_.base_freqs() );
            log_odds_aligner ali_score( ma->state_probs(), ma->gap_probs(), refs_.base_freqs(), refs_.per_column_gap_freq() );
            //log_odds_aligner_score_only ali_score( ma->state_probs(), ma->gap_probs(), refs_.base_freqs() );

            //std::cout << "ref: " << i << "\n";



            for( size_t j = 0; j < qs_.size(); ++j )
                //for( size_t x = 0, j = 9; x < 10; ++x, --j )
            {
                const std::vector<uint8_t> &b = qs_.get_recoded(j);


                cups += ma->state_probs().size2() * b.size();

                //double score2 = ali.align(b);
                double score = ali_score.align(b);

//                 double score_vit = vit.align(b);
                double score_vit = 0;
                
//                 std::cout << "score: " << score << " " << score_vit << "\n";
               
                if ( false ) {
                    tb_tmp.clear();
                    ali_score.traceback(&tb_tmp);
                    rtb_tmp.clear();
                    align_utils::realize_trace( b, tb_tmp, &rtb_tmp );
                
//                 os << std::setw(max_name_len+1) << std::left << qs.get_name(j);
                    std::copy( rtb_tmp.begin(), rtb_tmp.end(), std::ostream_iterator<char>(std::cout));
                    std::cout << "\n";
                }
                //assert( score == score2 );
//                 std::cout << "score: " << j << " " << i << " " << score << "\n";
                
                
                
                if( local_res.offer(j, i, score, score_vit ) ) {
                    //tb_tmp.clear();
                    //	ali.traceback(&tb_tmp);

                    //res_->set_trace(j, i, &tb_tmp );
                }

                //			if( j == 0 ) {
                //				std::cout << "score: " << score << "\n";
                //			}
                
//                 getchar();
            }



            //std::cout << "score: " << score << "\n";

        }
        
        res_->merge(local_res);
        std::cout << "time: " << t1.elapsed() << " " << cups / (t1.elapsed()*1e9) << " gncup/s\n";
    }

private:
    const queries &qs_;
    const references &refs_;
    scoring_results *res_;
    size_t rank_;
    size_t num_workers_;
};


size_t write_ref_seqs( std::ostream &os, const references &refs, size_t num_qs, size_t max_qs_name_len ) {
    os << refs.ref_seqs().size() + num_qs << " " << refs.ref_seqs().front().size() << "\n";


    size_t max_name_len = max_qs_name_len;
    for( size_t i = 0; i < refs.ref_names().size(); ++i ) {
        max_name_len = std::max( max_name_len, refs.ref_names()[i].size() );
    }


    for( size_t i = 0; i < refs.ref_names().size(); ++i ) {
        os << std::setw(max_name_len+1) << std::left << refs.ref_names()[i];

        std::copy( refs.ref_seqs()[i].begin(), refs.ref_seqs()[i].end(), std::ostream_iterator<char>(os));
        os << "\n";
    }

    return max_name_len;
}

int main( int argc, char *argv[] ) {
    namespace igo = ivy_mike::getopt;
    ivy_mike::getopt::parser igp;
    std::string opt_tree_name;
    std::string opt_alignment_name;
    std::string opt_qs_name;
    bool opt_use_cgap;
    int opt_num_threads;
    std::string opt_run_name;
    bool opt_write_testbench;
    bool opt_use_gpu;
    bool opt_force_overwrite;
    igp.add_opt('t', igo::value<std::string>(opt_tree_name));
    igp.add_opt('s', igo::value<std::string>(opt_alignment_name));
    igp.add_opt('q', igo::value<std::string>(opt_qs_name));
    igp.add_opt('c', igo::value<bool>(opt_use_cgap, true).set_default(false));
    igp.add_opt('j', igo::value<int>(opt_num_threads).set_default(1));
    igp.add_opt('n', igo::value<std::string>(opt_run_name).set_default("default"));
    igp.add_opt('b', igo::value<bool>(opt_write_testbench, true).set_default(false));
    igp.add_opt('g', igo::value<bool>(opt_use_gpu, true).set_default(false));
    igp.add_opt('f', igo::value<bool>(opt_force_overwrite, true).set_default(false));
    igp.parse(argc, argv);
    if(igp.opt_count('t') != 1 || igp.opt_count('s') != 1){
        std::cerr << "missing options -t and/or -s (-q is optional)\n";
        return 0;
    }
    ivy_mike::timer t;
    const char *qs_name = 0;
    if(!opt_qs_name.empty()){
        qs_name = opt_qs_name.c_str();
    }
    std::string log_filename = filename(opt_run_name, "log");
    if(opt_run_name != "default" && !opt_force_overwrite && file_exists(log_filename.c_str())){
        std::cout << "log file already exists for run '" << opt_run_name << "'\n";
        return 0;
    }
    
    papara::add_log_tee log_cout(std::cout);
    std::ofstream logs(log_filename.c_str());
    if(!logs){
        std::cout << "could not open logfile for writing: " << log_filename << std::endl;
        return 0;
    }
    
    papara::add_log_tee log_file(logs);
    
    
    std::shared_ptr<ln_pool> pool(new ln_pool(ln_pool::fact_ptr_type(new my_fact)));
    queries qs;
    if(qs_name != 0){
        qs.load_fasta(qs_name);
    }
    references refs(pool, opt_tree_name, opt_alignment_name);
    refs.preprocess(qs);
    qs.preprocess();



    //    std::vector<double> best_score(qs.size(), -std::numeric_limits<double>::infinity());
    //    std::vector<size_t> best_ref(qs.size(), size_t(-1));
    //    std::vector<std::vector<uint8_t> > best_tb(qs.size());


    uint64_t cups = 0;

    ivy_mike::timer t1;

    scoring_results res( qs.size() );

    ScoringWorker w0(qs, refs, &res, 0, opt_num_threads );

    ivy_mike::thread_group tg;
    for( int i = 1; i < opt_num_threads; ++i ) {
        std::cout << "starting additional thread: " << i << "\n";
        tg.create_thread(ScoringWorker(qs, refs, &res, i, opt_num_threads ));
    }

    w0();

    tg.join_all();

    std::cout << "time: " << t1.elapsed() << " " << cups / (t1.elapsed()*1e9) << " gncup/s\n";


//     std::vector<uint8_t> tb;
    
    std::string out_name(filename( opt_run_name, "alignment" ));
//     std::ofstream os( out_name.c_str());

    std::string qual_name(filename( opt_run_name, "quality" ));
    std::ofstream os_qual( qual_name.c_str());



//     const size_t max_name_len = write_ref_seqs( os, refs, qs.size(), qs.max_name_length() );


    double qual = 0.0;
    int num_qual = 0;
    
    size_t max_name_len = 0;
    
    std::vector<std::vector<uint8_t> > qs_traces(qs.size());
    std::vector<double > qs_scores(qs.size());
    for( size_t i = 0; i < refs.node_size(); ++i )
    {




        const lnode *a = refs.get_node(i);
        assert( a->towards_root );
        assert( a->m_data != 0 );

        const my_adata *ma = dynamic_cast<const my_adata *>(a->m_data.get());

        assert( res.best_ref_.size() == qs.size() );
        log_odds_aligner ali( ma->state_probs(), ma->gap_probs(), refs.base_freqs(), refs.per_column_gap_freq() );

        
        max_name_len = std::max( max_name_len, refs.ref_names().at(i).size() );
        
        for( size_t j = 0; j < qs.size(); ++j ) {

            if( i == 0 ) {
                // on the first round through the refs also search for the maximum qs name length
                max_name_len = std::max( max_name_len, qs.get_name(j).size() );
            }
            
            if( res.best_ref_[j] != i ) {
                continue;
            }

            const std::vector<uint8_t> &b = qs.get_recoded(j);

            double score = ali.align(b);
            qs_scores.at(j) = score;
            //assert( score == res.best_score_.at(j) );

//             tb.clear();
//             ali.traceback(&tb);
            ali.traceback( &qs_traces.at(j) );
        }
    }




    papara::output_alignment_phylip oa( out_name.c_str() );

    typedef sequence_model::model<sequence_model::tag_dna4> seq_model;
    
    // collect ref gaps introduiced by qs
    const size_t pad = max_name_len + 1;
    const bool ref_gaps = true;

    papara::ref_gap_collector rgc( refs.ref_len() );
    for( std::vector<std::vector<uint8_t> >::iterator it = qs_traces.begin(); it != qs_traces.end(); ++it ) {
        rgc.add_trace(*it);
    }



    const size_t num_refs = refs.ref_seqs().size();
    
    if( ref_gaps ) {
        oa.set_size(num_refs + qs.size(), rgc.transformed_ref_len());
    } else {
        oa.set_size(num_refs + qs.size(), refs.ref_len() );
    }
    oa.set_max_name_length( pad );
    
    // write refs (and apply the ref gaps)

    std::vector<char> tmp;
    for( size_t i = 0; i < refs.ref_seqs().size(); i++ ) {
        tmp.clear();
        
        

        if( ref_gaps ) {
            rgc.transform( refs.ref_seqs().at(i).begin(), refs.ref_seqs().at(i).end(), std::back_inserter(tmp), '-' );
        } else {
            std::transform( refs.ref_seqs().at(i).begin(), refs.ref_seqs().at(i).end(), std::back_inserter(tmp), seq_model::normalize);
        }

        oa.push_back( refs.ref_names().at(i), tmp, papara::output_alignment::type_ref );
        //std::transform( m_ref_seqs[i].begin(), m_ref_seqs[i].end(), std::ostream_iterator<char>(os), seq_model::normalize );
        
    }



    std::vector<uint8_t>out_qs_cstate;
    std::vector<char>out_qs_char;
    for( size_t i = 0; i < qs.size(); i++ ) {
        tmp.clear();
        const std::vector<uint8_t> &qp = qs.get_recoded(i);

        out_qs_cstate.clear();

        // realize+write the tracebacks. this is done in cstate (=canonical state) space of the dna4 model (=the
        // sequences don't contain gaps/ambiguous states)
        
        
        if( ref_gaps ) {
            papara::gapstream_to_alignment(qs_traces.at(i), qp, &out_qs_cstate, seq_model::gap_cstate(), rgc);
        } else {
            papara::gapstream_to_alignment_no_ref_gaps(qs_traces.at(i), qp, &out_qs_cstate, seq_model::gap_cstate() );
        }

        //os << std::setw(pad) << std::left << qs.name_at(i);
//         std::transform( out_qs_ps.begin(), out_qs_ps.end(), std::back_inserter(tmp), seq_model::p2s );
        
        
        // FIXME: this looks like a design quirk: why is papara::output_alignment based on char sequences while
        // everywhere else uint8_t seqeunces?
        //out_qs_char.assign( out_qs.begin(), out_qs.end() );
        
        // transform from canonical state space into (ascii) sequence space
        
        out_qs_char.clear();
        std::transform( out_qs_cstate.begin(), out_qs_cstate.end(), std::back_inserter(out_qs_char), seq_model::c2s );
        
        oa.push_back( qs.get_name(i), out_qs_char, papara::output_alignment::type_qs );




#if 0
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
#endif

    }

    std::cout << t.elapsed() << std::endl;
    papara::lout << "SUCCESS " << t.elapsed() << std::endl;

    std::cout << "mean quality: " << qual / num_qual << "\n";

    return 0;
}

