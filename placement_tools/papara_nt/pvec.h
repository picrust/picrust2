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


#ifndef __pvec_h
#define __pvec_h
//#define BOOST_UBLAS_NDEBUG
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <EigenvalueDecomposition.hpp>
#include <cassert>
#include <algorithm>
#include <iostream>

//#define USE_CBLAS
#ifdef USE_CBLAS
extern "C" {
#include <cblas.h>
}
#endif

#include "ivymike/tree_traversal_utils.h"
#include "parsimony.h"
#include "ivymike/stupid_ptr.h"
#include "ivymike/algorithm.h"
#include "sequence_model.h"
namespace {
//using ivy_mike::tip_case;
using ivy_mike::TIP_TIP;
using ivy_mike::TIP_INNER;
using ivy_mike::INNER_INNER;
}

class pvec_cgap {
    //     aligned_buffer<parsimony_state> v;
    std::vector<parsimony_state> v;
    std::vector<int> auxv;

    template<typename seq_model>
    class seq2aux {
    public:
        int operator()( parsimony_state c ) {

            if( seq_model::pstate_is_gap(c) ) {
                return AUX_CGAP;
            } else {
                return 0;
            }
        }
    };

public:
    
    
    void init( const std::vector<uint8_t> &seq ) {
//         assert( v.size() == 0 );
        v.resize(seq.size());
        auxv.resize( seq.size() );
        std::transform( seq.begin(), seq.end(), v.begin(), dna_parsimony_mapping::d2p );
        std::transform( seq.begin(), seq.end(), auxv.begin(), dna_parsimony_mapping::d2aux );
    }




    template<typename seq_model>
    void init2( const std::vector<uint8_t> &seq, const seq_model &sm ) {
        v.resize(seq.size());
        auxv.resize(seq.size());
        std::transform( seq.begin(), seq.end(), v.begin(), seq_model::s2p );
        std::transform( v.begin(), v.end(), auxv.begin(), seq2aux<seq_model>() );

        std::cerr << ">>>>>>>>>>>>>>>> WARNING: untested strange code!!!\n";
    }

    inline const std::vector<parsimony_state> &get_v() {
        return v;
    }
    inline const std::vector<int> &get_auxv() {
        return auxv;
    }

    static void newview( pvec_cgap &p, pvec_cgap &c1, pvec_cgap &c2, double /*z1*/, double /*z2*/, ivy_mike::tip_case tc ) {

    	if( c1.v.size() != c2.v.size() ) {
    		std::cout << "size2: " << c1.size() << " " << c2.size() << "\n";
    	}
        assert( c1.v.size() == c2.v.size() );

        if( c1.v.size() != c2.v.size() ) {
            std::cerr << "not equal: " << c1.size() << " " << c2.size() << "\n";
            throw std::runtime_error( "newview: vectors have different lengths (illegal incremetal newview on modified data?)" );
        }
        
//         p.v.resize(0);
        p.v.resize(c1.v.size());
        p.auxv.resize(c1.auxv.size());

        

        for( size_t i = 0; i < c1.v.size(); i++ ) {
            parsimony_state ps = c1.v[i] & c2.v[i];

            if( ps == 0 ) {
                ps = c1.v[i] | c2.v[i];
            }

            //p.v.push_back( ps );
            p.v[i] = ps;



            const int a1 = c1.auxv[i];
            const int a2 = c2.auxv[i];

            const bool cgap1 = (a1 & AUX_CGAP) != 0;
            const bool cgap2 = (a2 & AUX_CGAP) != 0;

//             const bool open1 = (a1 & AUX_OPEN) != 0;
//             const bool open2 = (a2 & AUX_OPEN) != 0;

            p.auxv[i] = 0;

            if( tc == TIP_TIP ) {
                if( cgap1 && cgap2 ) {
                    p.auxv[i] = AUX_CGAP;
                } else if( cgap1 != cgap2 ) {
                    p.auxv[i] = AUX_CGAP | AUX_OPEN;
                }
            } else if( tc == TIP_INNER ) {
                if( cgap1 && cgap2 ) {
                    p.auxv[i] = AUX_CGAP;
                } else if( cgap1 != cgap2 ) {
                    p.auxv[i] = AUX_CGAP | AUX_OPEN;
                }
            } else {
                if( a1 == AUX_CGAP && a2 == AUX_CGAP ) {
                    p.auxv[i] = AUX_CGAP;
                } else if( a1 == AUX_CGAP || a2 == AUX_CGAP ) {
                    p.auxv[i] = AUX_CGAP | AUX_OPEN;
                }
            }

        }
    }

    inline size_t size() {
        return v.size();
    }

    template<typename vt>
    inline void to_int_vec( std::vector<vt> &outv ) {

        outv.resize( v.size() );

        std::copy( v.begin(), v.end(), outv.begin() );
    }
    template<typename oiter_, size_t STRIDE>
    inline void to_int_vec_strided( oiter_ out ) {
        //outv.assign( v.begin(), v.end() );
        for( std::vector< parsimony_state >::iterator it = v.begin(); it != v.end(); ++it, out += STRIDE ) {
            *out = *it;
            
        }
        
//         outv.resize( v.size() );
//         std::copy( v.begin(), v.end(), outv.begin() );

    }
    
    
    template<typename vt>
    inline void to_aux_vec( std::vector<vt> &outv ) {
        //         std::cout << "v: " << v.size() << "\n";
        
        outv.resize( v.size() );
        std::copy( auxv.begin(), auxv.end(), outv.begin() );
        

//         std::for_each( auxv.begin(), auxv.end(), ostream_test(std::cout) );


    }
    
    template<typename oiter_, size_t STRIDE>
    inline void to_aux_vec_strided( oiter_ out ) {
        //         std::cout << "v: " << v.size() << "\n";
        
        for( std::vector< int >::iterator it = auxv.begin(); it != auxv.end(); ++it, out += STRIDE ) {
            if( *it == AUX_CGAP ) {
                *out = 0xFFFF;
            } else {
                *out = 0x0;
            }
        }
       // std::transform( auxv.begin(), auxv.end(), out );
        
        
        

//         std::for_each( auxv.begin(), auxv.end(), ostream_test(std::cout) );


    }
    
    bool operator==( const pvec_cgap &other ) const {
        return v == other.v && auxv == other.auxv;
    }
    
    bool operator!=( const pvec_cgap &other ) const {
        return !(*this == other);
    }
    
};

class probgap_model {

    boost::numeric::ublas::matrix<double> m_evecs;
    boost::numeric::ublas::vector<double> m_evals;
    boost::numeric::ublas::matrix<double> m_evecs_inv;

//    ublas::diagonal_matrix<double> m_evals_diag; // used temporarily.
//    ublas::matrix<double> m_prob_matrix;
//    ublas::matrix<double> m_temp_matrix;

    double m_gap_freq;
    bool m_valid;

    double calc_gap_freq ( const std::vector< std::vector< uint8_t > > &seqs ) {
        size_t ngaps = 0;
        size_t nres = 0;

        for( std::vector< std::vector< uint8_t > >::const_iterator it = seqs.begin(); it != seqs.end(); ++it ) {
            nres += it->size();
            ngaps += std::count_if( it->begin(), it->end(), dna_parsimony_mapping::is_gap );
        }

        double rgap = double(ngaps) / nres;
        std::cout << "gap rate: " << ngaps << " " << nres << "\n";
        std::cout << "gap rate: " << rgap << "\n";
        return rgap;
    }

public:
    probgap_model() : m_valid(false) {}

    probgap_model( const std::vector< std::vector<uint8_t> > &seqs ) : m_valid(false) {
    	reset( seqs );
    }

    probgap_model( double gap_freq ) : m_valid( false ) {
        reset( gap_freq );
    }
    
    void reset( const std::vector< std::vector<uint8_t> > &seqs ) {
	   // initialize probgap model from input sequences
    	reset( calc_gap_freq( seqs ) );
    }
    void reset( double gap_freq ) {
    	namespace ublas = boost::numeric::ublas;

    	m_gap_freq = gap_freq;
		double f[2] = {m_gap_freq, 1-m_gap_freq};
    	        //double f[2] = {1-m_gap_freq, m_gap_freq};
		ublas::matrix<double> rate_matrix(2,2);
		rate_matrix(0,0) = -f[0];
		rate_matrix(0,1) = f[0];
		rate_matrix(1,0) = f[1];
		rate_matrix(1,1) = -f[1];

		std::cout << "rate matrix: " << rate_matrix << "\n";

		ublas::EigenvalueDecomposition ed(rate_matrix);

		m_evecs = ed.getV();
		m_evals = ed.getRealEigenvalues();

		// use builtin ublas lu-factorization rather than jama
		{
			ublas::matrix<double> A(m_evecs);
			ublas::permutation_matrix<size_t> pm( A.size1() );
			size_t res = ublas::lu_factorize(A,pm);
			if( res != 0 ) {
				throw std::runtime_error( " ublas::lu_factorize failed" );
			}
			m_evecs_inv = ublas::identity_matrix<double>( A.size1());

			ublas::lu_substitute(A,pm,m_evecs_inv);
		}
		m_valid = true;
    }

    boost::numeric::ublas::matrix<double> setup_pmatrix( double t ) {
    	namespace ublas = boost::numeric::ublas;

    	ublas::diagonal_matrix<double> evals_diag(2);
		ublas::matrix<double> prob_matrix(2,2);
		ublas::matrix<double> temp_matrix(2,2);


        for( int i = 0; i < 2; i++ ) {
            evals_diag(i,i) = exp( t * m_evals(i));
        }

        prob_matrix = ublas::prod( ublas::prod( m_evecs, evals_diag, temp_matrix ), m_evecs_inv );
        //m_prob_matrix = ublas::prod( ublas::prod( m_evecs, m_evals_diag ), m_evecs_inv );

//         std::cout << "pmatrix: " << m_prob_matrix << "\n";



//         throw std::runtime_error( "xxx" );

        return prob_matrix;
    }


    inline double gap_freq() { return m_gap_freq; }

};

class pvec_pgap {
    //     aligned_buffer<parsimony_state> v;
    std::vector<parsimony_state> v;
    boost::numeric::ublas::matrix<double> gap_prob;

    
    
public:
    // WARNING WARNING WARNING: this is the stupid_pointer, used to inject a global probgap_model into class pvec_pgap
    static ivy_mike::stupid_ptr<probgap_model> pgap_model;

    inline const std::vector<parsimony_state> &get_v() const {
        return v;
    }

    inline const boost::numeric::ublas::matrix<double> &get_gap_prob() const {
        return gap_prob;
    }


    void init( const std::vector<uint8_t> &seq ) {
        //assert( v.size() == 0 );
        v.resize(seq.size());

        std::transform( seq.begin(), seq.end(), v.begin(), dna_parsimony_mapping::d2p );
        //std::transform( seq.begin(), seq.end(), auxv.begin(), dna_parsimony_mapping::d2aux );

        gap_prob.resize(2, seq.size());

        for( size_t i = 0; i < seq.size(); i++ ) {
            // ( 0 )
            // ( 1 ) means gap


            if( v[i] == 0xf ) {
                gap_prob( 0, i ) = 0.0;
                gap_prob( 1, i ) = 1.0;
            } else {
                gap_prob( 0, i ) = 1.0;
                gap_prob( 1, i ) = 0.0;
            }
        }

    }

    template<typename seq_model>
    void init2( const std::vector<uint8_t> &seq, seq_model sm ) {
        //assert( v.size() == 0 );
        v.resize(seq.size());

        std::transform( seq.begin(), seq.end(), v.begin(), seq_model::s2p );
        //std::transform( seq.begin(), seq.end(), auxv.begin(), dna_parsimony_mapping::d2aux );

        gap_prob.resize(2, seq.size());

        for( size_t i = 0; i < seq.size(); i++ ) {
            // ( 0 )
            // ( 1 ) means gap


            if( seq_model::pstate_is_gap(v[i]) ) {
                gap_prob( 0, i ) = 0.0;
                gap_prob( 1, i ) = 1.0;
            } else {
                gap_prob( 0, i ) = 1.0;
                gap_prob( 1, i ) = 0.0;
            }
        }



    }


    void newview( const pvec_pgap &c1, const pvec_pgap &c2, double z1, double z2, ivy_mike::tip_case tc ) {
    	newview( *this, c1, c2, z1, z2, tc );

    	//std::cout << "newview: " << gap_prob.size1() << "\n";
    }

    static void newview( pvec_pgap &p, const pvec_pgap &c1, const pvec_pgap &c2, double z1, double z2, ivy_mike::tip_case tc ) {
        namespace ublas = boost::numeric::ublas;
            
        assert( c1.v.size() == c2.v.size() );

        if( c1.v.size() != c2.v.size() ) {
            std::cerr << "not equal: " << c1.size() << " " << c2.size() << "\n";
            throw std::runtime_error( "newview: vectors have different lengths (illegal incremetal newview on modified data?)" );
        }

//         p.v.resize(0);
        p.v.resize(c1.v.size());
        //p.gap_prob.resize(2, c1.v.size());

        assert( pgap_model.is_valid_ptr() );

        ublas::matrix<double> p1 = pgap_model->setup_pmatrix(z1);
        ublas::matrix<double> p2 = pgap_model->setup_pmatrix(z2);



#if 1
        ublas::matrix<double> t1 = ublas::prod(p1, c1.gap_prob);
        ublas::matrix<double> t2 = ublas::prod(p2, c2.gap_prob);
#else
        ublas::matrix<double> t1(c1.gap_prob.size1(), c1.gap_prob.size2(),0 );
        ublas::matrix<double> t2(c2.gap_prob.size1(), c2.gap_prob.size2(),0 );    
        
        {
            double *p1x = p1.data().begin();
            double *p2x = p2.data().begin();
            
            double *c1x = c1.gap_prob.data().begin();
            double *c2x = c2.gap_prob.data().begin();
            
            double *t1bx = t1.data().begin();
            double *t2bx = t2.data().begin();
            
            size_t m = p1.size2();
            size_t n = c1.gap_prob.size2();
            size_t k = p1.size1();
            
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1, p1x, m, c1x, n, 1, t1bx, n);
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1, p2x, m, c2x, n, 1, t2bx, n);
        }
#endif   
        p.gap_prob = ublas::element_prod( t1, t2 );


        const static double twotothe256 = 115792089237316195423570985008687907853269984665640564039457584007913129639936.0;
                                                     /*  2**256 (exactly)  */

        const static double minlikelihood  = 1.0/twotothe256;

        for( ublas::matrix<double>::iterator2 it = p.gap_prob.begin2(); it != p.gap_prob.end2(); ++it ) {
        	double &p1 = *it;
        	double &p2 = *(it.begin() + 1);

        	if( fabs(p1) < minlikelihood && fabs(p2) < minlikelihood ) {
//        		std::cout << "scale: " << p1 << " " << p2;

        		p1 *= twotothe256;
        		p2 *= twotothe256;


//        		std::cout << " -> " << p1 << " " << p2 << "\n";
        	}


        }

        //ublas::matrix< double >::const_iterator1 xxx = p.gap_prob.begin1();
       // xxx.

        //std::cout << "pvec: " << *xxx << "\n";

//         ublas::matrix< double >::iterator1 tit1 = p.gap_prob.begin1();
//
//
//         boost::io::ios_all_saver ioss(std::cout);
//         std::cout << std::setprecision(2);
//
//         std::cout << std::left;
//
// //         std::transform( tit1.begin(), tit1.end(), tit1.begin(), log );
//         std::copy( tit1.begin(), tit1.end(), std::ostream_iterator<double>(std::cout, " "));
//         std::cout << "\n";
//         ++tit1;
// //         std::transform( tit1.begin(), tit1.end(), tit1.begin(), log );
//         std::copy( tit1.begin(), tit1.end() + 4, std::ostream_iterator<double>(std::cout, " "));
//         std::cout << "\n";


        for( size_t i = 0; i < c1.v.size(); i++ ) {
            parsimony_state ps = c1.v[i] & c2.v[i];

            if( ps == 0 ) {
                ps = c1.v[i] | c2.v[i];
            }

            //p.v.push_back( ps );
            p.v[i] = ps;
        }
    }

    inline size_t size() const {
        return v.size();
    }
    
    
    
    template<typename vt>
    inline void to_int_vec( std::vector<vt> &outv ) {
        outv.assign( v.begin(), v.end() );
//         outv.resize( v.size() );
//         std::copy( v.begin(), v.end(), outv.begin() );

    }
    

    static inline double gap_posterior( double v1, double v2 ) {
    	assert( pgap_model.is_valid_ptr() );
    	//return v1 / (v1 + v2);

    	v1 *= 1 - pgap_model->gap_freq();
    //	throw std::runtime_error( "i think there is an error in this function. why v1 in the next line?");
    	v2 *= pgap_model->gap_freq();

    	float v = float(v1 / (v1 + v2));

    	if( v != v ) {
     		std::cerr << "meeeeep: " << v1 << " " << v2 << "\n";

     		throw std::runtime_error("bailing out.");
    	}

    	return v;
    }

    static inline double gap_ancestral_probability( double v1, double v2 ) {
		assert( pgap_model.is_valid_ptr() );
		//return v1 / (v1 + v2);

		const double freq_ngap = 1 - pgap_model->gap_freq();
		const double freq_gap = pgap_model->gap_freq();


		return v2 * freq_gap / (v1 * freq_ngap + v2 * freq_gap);
	}



    inline void to_gap_post_vec( std::vector<double> &outv ) {

    	const boost::numeric::ublas::matrix<double> &t = get_pgap();

    	outv.clear();
    	outv.resize(v.size());

    	ivy_mike::binary_twizzle(t.begin2(), t.end2(), (t.begin1() + 1).begin(), outv.begin(), gap_posterior );

    }
   
    template<typename oiter>
    inline void to_ancestral_gap_prob( oiter it ) {
    	const boost::numeric::ublas::matrix<double> &t = get_pgap();
    	std::transform(t.begin2(), t.end2(), (t.begin1() + 1).begin(), it, gap_ancestral_probability );
    }


    template<typename vt>
    inline void to_aux_vec( std::vector<vt> &outv ) {
    
//         std::cout << "v: " << v.size() << "\n";


//    	{
//    		std::vector<double> postv;
//    		to_gap_post_vec(postv);
//
//
//
//    		//std::copy( postv.begin(), postv.end(), std::ostream_iterator<double>(std::cout, " " ));
//
//    		std::transform( postv.begin(), postv.end(), std::ostream_iterator<int>(std::cout), scaler<double>(10,0,9) );
//    		std::cout << "\n";
//    	}

        outv.clear();
        outv.reserve(v.size());
//         outv.resize( v.size() );
//         std::fill( outv.begin(), outv.end(), 0 );

//         std::for_each( auxv.begin(), auxv.end(), ostream_test(std::cout) );

        namespace ublas = boost::numeric::ublas;


        ublas::matrix<double> t = get_pgap();

        // yeah! metaprogramming massacre!

        ublas::matrix< double >::iterator1 tit1 = t.begin1();
#if 0
        std::vector<double> odds;
        odds.reserve(t.size2());

//         NOTE: t.begin1() + 1).begin() is the boost::ublas way of saying "iterator over the second row"
        ivy_mike::binary_twizzle( t.begin2(), t.end2(), (t.begin1() + 1).begin(), std::back_inserter(odds), std::divides<double>() );
        
        std::transform( odds.begin(), odds.end(), odds.begin(), std::ptr_fun<double>(log) );

        //std::transform( odds.begin(), odds.end(), std::back_inserter(outv), std::bind1st( std::greater<double>(), -1.0 ) );
        
        const double odds_threshold = 0.0;
        std::transform( odds.begin(), odds.end(), std::back_inserter(outv), std::bind2nd( std::less<double>(), odds_threshold ) );
        //         std::fill( outv.begin(), outv.end(), 1 );
#else
        // set outv[i] == 1 if, prob(non-gap) is less than prob(gap)
        // this is equivalent to the above code with odds_threshold == 0.0
        ivy_mike::binary_twizzle( t.begin2(), t.end2(), (t.begin1() + 1).begin(), std::back_inserter(outv), std::less<double>() );
        
                
#endif


#if 0
//         std::transform( t.begin2(), t.end2(), odds.begin(), calc_odds<ublas::matrix< double > > );

        //boost::io::ios_all_saver ioss(std::cout);
        std::cout << std::setprecision(0);
//         std::cout << std::setw(5);
        std::cout << std::fixed;
        std::cout << std::left;

//         std::cout << t(0,0) << t(0,1) << t(0,2) <<t(0,3) <<t(0,4) << "\n";


        //std::copy( odds.begin(), odds.end(), std::ostream_iterator<double>(std::cout, " "));
        std::copy( bias.begin(), bias.end(), std::ostream_iterator<int>(std::cout));
        std::cout << "\n";

        std::copy( lomag.begin(), lomag.end(), std::ostream_iterator<unsigned char>(std::cout));
        std::cout << "\n";



        std::transform( tit1.begin(), tit1.end(), tit1.begin(), log );
        std::copy( tit1.begin(), tit1.end(), std::ostream_iterator<double>(std::cout, " "));
        std::cout << "\n";
        ++tit1;
        std::transform( tit1.begin(), tit1.end(), tit1.begin(), log );
        std::copy( tit1.begin(), tit1.end(), std::ostream_iterator<double>(std::cout, " "));
        std::cout << "\n";
#endif

    }
//     template<typename oiter_, size_t STRIDE>
//     class strided_inserter {
//         typedef std::iterator_traits<oiter_> traits;
//         typedef typename traits::value_type reference_type;
//         
//         oiter_ out;
//     public:
//         strided_inserter( oiter_ out_ ) : out(out_) {}
//         
//         strided_inserter<oiter_,STRIDE> &operator=( reference_type v ) {
//             *out = v;
//         }
//         strided_inserter<oiter_,STRIDE> &operator*() {
//             return *this;
//         }
//         
//         strided_inserter<oiter_,STRIDE> &operator++() {
//             out+=STRIDE;
//             return *this;
//         }
//         
//         strided_inserter<oiter_,STRIDE> operator++(int w) {
//             out+=STRIDE * w;
//             return *this;
//         }
//         
//         
//     };
//     
//     template<typename oiter_, size_t STRIDE>
//     inline void to_aux_vec( oiter_ out ) {
// 
//         ublas::matrix<double> t = get_pgap();
// 
//         // yeah! metaprogramming massacre!
// 
//         ublas::matrix< double >::iterator1 tit1 = t.begin1();
// 
//         // set outv[i] == 1 if, prob(non-gap) is less than prob(gap)
//         // this is equivalent to the above code with odds_threshold == 0.0
//         ivy_mike::binary_twizzle( t.begin2(), t.end2(), (t.begin1() + 1).begin(), strided_inserter<oiter_,STRIDE>(out), std::less<double>() );
//     }
    
    const boost::numeric::ublas::matrix<double> &get_pgap() {
        return gap_prob;
    }


};


#endif
