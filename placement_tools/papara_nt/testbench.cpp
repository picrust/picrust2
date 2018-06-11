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
#include <stdexcept>
#include <vector>
#include <stdint.h>
#include <climits>
#include <algorithm>
#include <iterator>
#include <utility>
#include <map>
#include <cassert>

#include <sys/time.h>

#include <boost/thread.hpp>
#include <boost/shared_ptr.hpp>

#include "align_pvec_vec.h"

static double gettime(void )
{
    struct timeval ttime;
    gettimeofday(&ttime , 0);
    return ttime.tv_sec + ttime.tv_usec * 0.000001;
}


class timer {
    double m_st;
    
public:
    timer() : m_st(gettime()) {};
    
    double elapsed() const {
        return gettime() - m_st;
    }
    
    
};



typedef int score_t;


struct align_arrays {
    std::vector<score_t> s;
};


static score_t align_pvec_score_seq( std::vector<int> &a, std::vector<unsigned int> &a_aux, const std::vector<uint8_t> &b, score_t mismatch_score, score_t match_cgap, score_t gap_open, score_t gap_extend, align_arrays &arr ) {

    if( arr.s.size() < a.size()  ) {
        arr.s.resize( a.size() );
    }
    assert( a.size() > b.size() );
    const score_t aux_cgap = 0x1;
    
    const size_t band_width = a.size() - b.size();
    std::fill( arr.s.begin(), arr.s.end(), 0 );
    
    const score_t LARGE = 32000;
//     std::fill( arr.si.begin(), arr.si.end(), LARGE );
    arr.s.at(band_width+1) = LARGE;

//     score_t min_score = LARGE;
    
    for( size_t ib = 0; ib < b.size(); ib++ ) {
        int bc = b[ib];
        
        score_t last_sl = LARGE;
        score_t last_sc = LARGE;
        score_t last_sdiag = 0.0;
        
        size_t astart = ib;
                
        score_t * __restrict s_iter = &arr.s[0];

        
//        bool lastrow = ib == (b.size() - 1);
        
       
        last_sdiag = *s_iter;
        
        for( size_t i = 0; i < astart; ++i ) {
            std::cout << "   ";
        }
        for( size_t ia = astart; ia <= ib + band_width; ++ia, ++s_iter ) {  
            score_t ac = a[ia];
            const bool cgap = int(a_aux[ia]) == aux_cgap;
            
            // determine match or mis-match according to parsimony bits coming from the tree.
            bool is_match = ( ac & bc ) != 0;

            score_t sm = last_sdiag;
            
            if( !is_match ) {
                sm += mismatch_score;
            }
            
            last_sdiag = *(s_iter+1);

            score_t last_sc_OPEN; 
            score_t sl_score_stay;
            
            if( cgap ) {
                last_sc_OPEN = last_sc; 
                sl_score_stay = last_sl;
                sm += match_cgap;
            } else {
                last_sc_OPEN = last_sc + gap_open + gap_extend;  
                sl_score_stay = last_sl + gap_extend;
            }
            
            score_t sl = std::min( sl_score_stay, last_sc_OPEN );
                        
            last_sl = sl;
            

            score_t su = last_sdiag + mismatch_score;
            
            //std::cout << "x: " << sm << " " << su << " " << sl << "  ";
            
            score_t sc = std::min( sm, std::min( su, sl ) );
            std::cout << std::setw(3) << sc;
            last_sc = sc;
            *s_iter = sc;
            
        
//             if( lastrow ) {
//                 if( sc < min_score ) {
//                     
//                     min_score = sc;
//                 }
//             }
           
        }
        std::cout << "\n";
        
        
        
//         std::cout << "\n";
//         getchar();
    }
    
    return *std::min_element( arr.s.begin(), arr.s.begin() + band_width + 1 );
    //std::cout << "min2: " << min_score << " " << min2 << "\n";

    //return min_score;
    
}




class testbench {
    const static int score_gap_open = 1;
    const static int score_gap_extend = 1;
    const static int score_mismatch = 1;
    const static int score_match_cgap = 4;
    
    std::vector<std::vector <int> > m_ref_pvecs;
    std::vector<std::vector <unsigned int> > m_ref_aux;
    
    std::vector<std::vector <uint8_t> > m_qs_pvecs;

    std::map<std::pair<size_t,size_t>,int> m_seq_res;
    std::map<std::pair<size_t,size_t>,int> m_vec_res;
    
    
    void load_qs( std::istream &is ) {
        size_t num_qs;
        is >> num_qs;
        
        m_qs_pvecs.resize( num_qs );
        
        for( size_t i = 0; i < m_qs_pvecs.size(); ++i ) {
            size_t len;
            is >> len;
            while( isspace( is.get() ) && !is.eof() ) {}
            is.unget();
            
            if( is.eof() ) {
                throw std::runtime_error( "unexpected end of file" );
            }
            
            m_qs_pvecs[i].resize(len);
            is.read( (char*)&m_qs_pvecs[i][0], len );
//             std::cout << "qs: " << i << " " << m_qs_pvecs[i].size() << "\n";
        }
    }
    
    void load_ref( std::ifstream &is ) {
        size_t num_ref;
        is >> num_ref;
        m_ref_pvecs.resize( num_ref );
        m_ref_aux.resize( num_ref );
        
        for( size_t i = 0; i < m_ref_pvecs.size(); ++i ) {
            size_t len;
            is >> len;
            while( isspace( is.get() ) && !is.eof() ) {}
            is.unget();
            
            if( is.eof() ) {
                throw std::runtime_error( "unexpected end of file" );
            }
            
            m_ref_pvecs[i].resize(len);
            is.read( (char*)&m_ref_pvecs[i][0], len * sizeof(int) );
            
            m_ref_aux[i].resize(len);
            is.read( (char*)&m_ref_aux[i][0], len * sizeof(unsigned int) );
            
//             std::cout << "ref: " << i << " " << m_ref_pvecs[i].size() << "\n";
        }
        
    }
    
    template<typename T, size_t min_len>
    static void trim_ref( std::vector<T> &v ) {
        
    	const size_t max_len = min_len;
        while( v.size() < min_len ) {
            std::vector<T> v_new = v;
            
            
            std::copy( v.begin(), v.end(), std::back_inserter(v_new) );
            
            v.swap(v_new);
            
        }
        

        if( v.size() > max_len ) {
            v.resize( max_len );
        }
        
        
    }
    template<typename T, size_t max_len>
    static void trim_max( std::vector<T> &v ) {

        if( v.size() > max_len ) {
            v.resize( max_len );
        }
        
        
    }
public:  
    testbench( const char *qs_name, const char *ref_name, size_t num_ref, size_t num_qs ) {
        {
            std::ifstream is( qs_name );
            load_qs( is );
        }
        {
            std::ifstream is( ref_name );
            load_ref( is );
        }
        
        
        // reduce size of testset
//        size_t num_ref = 192;
//        size_t num_qs = 200;
      
#if 0
        const size_t ref_len = 128;
        const size_t qs_len = ref_len - 1;
        std::for_each( m_ref_pvecs.begin(), m_ref_pvecs.end(), &trim_ref<int,ref_len> );
        std::for_each( m_ref_aux.begin(), m_ref_aux.end(), &trim_ref<unsigned int,ref_len> );
        std::for_each( m_qs_pvecs.begin(), m_qs_pvecs.end(), &trim_max<uint8_t,qs_len> );
#endif

//        for( size_t i = 0; i < m_ref_pvecs.size(); ++i ) {
//        	const size_t ref_len = 4000;
//
//        	assert( m_ref_pvec)
//
//        }

        num_ref = std::min( num_ref, m_ref_pvecs.size() );
        num_qs = std::min( num_qs, m_qs_pvecs.size() );
        
        m_ref_pvecs.resize(num_ref);
        m_ref_aux.resize(num_ref);
        m_qs_pvecs.resize(num_qs);
        
    }
    
    
    void run_seq() {
        align_arrays arr;
        
        uint64_t ncup = 0;
        
        timer t1;
        
        for( size_t i = 0; i < m_ref_pvecs.size(); ++i ) {
            
            const size_t reflen = m_ref_pvecs[i].size();
            
            for( size_t j = 0; j < m_qs_pvecs.size(); ++j ) {
                

//                int score2 = 0;
                int score = align_pvec_score_seq( m_ref_pvecs[i], m_ref_aux[i], m_qs_pvecs[j], score_mismatch, score_match_cgap, score_gap_open, score_gap_extend, arr );
                //std::cout << "score " << i << " " << j << " = " << score << " " << score2 << "\n";
//                 if( score != score2 ) {
   //                 std::cout << "score " << i << " " << j << " = " << score << "\n";
//                 }
                ncup += reflen * m_qs_pvecs[i].size();
                
                m_seq_res[std::pair<size_t,size_t>(i,j)] = score;
            }
            
            
        }
        
        std::cerr << "aligned " << m_ref_pvecs.size() << " x " << m_qs_pvecs.size() << " seqs in " << t1.elapsed() << " s\n";
        std::cerr << "GNcups: " << ncup / (t1.elapsed() * 1e9) << "\n";
    }
    
    template<typename score_t, const size_t W>
    void run_vec() {
        align_arrays arr;
        
        uint64_t ncup = 0;
        
        timer t1;
        
        
        for( size_t i = 0; i < m_ref_pvecs.size(); i += W ) {
            const int *seqptrs[W];
            const unsigned int *auxptrs[W];
            size_t num_valid = 0;
            for( size_t j = 0; j < W; ++j ) {
                if( i + j < m_ref_pvecs.size() ) {
                    seqptrs[j] = m_ref_pvecs[i+j].data();
                    auxptrs[j] = m_ref_aux[i+j].data();
                    ++num_valid;
                } else {
                    seqptrs[j] = seqptrs[0];
                    auxptrs[j] = auxptrs[0];
                }
            }
            
            assert( num_valid > 0 );
            
            const size_t reflen = m_ref_pvecs[i].size();
            
            const align_pvec_score<score_t,W> aligner(seqptrs, auxptrs, reflen, score_mismatch, score_match_cgap, score_gap_open, score_gap_extend );
            
            
//             std::cerr << "ref: " << i << "\n";
            for( size_t j = 0; j < m_qs_pvecs.size(); ++j ) {
                
                aligner.align( m_qs_pvecs[j] );
                const score_t *scores = aligner.get_scores();
                
                for( size_t k = 0; k < num_valid; ++k ) {
                
           //         std::cout << "score " << i + k << " " << j << " = " << scores[k] << "\n";
                    m_vec_res[std::pair<size_t,size_t>(i+k,j)] = scores[k];
                }
                ncup += num_valid * reflen * m_qs_pvecs[i].size();
            }
            
            
        }
        
        std::cerr << "aligned " << m_ref_pvecs.size() << " x " << m_qs_pvecs.size() << " seqs in " << t1.elapsed() << " s\n";
        std::cerr << "GNcups: " << ncup / (t1.elapsed() * 1e9) << "\n";
    }
    
    

    template<typename score_t, const size_t W>
    class worker {
    	size_t m_num_threads;
		size_t m_rank;

    	std::vector< std::vector<int> > &m_ref_pvecs;
    	std::vector< std::vector<unsigned int> > &m_ref_aux;
    	std::vector<std::vector <uint8_t> > &m_qs_pvecs;
    	std::map<std::pair<size_t,size_t>,int> &m_vec_res;


    public:
    	worker( size_t num_threads, size_t rank, std::vector< std::vector<int> > &ref_pvecs, std::vector< std::vector<unsigned int> > &ref_aux, std::vector<std::vector <uint8_t> > &qs_pvecs, std::map<std::pair<size_t,size_t>,int> &vec_res )
    	:
          m_num_threads( num_threads ),
          m_rank( rank ),
          m_ref_pvecs( ref_pvecs ),
          m_ref_aux( ref_aux ),
    	  m_qs_pvecs( qs_pvecs ),
    	  m_vec_res( vec_res )
    	{}

    	void pin() {

    		cpu_set_t cpuset;

    		CPU_ZERO(&cpuset);
    		CPU_SET(m_rank, &cpuset);

    		if(pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset) != 0)
			{
    			throw std::runtime_error( "thread pinning failed\n");
			}
    	}

    	void operator()() {

    		const size_t start = W * m_rank;

    		for( size_t i = start; i < m_ref_pvecs.size(); i += (W * m_num_threads) ) {
				const int *seqptrs[W];
				const unsigned int *auxptrs[W];
				size_t num_valid = 0;
				for( size_t j = 0; j < W; ++j ) {
					if( i + j < m_ref_pvecs.size() ) {
						assert(m_ref_pvecs[i+j].size() > 1000 );
						seqptrs[j] = m_ref_pvecs[i+j].data();
						auxptrs[j] = m_ref_aux[i+j].data();
						++num_valid;
					} else {
						seqptrs[j] = seqptrs[0];
						auxptrs[j] = auxptrs[0];
					}
				}

				assert( num_valid > 0 );

				const size_t reflen = m_ref_pvecs[i].size();

//				std::cout << "reflen: " << reflen << "\n";

				const align_pvec_score<score_t,W> aligner(seqptrs, auxptrs, reflen, score_mismatch, score_match_cgap, score_gap_open, score_gap_extend );


	//             std::cerr << "ref: " << i << "\n";
				for( size_t j = 0; j < m_qs_pvecs.size(); ++j ) {

					aligner.align( m_qs_pvecs[j].begin(), m_qs_pvecs[i].end() );
					const score_t *scores = aligner.get_scores();

					for( size_t k = 0; k < num_valid; ++k ) {

			   //         std::cout << "score " << i + k << " " << j << " = " << scores[k] << "\n";
						m_vec_res[std::pair<size_t,size_t>(i+k,j)] = scores[k];
					}
					//ncup += num_valid * reflen * m_qs_pvecs[i].size();
				}


			}

    	}

    };

    template<typename score_t, const size_t W>
	void run_vec_mc( size_t num_threads ) {



    	//typedef worker<typename score_t, Wx> worker_t;


    	timer t1;

    	boost::thread_group tg;

    	std::vector<std::map<std::pair<size_t,size_t>,int> > res(num_threads);
    	if( num_threads > 1 ) {
			for( size_t i = 0; i < num_threads; ++i ) {


				tg.create_thread(worker<score_t, W>(num_threads, i, m_ref_pvecs, m_ref_aux, m_qs_pvecs, res[i]));

			//	tg.create_thread(worker<score_t, W>(num_threads, i, &m_ref_pvecs, &m_ref_aux, &m_qs_pvecs));
			}

			tg.join_all();
    	} else {
    		worker<score_t, W> w(1, 0, m_ref_pvecs, m_ref_aux, m_qs_pvecs, res[0]);
    		w();
    	}
    	m_vec_res.clear();

    	for( size_t i = 0; i < num_threads; ++i ) {
//    		std::cout << res[i].size() << "\n";

    		m_vec_res.insert(res[i].begin(), res[i].end() );
    	}


		std::cerr << "aligned " << m_ref_pvecs.size() << " x " << m_qs_pvecs.size() << " seqs in " << t1.elapsed() << " s\n";
	}


    void compare() {
        
        for( std::map< std::pair< size_t, size_t >, int >::iterator it = m_seq_res.begin(); it != m_seq_res.end(); ++it ) {
            int seq_score = it->second;
            int vec_score = m_vec_res[it->first];
            
            if( seq_score != vec_score ) {
                std::cerr << "error: " << seq_score << " != " << vec_score << "\n";
            }
            
        }
        std::cout << "checked\n";

    }
};


void mini_test() {
    align_arrays arr;
//     const static int score_gap_open = 1;
//     const static int score_gap_extend = 1;
//     const static int score_mismatch = 1;
//     const static int score_match_cgap = 1;    
        
// //     const static int score_gap_open = 3;
// //     const static int score_gap_extend = 1;
// //     const static int score_mismatch = 3;
// //     const static int score_match_cgap = 0;
            
    
    
    const static int score_gap_open = 3;
    const static int score_gap_extend = 1;
    const static int score_mismatch = 3;
    const static int score_match_cgap = 3;
//     
            
    std::vector<int> ref_pvec = { 2, 4, 1, 5, 8, 4, 4, 2, 4, 1, 8, 1, 1, 2, 2 };
    std::vector<uint> ref_aux = {  0, 0, 0, 1, 0, 0, 0, 0 ,0 ,0, 0, 0, 0, 0, 0 };
    
    std::vector<uint8_t> qs_pvec = { 4, 2, 1, 8, 2 };
//     std::vector<uint8_t> qs_pvec = { 1, 4, 8, 1, 4 };
    //                int score2 = 0;
    int score = align_pvec_score_seq( ref_pvec, ref_aux, qs_pvec, score_mismatch, score_match_cgap, score_gap_open, score_gap_extend, arr );
    //std::cout << "score " << i << " " << j << " = " << score << " " << score2 << "\n";
    //                 if( score != score2 ) {
        //                 std::cout << "score " << i << " " << j << " = " << score << "\n";
    //                 }
    
    
}

int main( int argc, char *argv[] ) {

    if( true ) {
        mini_test();
        return 0;
    }
    
	size_t num_refs = atoi(argv[1] );
	size_t num_qs = atoi(argv[2] );
	size_t num_threads = atoi(argv[3] );
	size_t width = atoi(argv[4] );


    testbench tb( "qs.bin", "ref.bin", num_refs, num_qs );
//    tb.run_vec<short,8>();
//    tb.run_vec<int,4>();

    std::cout << "width: " << width << " threads: " << num_threads << " : ";

    if( width == 1 ) {
    	tb.run_seq();
    } else if( width == 4 ) {
    	tb.run_vec_mc<int,4>(num_threads);
    } else if( width == 8 ) {
    	tb.run_vec_mc<short,8>(num_threads);
    } else {
    	std::cerr << "unsupported vector unit width: " << width << "\n";
    }

    
    if( !true ) {
    	tb.run_seq();

    	tb.compare();
    }
}
