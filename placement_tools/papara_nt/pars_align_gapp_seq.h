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
#include <iomanip>
#include <cstdio>
#include <cstddef>

class pars_align_gapp_seq {
	typedef unsigned int pars_state_t;

    typedef double score_t;
    const static score_t LARGE_VALUE;// = 32000; //  score_t must be able to keep LARGE_VALUE + GAP_OPEN + GAP_EXTEND without overflowing!. waiting for c++0x to be able to use std::numeric_limits at compile-time...
    const score_t GAP_OPEN; //= 1;
    const score_t GAP_EXTEND; // = 1;
    const score_t GAP_OPEN_EXTEND;
    const score_t MISMATCH_PENALTY;// = 4;
    const score_t MATCH_CGAP;
public:

    bool g_dump;
    enum {
        BT_STAY = 0x1,
        BT_STAY_L = 0x2,
        BT_UP = 0x4
    };

    class arrays {
    public:

        size_t m_s;

        std::vector<score_t> score;
        std::vector<score_t> scoreL;
        std::vector<uint8_t> dir;


        bool m_haveDir;



        arrays( bool haveDir = false ) : m_s(0), m_haveDir(haveDir) {

        }



        void size( size_t s ) {
            if ( m_s < s ) {
                std::cout << "arr resize: " << m_s << " -> " << s << "\n";

                m_s = s;
                score.resize(m_s);
                scoreL.resize(m_s);

                if ( m_haveDir ) {
                    dir.resize(m_s);
                }
            }
        }
        
        size_t size() {
            return m_s;
        }
        
        
        
    };
    size_t m_ncups;
private:
    size_t m_na;
    size_t m_nb;

    size_t m_ma;
    size_t m_mb;

    size_t m_msize;

    const int *m_a;
    const unsigned char *m_b;
    const double *m_gapp;

    size_t m_aStride;


    arrays &m_arr;
    ptrdiff_t m_tbStartA;
    ptrdiff_t m_tbStartB;

    const unsigned int *m_bvtrans;

    inline size_t addr( size_t a, size_t b ) {
#if 1
        return a + b * m_ma;
#else
        return b + a * m_mb;
#endif
    }

    inline size_t saddr( ptrdiff_t a, ptrdiff_t b ) {
        return addr( a + 1, (b + 1));
    }

//     inline int xsaddr( int a, int b ) {
//      return addr( a + 1, (b + 1) % 2 );
//     }

    inline pars_state_t getSeqA ( size_t a ) {
        return pars_state_t(m_a[a * m_aStride]);
    }


    inline double get_gapp ( size_t a ) {
        return m_gapp[a];
    }

public:

    pars_align_gapp_seq( const int* seqA, const unsigned char* seqB, size_t n_a, size_t n_b, size_t aStride, const double *a_gapp, size_t aAuxStride, arrays &arr, const unsigned int *bvtrans = 0,
               score_t gapOpen = 1, score_t gapExtend = 1, score_t mismatch = 3, score_t matchCGap = 10 )
    // : LARGE_VALUE(std::numeric_limits<score_t>::max() - 100)
    //: LARGE_VALUE(32000),
    ;

    virtual ~pars_align_gapp_seq() ;



//      void alignFreeshiftS1 () {
//
//
//
//          for ( int ia = 0; ia <= m_na - m_nb - 1; ia++ ) {
//              m_arr->scoreL[saddr ( ia, -1 ) ] = 0;
//              m_arr->score[saddr ( ia, -1 ) ] = 0;
//          }
//
//
//          m_arr->score[saddr ( -1, -1 ) ] = 0;
//          m_arr->scoreL[saddr ( -1, -1 ) ] = 0;
//
//          // NM += pa->lenB * (pa->lenA - pa->lenB);
//
//
//      }

    void alignFreeshiftS11 () ;



    




    
    void tracebackCompressed( std::vector< uint8_t>& bvec ) ;


    inline int alignFreeshiftS3Sep(int highCutoff) {
        
#define DUMP_MATRIX 1

        bool g_dump = !true;

#define AUX_CGAP ( 0x1 )
        score_t best = LARGE_VALUE;
        const size_t band_width = m_na - m_nb;
        score_t * __restrict sp = m_arr.score.data();
        score_t * __restrict sLp = m_arr.scoreL.data();
        uint8_t * __restrict dir = m_arr.dir.data();

//        assert( dir != 0 );
        const bool fill_dir = m_arr.m_haveDir;


        if( fill_dir ) {
			for ( size_t ib = 0; ib < m_nb; ib++ ) {
				pars_state_t cb;// = m_b[ib];



				if ( m_bvtrans == 0 ) {
					cb = m_b[ib];
				} else {
					cb = m_bvtrans[m_b[ib]];
				}

				size_t astart = ib;
				score_t last_sp = sp[saddr(astart - 1, ib)];
				score_t last_sLp = sLp[saddr(astart - 1, ib)];



				//          int saddr_0_0 = saddr( astart, ib ); // score address of the current cell
				//             int saddr_1_1 = saddr( astart-1, ib-1 ); // score address of the upper-left-diagonal cell

				// /*       score_t *sp_0_0 = &sp[saddr_0_0];
				//          score_t *sp_1_1 = &sp[saddr_1_1];*/

				for ( size_t ia = astart; ia <= ib + band_width; ia++ /*, saddr_0_0++, saddr_1_1++, sp_0_0++, sp_1_1++*/ ) {
					const pars_state_t ca = getSeqA( ia );


					const double p_nongap = get_gapp( ia );
					const double p_gap = 1 - p_nongap;

					// determine match or mis-match according to parsimony bits coming from the tree.
					const bool match = ( ca & cb ) != 0;
					const size_t addr = saddr( ia, ib );


					// calculate match score ('go diagonal')
					// penalize based on parsimony and gap probability

					score_t sd = sp[saddr( ia-1, ib-1 )];
					score_t su = sp[saddr( ia, ib-1 )];


					su += GAP_OPEN_EXTEND;
					if ( !match ) {
						sd += MISMATCH_PENALTY;
					}

					score_t scoreExt = last_sLp + p_nongap * GAP_EXTEND;
					score_t scoreOpen = last_sp + p_nongap * GAP_OPEN_EXTEND;

					sd += MATCH_CGAP * p_gap;
					//su += MATCH_CGAP * p_gap;

					if ( scoreExt < scoreOpen ) {
						if( fill_dir ) {
							dir[addr] = BT_STAY_L;
						}
						last_sLp = scoreExt;
					} else {

						if( fill_dir ) {
							dir[addr] = 0;
						}
						last_sLp = scoreOpen;
					}


					int dsd = 0;
					if ( su < sd ) {
						sd = su;
						//dir[addr] |= BT_UP;
						dsd = BT_UP;
					} else {
						dsd = BT_STAY;
					}

					if ( last_sLp < sd ) {
						sp[addr] = last_sp = last_sLp;
					} else {
						sp[addr] = last_sp = sd;

						if( fill_dir ) {
							dir[addr] |= dsd;
						}
					}

					//              rowmin = std::min( last_sp, rowmin );
					//          ct++;
					//              m_ncups++;
				}


				//              if ( rowmin >= highCutoff ) {
				//                  return INT_MAX;
				//              }
			}

			if( g_dump ) {
				for( int ib = -1; ib < int(m_nb); ib++ ) {

					for( int ia = -1; ia < int(m_na); ia++ ) {

						//printf( "\t%d", sp[saddr(ia, ib)] );

						int d = 0;
						if( fill_dir ) {
							d = int(dir[saddr(ia,ib)]);
						}

						std::cout << std::setw(10) << sp[saddr(ia, ib)] << "," << d;


					}
					std::cout << "\n";



				}

			}

			//      printf( "ct: %zd\n", ct );
			//      exit(0);
			m_tbStartB = m_nb - 1;
			m_tbStartA = -1;

			for ( size_t a = m_na - 1; a >= m_nb - 1; a-- ) {


				score_t s = m_arr.score[saddr ( a, m_tbStartB ) ];
				if ( s < best ) {
					best = s;
					m_tbStartA = a;
				}
			}





			return int(m_arr.score[saddr ( m_tbStartA, m_tbStartB ) ]);

        } else {
        	// back-port to sequential from align_pvec_vec.h

        	std::fill( sp, sp + m_na, 0 );
        	sp[band_width+1] = LARGE_VALUE;


        	for( size_t ib = 0; ib < m_nb; ib++ ) {
				pars_state_t bc;

				if ( m_bvtrans == 0 ) {
					bc = m_b[ib];
				} else {
					bc = m_bvtrans[m_b[ib]];
				}

				double last_sl = LARGE_VALUE;
				double last_sc = LARGE_VALUE;
				double last_sdiag = 0;

				size_t astart = ib;

				double * __restrict s_iter = &sp[0];
				double * __restrict s_iter_next = &sp[1];

				last_sdiag = *s_iter;


				for( size_t ia = astart; ia <= ib + band_width; ++ia, ++s_iter, ++s_iter_next ) {
					const pars_state_t ac = getSeqA(ia);


					const double p_nongap = get_gapp( ia );
					const double p_gap = 1 - p_nongap;


					const bool match = ( ac & bc ) != 0;

					const double sm = ((!match) ? last_sdiag + MISMATCH_PENALTY : last_sdiag) + p_gap * MATCH_CGAP;


					const double sl_ext = last_sl + p_nongap * GAP_EXTEND;
					const double sl_open = last_sc + p_nongap * GAP_OPEN_EXTEND;
					last_sl = std::min( sl_ext, sl_open );

					const double min_sm_sl = std::min( sm, last_sl );

					last_sdiag = *s_iter_next;
					const double su = last_sdiag + GAP_OPEN_EXTEND;
					const double sc = std::min( min_sm_sl, su );

					last_sc = sc;

					*s_iter = sc;

				}
			}
        	double minscore = LARGE_VALUE;
			for( size_t i = 0; i < band_width + 1; ++i ) {
				minscore = std::min( minscore, sp[i]);
			}

			return int(minscore);
//			for ( size_t ib = 0; ib < m_nb; ib++ ) {
//				pars_state_t cb;// = m_b[ib];
//
//				if ( m_bvtrans == 0 ) {
//					cb = m_b[ib];
//				} else {
//					cb = m_bvtrans[m_b[ib]];
//				}
//
//				size_t astart = ib;
//				score_t last_sp = sp[saddr(astart - 1, ib)];
//				score_t last_sLp = sLp[saddr(astart - 1, ib)];
//
//				for ( size_t ia = astart; ia <= ib + band_width; ia++ /*, saddr_0_0++, saddr_1_1++, sp_0_0++, sp_1_1++*/ ) {
//					const pars_state_t ca = getSeqA( ia );
//
//					const double p_nongap = get_gapp( ia );
//					const double p_gap = 1 - p_nongap;
//
//					// determine match or mis-match according to parsimony bits coming from the tree.
//					const bool match = ( ca & cb ) != 0;
//					const size_t addr = saddr( ia, ib );
//
//					// calculate match score ('go diagonal')
//					// penalize based on parsimony and gap probability
//
//					score_t sd = sp[saddr( ia-1, ib-1 )];
//					score_t su = sp[saddr( ia, ib-1 )];
//
//					su += GAP_OPEN_EXTEND;
//					if ( !match ) {
//						sd += MISMATCH_PENALTY;
//					}
//
//					score_t scoreExt = last_sLp + p_nongap * GAP_EXTEND;
//					score_t scoreOpen = last_sp + p_nongap * GAP_OPEN_EXTEND;
//
//					sd += MATCH_CGAP * p_gap;
//					//su += MATCH_CGAP * p_gap;
//
//					last_sLp = std::min( scoreExt, scoreOpen);
//
//
//					int dsd = 0;
//
//					sd = std::min( sd, su );
//					sp[addr] = last_sp = std::min( last_sLp, sd );
//
//				}
//
//			}

        }
        

    }
    inline int alignFreeshift(int highCutoff) {
        alignFreeshiftS11();


        //     return alignFreeshiftS3(highCutoff);
        return alignFreeshiftS3Sep(highCutoff);
    }  
    
};
