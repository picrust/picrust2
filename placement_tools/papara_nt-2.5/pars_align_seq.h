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


#ifndef __pars_align_seq_h
#define __pars_align_seq_h
#include <iostream>
#include <cstdio>
#include <cstddef>


template<typename qs_state_t = uint8_t>
class pars_align_seq {
    typedef unsigned int pars_state_t;

    typedef short score_t;
    const static score_t LARGE_VALUE = 32000; //  score_t must be able to keep LARGE_VALUE + GAP_OPEN + GAP_EXTEND without overflowing!. waiting for c++0x to be able to use std::numeric_limits at compile-time...
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
    	arrays( const arrays & ) {}
    	const arrays& operator=(const arrays &) {return *this;}


    public:

        size_t m_s;
        score_t *score;
        score_t *scoreL;

        bool m_haveDir;
        uint8_t *dir;


        arrays( bool haveDir = false ) : m_s(0), score(0), scoreL(0), m_haveDir(haveDir), dir(0) {

        }

        ~arrays() {
            delete [] score;
            delete [] scoreL;
            delete [] dir;
        }

        void size( size_t s ) {
            if ( m_s < s ) {
                std::cout << "arr resize: " << m_s << " -> " << s << "\n";
                
                delete [] score;
                delete [] scoreL;


                m_s = s;
                score = new score_t[m_s];
                scoreL = new score_t[m_s];

                if ( m_haveDir ) {
                    delete [] dir;
                    dir = new uint8_t[m_s];
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
    const qs_state_t *m_b;
    const unsigned int *m_aAux;

    size_t m_aStride;
    size_t m_aAuxStride;

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


    inline int getAuxA ( size_t a ) {
        return m_aAux[a * m_aAuxStride];
    }

    // FIXME: a private copy constructor would still have to initialize the const data members. it's not worth the effort
    // wait for c++0x explicit constructor delete...
    //pars_align_seq( const pars_align_seq &other ) {}
    const pars_align_seq &operator=( const pars_align_seq &other ) { return *this; }

public:

    pars_align_seq( const int* seqA, const qs_state_t* seqB, size_t n_a, size_t n_b, size_t aStride, const unsigned int *aAux, size_t aAuxStride, arrays &arr, const unsigned int *bvtrans = 0,
               score_t gapOpen = 1, score_t gapExtend = 1, score_t mismatch = 3, score_t matchCGap = 10 )
          : GAP_OPEN(gapOpen), GAP_EXTEND(gapExtend), GAP_OPEN_EXTEND(
                   GAP_OPEN + GAP_EXTEND), MISMATCH_PENALTY(mismatch), MATCH_CGAP(
                   matchCGap), g_dump(false), m_ncups(0), m_arr(arr), m_tbStartA(
                   -1), m_tbStartB(-1), m_bvtrans(bvtrans)

   {
       if (n_a <= n_b) {
           throw std::runtime_error("n_a <= n_b\n");
       }

       m_na = n_a;
       m_nb = n_b;

       m_ma = n_a + 1;
       m_mb = n_b + 1;

       m_a = seqA;
       m_b = seqB;

       m_aAux = aAux;
       m_aStride = aStride;
       m_aAuxStride = aAuxStride;

       m_msize = m_ma * m_mb;

       m_arr.size(m_msize);

   }
    // : LARGE_VALUE(std::numeric_limits<score_t>::max() - 100)
    //: LARGE_VALUE(32000),
    ;

    virtual ~pars_align_seq() {}



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

    void alignFreeshiftS11() {

        if (m_na > m_nb) {
            for (size_t ia = 0; ia <= m_na - m_nb - 1; ia++) {
                //  std::cout << "init: " << ia << " " << m_na << " " << m_nb << "\n";

                //size_t xxx = saddr ( ia, -1 );
                //         if( xxx >= m_arr.size() ) {
                //             std::cout << "meeeep: " << xxx << " " << ia << "\n";
                //         }
                m_arr.scoreL[saddr(ia, -1)] = 0;
                m_arr.score[saddr(ia, -1)] = 0;
            }
        }

        m_arr.score[saddr(-1, -1)] = 0;
        m_arr.scoreL[saddr(-1, -1)] = 0;

        for (size_t i = 0; i < m_nb; i++) {
            m_arr.score[saddr(i - 1, i)] = LARGE_VALUE;
            m_arr.scoreL[saddr(i - 1, i)] = LARGE_VALUE;
            m_arr.scoreL[saddr(i, i)] = LARGE_VALUE;
            m_arr.score[saddr((m_na - m_nb) + i, i - 1)] = LARGE_VALUE;
        }

        for (size_t i = 0; i < m_na; i++) {
            m_arr.score[saddr(i, -1)] = 0;
        }

        if (m_arr.dir != 0) {

            if (m_na > m_nb) {
                for (size_t ia = 0; ia <= m_na - m_nb - 1; ia++) {
                    m_arr.dir[saddr(ia, -1)] = BT_STAY_L;
                }
            }
        }
    }




    




    
    void tracebackCompressed(std::vector<uint8_t>& bvec) {
        //     int idx_b = m_na-1;
        ptrdiff_t ba = m_tbStartA;
        ptrdiff_t bb = m_tbStartB;
        bool inL = false;

        if (m_tbStartA < ptrdiff_t(m_na) - 1) {
            for (ptrdiff_t ba = m_na - 1; ba > m_tbStartA; ba--) {
                bvec.push_back(1);

            }
        }

        // traceback

        while ((ba >= 0 || bb >= 0)) {
            if (!inL) {

                // we are currently not in a gap. see if we should match ('go diagonal') or open a gap ('jump to gap-matrix')
                if ((m_arr.dir[saddr(ba, bb)] & BT_STAY) != 0) {
                    bvec.push_back(0);
    //                 idx_b--;
                    ba--;
                    bb--;
                } else if ((m_arr.dir[saddr(ba, bb)] & BT_UP) != 0) {
                    bb--;
                    bvec.push_back(2);

                } else {
                    // open gap. keep position and switch to gap-matrix
                    inL = true;
                }
            } else {
                // we are in a gap

                // output gap
                bvec.push_back(1);

    //             idx_b--;
                // check if we should stay in the gap-matrix
                inL = (m_arr.dir[saddr(ba, bb)] & BT_STAY_L) != 0;
                ba--;
            }
        }

    //     for( int i = 0; i <  bvec.size(); i++ ) {
    //         printf( "tb %d: %d\n", i, bvec[i] );
    //     }

    }


    inline int alignFreeshiftS3Sep(int highCutoff) {
        
#define DUMP_MATRIX 1

        bool g_dump = !true;

        const static int AUX_CGAP = 0x1;
        score_t best = LARGE_VALUE;
        const size_t band_width = m_na - m_nb;
        score_t * __restrict sp = m_arr.score;
        score_t * __restrict sLp = m_arr.scoreL;
        uint8_t * __restrict dir = m_arr.dir;

        assert( dir != 0 );


        for ( size_t ib = 0; ib < m_nb; ib++ ) {
            pars_state_t cb;// = m_b[ib];


#ifdef DUMP_MATRIX
            if ( g_dump ) {
                std::cout << ib << " ";

            }
#endif

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

                // cgap == true means that we are at a position with a 'constant gap', according to the tree phylogeny (this is
                // related to the gap flagging in prank+F).
                // At such a position matches are penalized and gaps are free.

                const bool cgap = getAuxA( ia ) == AUX_CGAP;
                //const bool cgap = false; // FIXME: vector version debugging!
                //const bool cgap = get_aux_a( pa, ia ) == AUX_CGAP;

                // determine match or mis-match according to parsimony bits coming from the tree.
                const bool match = ( ca & cb ) != 0;
                const size_t addr = saddr( ia, ib );


                // calculate match score ('go diagonal')
                // penalize based on parsimony and fixed gaps

                score_t sd = sp[saddr( ia-1, ib-1 )];
                score_t su = sp[saddr( ia, ib-1 )];

                
                su += MISMATCH_PENALTY;
//                      su = LARGE_VALUE;
                if ( !match ) {
                    sd += MISMATCH_PENALTY;
                }

                score_t scoreExt;
                score_t scoreOpen;

                if ( cgap ) {
                    scoreExt = last_sLp;
                    scoreOpen = last_sp;
//                         scoreExt = last_sLp + GAP_EXTEND;
//                         scoreOpen = last_sp + GAP_OPEN_EXTEND;
                    sd += MATCH_CGAP;
                } else {
                    scoreExt = last_sLp + GAP_EXTEND;
                    scoreOpen = last_sp + GAP_OPEN_EXTEND;
                }


                if ( scoreExt < scoreOpen ) {
                    dir[addr] = BT_STAY_L;
                    last_sLp = scoreExt;
                } else {
                    dir[addr] = 0;
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
                    dir[addr] |= dsd;
                }
#ifdef DUMP_MATRIX
                if ( g_dump ) {

                    printf( "%d(%d,%d) ", last_sp, dir[addr], sd );
                    
                }
#endif


                //              rowmin = std::min( last_sp, rowmin );
                //          ct++;
                //              m_ncups++;
            }

#ifdef DUMP_MATRIX
            if ( g_dump ) {
                printf( "\n" );
            }
#endif


            //              if ( rowmin >= highCutoff ) {
            //                  return INT_MAX;
            //              }
        }
        
        
        if( g_dump ) {
            for( int ib = -1; ib < int(m_nb); ib++ ) {
                
                for( int ia = -1; ia < int(m_na); ia++ ) {
                        
                    printf( "\t%d", sp[saddr(ia, ib)] );
                    
                    
                }
                printf( "\n" );
                
                
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





        return m_arr.score[saddr ( m_tbStartA, m_tbStartB ) ];
    }
    inline int alignFreeshift(int highCutoff) {
        alignFreeshiftS11();


        //     return alignFreeshiftS3(highCutoff);
        return alignFreeshiftS3Sep(highCutoff);
    }  
    
};
#endif
