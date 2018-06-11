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

#ifndef __parsimony_h
#define __parsimony_h

#include <algorithm>
#include <stdexcept>
#include <vector>
#include <set>
#include <cassert>
#include <stdint.h>


enum aux_data {
    AUX_CGAP = 0x1,
    AUX_OPEN = 0x2
};




//
// parsimony related stuff
//

typedef int parsimony_state;


// this class is so strange, we could name a new design-pattern after it.
// it does something like a static {} block in java...
class dna_parsimony_mapping_real {
    std::vector<uint8_t> m_p2d;
    std::vector<parsimony_state> m_d2p;
	char my_tolower( char x ) {
		// std::tolower needs std::locale. I will not touch i18n

		if( x >= 'A' && x <= 'Z' ) {
			return x - ('A' - 'a');
		} else {
			return x;
		}
		
	}
public:
    dna_parsimony_mapping_real() : m_p2d(16), m_d2p(256, -1)
    {
        const uint8_t  pd[18] =  {'_', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'U', 'W', 'Y', 'H', 'K', 'D', 'B', '-', 'N'};
        const uint32_t pds[18] = { 0,   1,   2,   3,   4,   5,   6,   7,   8,    8,   9,  10,  11,  12,  13,  14,  15,  15};
        //m_p2d.assign( pd, pd + 16 );


        // this is some weird code, which is supposed tp setup the dna->pstate and pstate->dna maps from ps and pds...
        for( int i = 0; i < 18; i++ ) {
            if( pds[i] >= m_p2d.size() ) {
                throw std::runtime_error( "dna/parsimony map initialization fsck'ed up" );
            }

            m_d2p[pd[i]] = pds[i];
            m_d2p[my_tolower(pd[i])] = pds[i];

            if( m_p2d[pds[i]] == 0 ) {
                m_p2d[pds[i]] = pd[i];
            }
        }

    }

    static dna_parsimony_mapping_real s_pdm;

    static uint8_t p2d( parsimony_state c ) {
        if( c < 0 || c > 15 ) {
            throw std::runtime_error( "illegal parsimony state" );
        }

        return s_pdm.m_p2d[c];
    }

    static parsimony_state d2p( uint8_t c ) {
        if( s_pdm.m_d2p[c] == -1 ) {

//             std::cerr << "illegal: " << c << "\n";
//             throw std::runtime_error( "illegal dna character" );
            return 0xf; // default to undefined pars state
        }

        return s_pdm.m_d2p[c];
    }

    static int d2aux( uint8_t c ) {
        parsimony_state ps = d2p(c);
        if( ps == 0xf ) {
            return AUX_CGAP;
        } else {
            return 0;
        }
    }

    static bool is_gap( uint8_t c ) {
//        return d2p(c) == 0xf;
    	return c == '-';
    }
};




struct dna_parsimony_mapping_simple {
    static parsimony_state d2p( uint8_t c ) {

        switch( c ) {
        case 'A':
        case 'a':
            return 0x1;

        case 'C':
        case 'c':
            return 0x2;

        case 'G':
        case 'g':
            return 0x4;

        case 'U':
        case 'u':
        case 'T':
        case 't':
            return 0x8;

        default:
            return 0xf;
        };
    }

    static uint8_t p2d( parsimony_state c ) {
        switch( c ) {
        case 0x1:
            return 'A';
        case 0x2:
            return 'C';
        case 0x4:
            return 'G';
        case 0x8:
            return 'T';
        case 0xf:
            return '-';

        default:
            return 'X';
        };
    }
    static int d2aux( uint8_t c ) {
        parsimony_state ps = d2p(c);
        if( ps == 0xf ) {
            return AUX_CGAP;
        } else {
            return 0;
        }
    }
};

typedef dna_parsimony_mapping_real dna_parsimony_mapping;


#endif
