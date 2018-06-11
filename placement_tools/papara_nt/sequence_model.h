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

#ifndef SEQUENCE_MODEL_H_
#define SEQUENCE_MODEL_H_

#include <stdint.h>
#include <cctype>
#include <cstddef>
#include <cassert>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "ivymike/algorithm.h"

#ifndef _MSC_VER
template<typename T>
size_t popcount( T v ) {
	return __builtin_popcount(v);
}
#else
#include <intrin.h>

inline size_t popcount( unsigned short v ) {
	return __popcnt16(v);
}

inline size_t popcount( unsigned int v ) {
	return __popcnt(v);
}

//inline size_t popcount( unsigned __int64 v ) {
//	return __popcnt64(v);
//}


#endif


namespace sequence_model {

class illegal_character : public std::runtime_error {
public:
    illegal_character( const char *msg, int c ) : runtime_error(msg), c_(c) {}
    int c_;
};

class tag_dna;
class tag_dna4;
class tag_aa;


template<typename TAG>
class model {
    //static uint8_t normalize( uint8_t c );
};

template<>
class model<tag_dna> {
public:

    typedef uint8_t pars_state_t;

    const static std::vector<char> inverse_meaning;
//    const static std::vector<pars_state_t> bit_vector;

    static uint8_t normalize( size_t xc ) {
		assert( xc <= 255 );

		int c = int(xc);
        c = std::toupper(c);

        switch( c ) {
        case 'U':
            return 'T';
        case 'N':
        case '?':
        case '.':
            return '-';
        default:
            return c;
        }

    }
    static bool is_known_sstate( size_t c ) {
        return std::find(inverse_meaning.begin(), inverse_meaning.end(), c ) != inverse_meaning.end();
    }

    static pars_state_t s2p( size_t c ) {
        c = normalize(c);
        ptrdiff_t idx = std::distance(inverse_meaning.begin(),
                                   std::find(inverse_meaning.begin(), inverse_meaning.end(), c ) );

        assert( idx >= 0 );

        if( size_t(idx) >= inverse_meaning.size() ) {
            //std::cerr << "illegal character: " << int(c) << "\n";
            throw illegal_character( "illegal character in DNA/RNA sequence", int(c));
        }

        return pars_state_t(idx);
    }

    static uint8_t s2c( size_t c ) {
        c = normalize(c);
        ptrdiff_t idx = std::distance(inverse_meaning.begin(),
                                   std::find(inverse_meaning.begin(), inverse_meaning.end(), c ) );

        assert( idx >= 0 );

        if( size_t(idx) >= inverse_meaning.size() ) {
            throw illegal_character( "illegal character in DNA/RNA sequence", int(c));
        }

        return uint8_t(idx); // safe because of the check above
    }

    static pars_state_t c2p( size_t c ) {
        return pars_state_t(c);
    }


    static uint8_t p2s( pars_state_t c ) {
        return inverse_meaning.at(c);
    }


    static inline bool pstate_is_single(pars_state_t ps) {
        return !pstate_is_gap(ps) && ps != 0;
    }

    static inline bool pstate_is_gap(pars_state_t ps) {
        return ps == gap_pstate();
    }

    static inline bool cstate_is_gap( uint8_t cs) {
        return pstate_is_gap(c2p(cs));
    }
    
    static inline bool cstate_is_single( size_t cs) {
        return pstate_is_single(c2p(cs));
    }
    
    static inline pars_state_t gap_pstate() {
        return pars_state_t(inverse_meaning.size() - 1);
    }
    
    static inline size_t num_cstates() {
        return inverse_meaning.size();
    }

};



template<>
class model<tag_dna4>  {
public:

    typedef uint8_t pars_state_t;

    const static std::vector<char> inverse_meaning;
//    const static std::vector<pars_state_t> bit_vector;

    static uint8_t normalize( size_t xc ) {
                assert( xc <= 255 );

                int c = int(xc);
        c = std::toupper(c);

        switch( c ) {
        case 'U':
            return 'T';
        case 'N':
        case '?':
        case '.':
            return '-';
        default:
            return c;
        }

    }

    
    static bool sstate_is_character( uint8_t c ) {
        c = normalize(c);
        ptrdiff_t idx = std::distance(inverse_meaning.begin(),
                                   std::find(inverse_meaning.begin(), inverse_meaning.end(), c ) );
        
        return idx >= 0 && idx <= gap_cstate();
        
        
    }

    static uint8_t s2c( size_t c ) {
        c = normalize(c);
        ptrdiff_t idx = std::distance(inverse_meaning.begin(),
                                   std::find(inverse_meaning.begin(), inverse_meaning.end(), c ) );

        assert( idx >= 0 );

        if( size_t(idx) >= inverse_meaning.size() ) {
            throw illegal_character( "illegal character in DNA/RNA sequence", int(c));
        }

        return uint8_t(idx); // safe because of the check above
    }

    static uint8_t c2s( size_t c ) {
        return inverse_meaning.at(c);
    }
    

    
    static inline uint8_t gap_cstate() {
        return uint8_t(inverse_meaning.size() - 1);
    }
    
    static inline bool cstate_is_gap( uint8_t cs) {
        return cs == gap_cstate();
    }

    

    static inline size_t num_cstates() {
        return inverse_meaning.size();
    }

};


// FIXME: fake aa model
template<>
class model<tag_aa> {
public:


    typedef uint32_t pars_state_t;
//    const static char inverseMeaningPROT[23];
//    const static unsigned int bitVectorAA[23];

    const static std::vector<char> inverse_meaning;
    const static std::vector<unsigned int> bit_vector;

    static inline uint8_t normalize( size_t xc ) {
		assert( xc <= 255 );

		int c = int(xc);

        return std::toupper(c);
    }

    static bool is_known_sstate( size_t c ) {
        return std::find(inverse_meaning.begin(), inverse_meaning.end(), c ) != inverse_meaning.end();
    }
    
    static pars_state_t s2p( size_t c ) {
        c = normalize(c);
        ptrdiff_t idx = std::distance(inverse_meaning.begin(),
                                   std::find(inverse_meaning.begin(), inverse_meaning.end(), c ) );


//        std::cout << idx << "\n";

        // TODO: is there any reason to use more verbose range checking than '.at'?
        return bit_vector.at(idx);

    }

    static uint8_t s2c( size_t c ) {
        c = normalize(c);
        ptrdiff_t idx = std::distance(inverse_meaning.begin(),
                                   std::find(inverse_meaning.begin(), inverse_meaning.end(), c ) );

        assert( idx >= 0 );

        if( size_t(idx) >= inverse_meaning.size() ) {
            throw illegal_character( "illegal character in DNA/RNA sequence", int(c));
        }

        return uint8_t(idx);
    }

    static pars_state_t c2p( size_t c ) {
        return bit_vector.at(c);
    }

    static uint8_t p2s( pars_state_t c ) {

        ptrdiff_t idx = std::distance(bit_vector.begin(),
                                           std::find(bit_vector.begin(), bit_vector.end(), c ) );



        assert( idx >= 0 );
        if( idx >= ptrdiff_t(inverse_meaning.size()) ) {
            return 'X'; // parsimony state not representable as sequence character.
        }

        return inverse_meaning.at(idx);
    }


    static inline bool pstate_is_single(pars_state_t ps) {
        return popcount(ps) == 1;
    }

    static inline bool pstate_is_gap(pars_state_t ps) {
        return ps == gap_pstate();
    }

    static inline bool cstate_is_gap( size_t cs) {
        
        return pstate_is_gap(c2p(cs));
//         return cs == inverse_meaning.size() - 1;
    }
    
    static inline bool cstate_is_single( size_t cs) {
       return pstate_is_single(c2p(cs));
    }
    static inline pars_state_t gap_pstate() {
        return bit_vector.back(); // by convention the last element of bit_vector is the gap state
    }
    
    static inline size_t num_cstates() {
        return inverse_meaning.size();
    }
};


}



#endif /* SEQUENCE_MODEL_H_ */
