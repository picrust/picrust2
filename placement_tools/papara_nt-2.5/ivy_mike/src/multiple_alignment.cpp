/*
 * Copyright (C) 2009-2012 Simon A. Berger
 * 
 * This file is part of ivy_mike.
 * 
 *  ivy_mike is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ivy_mike is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with ivy_mike.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <fstream>
#include <cassert>
#include <stdexcept>
#include <iostream>
#include <cctype>
#include <algorithm>
#include <iterator>

#include "ivymike/multiple_alignment.h"

using namespace ivy_mike;

bool multiple_alignment::load_phylip( std::istream &is ) {
    size_t nTaxon;
    size_t seqLen;
    
    if( !is.good() ) {
        throw std::runtime_error( "cannot read phylip file" );
    }

    is >> nTaxon;
    is >> seqLen;

    while( isspace(is.peek() )) { is.get(); }
    
    
//     printf( "header: %zd %zd\n", nTaxon, seqLen );

    size_t n = 0;


    for( size_t i = 0; i < nTaxon; ++i ) {
        std::string line;
        std::getline( is, line ); 
        
        std::string::iterator ws_it = std::find_if( line.begin(), line.end(), isspace );
        
        if( std::distance( line.begin(), ws_it ) == 0 ) {
            throw std::runtime_error( "could not read taxon name in phylip file: line starts with whitespace\n" );
        }
        
        if( std::distance( ws_it, line.end() ) == 0 ) {
            throw std::runtime_error( "could not read taxon name in phylip file: no whitespace before end of line\n" );
        }
        
        
        std::string name( line.begin(), ws_it );
        
        std::vector<uint8_t> seq;
        std::remove_copy_if( ws_it, line.end(), std::back_inserter(seq), isspace );
        
        names.push_back(name);
        data.push_back( seq );
        
//         printf( "name: %s\n", name.c_str() );
        n++;
    }

    while( data.front().size() < seqLen ) {
        {
            // skip one empty line
            std::string line;
            std::getline( is, line ); 
        }
        assert( data.size() == nTaxon );
        for( size_t i = 0; i < nTaxon; ++i ) {
            std::string line;
            std::getline( is, line ); 
            
            if( !is.good() ) {
                throw std::runtime_error( "early end of (interleaved) phylip file." );
                
            }
            
            std::string::iterator ws_it = std::find_if( line.begin(), line.end(), isspace );
            
            std::vector<uint8_t> seq;
            std::remove_copy_if( line.begin(), line.end(), std::back_inserter(data.at(i)), isspace );

        }
        
    }
    
    assert( n == nTaxon );
//     printf( "n: %zd\n", n );

    return n == nTaxon;
}



// bool multiple_alignment::load_phylip( std::istream &is ) {
//     size_t nTaxon;
//     size_t seqLen;
//     
//     if( !is.good() ) {
//         throw std::runtime_error( "cannot read phylip file" );
//     }
// 
//     is >> nTaxon;
//     is >> seqLen;
// 
// //     printf( "header: %zd %zd\n", nTaxon, seqLen );
// 
//     size_t n = 0;
// 
// 
//     while ( !is.eof() ) {
//         std::string name;
//         std::string seq;
// 
//         is >> name;
//         is >> seq;
// 
//         if( is.eof() ) {
//             
//             break;
//         }
// //         std::cout << "name: " << name << "\n";
//         if( seq.size() != seqLen ) {
//             std::cerr << "name: " << name << "\n";
//             std::cerr << "data: " << seq.size() << "\n";
//             
//             throw std::runtime_error( "bad sequence in phylip file\n" );
//         }
//         
//         names.push_back(name);
//         data.push_back( std::vector<uint8_t>(seq.begin(), seq.end()));
//         
//         
// 		
// 		
// //         printf( "name: %s\n", name.c_str() );
//         n++;
//     }
// 
//     assert( n == nTaxon );
// //     printf( "n: %zd\n", n );
// 
//     return n == nTaxon;
// }

bool multiple_alignment::load_phylip( const char *name ) {
	std::ifstream is( name );
	
	if( is.good() ) {
		return load_phylip( is );
	} else {
		//throw std::runtime_error( "cannot open phylip file" );
		return false;
	}
	
}
