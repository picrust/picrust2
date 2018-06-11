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
#include <iterator>
#include <iomanip>
#include <algorithm>
#include "ivymike/fasta.h"


int main( int argc, char *argv[] ) {
    
    
    
    std::vector<std::string> names;
    std::vector<std::vector<uint8_t> > data;

//     while( std::cin.good() ) {
//         std::cout << char(std::cin.get());
//     }
//     return 0;
    
//     std::cerr << "goodx: " << std::cin.good() << std::endl;
    ivy_mike::read_fasta( std::cin, names, data, false );
    std::cerr << "num: " << data.size() << "\n";
    size_t max_name_len = 0;
    const size_t num_col = data.at(0).size();
    std::vector<size_t> col_nongap_count( num_col );
    for( size_t i = 0; i < names.size(); ++i ) {
        
        const std::vector< uint8_t > &seq = data.at(i);
        assert( seq.size() == num_col ); // TODO: make it a real check!
        
        for( size_t j = 0; j < num_col; ++j ) {
            if( seq[j] != '-' ) {
                ++col_nongap_count[j];
            }
        }
        max_name_len = std::max(names[i].size(), max_name_len);
    }

    const size_t min_nongap = names.size() / 2;
    std::vector<size_t> selected_cols;
    for( size_t i = 0; i < col_nongap_count.size(); ++i ) {
        if( col_nongap_count[i] >= min_nongap ) {
            selected_cols.push_back(i);
        }
    }
    
    
    std::cout << names.size() << " " << selected_cols.size() << "\n";
    
    for( size_t i = 0; i < names.size(); ++i ) {
        std::cout << std::setw(max_name_len + 1) << std::left << names[i];
        //std::copy( data[i].begin(), data[i].end(), std::ostream_iterator<char>(std::cout) );
        for( size_t j = 0; j < selected_cols.size(); ++j ) {
            size_t col = selected_cols[j];
            std::cout << data[i].at(col);
        }
        
        std::cout << "\n";
    }
    
}

