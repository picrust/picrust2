/*
 * Copyright (C) 2013 Simon A. Berger
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
#include <iomanip>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <iterator>

#include "ivymike/large_phylip.h"


void collect_name_prefixes( std::map<std::string, std::vector<std::pair<size_t,size_t> > > &name_prefixes, std::set<std::string> &non_prefixed_names, const ivy_mike::large_phylip &lp, size_t lp_num ) {
    
    const size_t lp_size = lp.size();
    
    for( size_t i = 0; i < lp_size; ++i ) {
        std::string name = lp.name_at(i);
        
        if( name.find("motu_cluster_") == 0 ) { 
            size_t epos = name.find( "." );
            
            assert( epos != name.npos );
            std::string prefix = name.substr(0, epos);
            
            //name_prefixes.insert( prefix );
            std::vector<std::pair<size_t,size_t> > &location_list = name_prefixes[prefix];
            
            location_list.push_back(std::make_pair( lp_num, i ) );
        } else {
            non_prefixed_names.insert( name );
        }
    }
}


std::string merge_nongap_columns( const std::vector<ivy_mike::large_phylip *> &phys, const std::vector<std::pair<size_t,size_t> > &matches ) {
    const size_t seq_len = phys.front()->sequence_len();
    
    std::string out_seq;
    
    
    out_seq.reserve( seq_len );
    
    
    
    for( size_t i = 0; i < seq_len; ++i ) {
        char c = '-';
        
        
        for( std::vector< std::pair< size_t, size_t > >::const_iterator it = matches.begin(); it != matches.end(); ++it ) {
            const ivy_mike::large_phylip &lp = *phys.at( it->first );
            char *base_ptr = (char *)lp.sequence_begin_at( it->second );
            
            if( base_ptr[i] != '-' ) {
                if( c != '-' ) {
                    std::cerr << "meeeeeeeeeep: conflict!" << std::endl;
                    throw std::runtime_error( "bailing out\n" );
                }
                
                c = base_ptr[i];
            }
        
        }
        
        out_seq.push_back( c );
        
    }
    
    return out_seq;
}

int main( int argc, char *argv[] ) {
    std::vector<ivy_mike::large_phylip *> phys;
    std::map<std::string, std::vector<std::pair<size_t,size_t> > > name_prefixes;
    std::set<std::string> non_prefixed_names;
    
    for( size_t i = 1; i < size_t(argc); ++i ) {
        phys.push_back( new ivy_mike::large_phylip(argv[i]));

        collect_name_prefixes( name_prefixes, non_prefixed_names, *phys.back(), phys.size() - 1 );
    }
    
    size_t max_name_len = 0;
    
    std::cout << "non prefix: " << non_prefixed_names.size() << "\n";
    for( std::set< std::string >::iterator it = non_prefixed_names.begin(); it != non_prefixed_names.end(); ++it ) {
        max_name_len = std::max( max_name_len, it->size() );
    }
    
    
    
    for( std::map< std::string, std::vector< std::pair< size_t, size_t > > >::iterator it = name_prefixes.begin(); it != name_prefixes.end(); ++it ) {
        max_name_len = std::max( max_name_len, it->first.size() );
        
        std::cout << it->first << ":";
        for( std::vector< std::pair< size_t, size_t > >::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2 ) {
            std::cout << " " << it2->first << "/" << it2->second;
        }
        
        std::cout << "\n";
    }

    std::ofstream os( "merged.phy" );
    
    
    
    {
        ivy_mike::large_phylip &temp = *phys.front();
        
        os << non_prefixed_names.size() + name_prefixes.size() << " " << temp.sequence_len() << "\n";
    
        for( std::set< std::string >::iterator it = non_prefixed_names.begin(); it != non_prefixed_names.end(); ++it ) {
            size_t idx = temp.getIdx( it->c_str() );
            
            os << std::setw( max_name_len + 1 ) << std::left << *it;
            std::copy( temp.sequence_begin_at(idx), temp.sequence_end_at(idx), std::ostream_iterator<char>(os) );
            os << "\n";
        }
    }
    
    for( std::map< std::string, std::vector< std::pair< size_t, size_t > > >::iterator it = name_prefixes.begin(); it != name_prefixes.end(); ++it ) {
        os << std::setw( max_name_len + 1 ) << std::left << it->first;
        
        std::string merged_seq = merge_nongap_columns( phys, it->second );
        std::copy( merged_seq.begin(), merged_seq.end(), std::ostream_iterator<char>(os) );
        os << "\n";
    }
    
}
