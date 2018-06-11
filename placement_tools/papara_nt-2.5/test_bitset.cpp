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


#include <boost/dynamic_bitset.hpp>
#include <vector>
//#include <tr1/unordered_set>
#include <boost/tr1/unordered_set.hpp>
#include <cstdlib>
#include <iostream>
#include <algorithm>

#include <ivymike/time.h>

template<typename Block>
class bitset_hash_iterator : public std::iterator<std::output_iterator_tag,void,void,void,void> {
    Block &hash;
    size_t i;
public:
    bitset_hash_iterator ( Block &out_hash ) : hash(out_hash), i(1) { hash = 1234; }
    
    inline bitset_hash_iterator<Block>& operator= (const Block &v ) {
        hash ^= v * i++;
        return *this;
    }
    
    inline bitset_hash_iterator<Block>& operator* ()
    { return *this; }
    inline bitset_hash_iterator<Block>& operator++ ()
    { return *this; }
    inline bitset_hash_iterator<Block>& operator++ (int)
    { return *this; }
};

class bitset_hash {
public:
    inline size_t operator()( const boost::dynamic_bitset<> &bs ) const {
		// on 64bit windows 'block_type' appears to be 32bit vs. 64bit size_t.
		// I guess the >= in the assert should work, even though the hashing
		// might be suboptimal if only half of the 64bit value is used. don't
		// know what to do. The type chaos (aka. LLP64) on windows sucks.
		
		BOOST_STATIC_ASSERT( sizeof( size_t ) >= sizeof( boost::dynamic_bitset<>::block_type ) );
        
        boost::dynamic_bitset<>::block_type hash = 0;
        
        to_block_range( bs, bitset_hash_iterator<boost::dynamic_bitset<>::block_type>(hash));
        
        return size_t(hash);
    }
};

int main() {
    
    std::vector<boost::dynamic_bitset<> > temp;
    std::vector<boost::dynamic_bitset<> > temp_search;
    
    
    for( int i = 0; i < 1024; i++ ) {
        temp.push_back(boost::dynamic_bitset<>());
        for( int j = 0; j < 1024; j++ ) {
            int b = rand() >= RAND_MAX / 2;
            
            temp.back().push_back(b);
        }
        
    }
    
    temp_search = temp;
    std::random_shuffle( temp_search.begin(), temp_search.end());
    for( int i = 0; i < 1024; i++ ) {
        temp_search.push_back(boost::dynamic_bitset<>());
        for( int j = 0; j < 1024; j++ ) {
            int b = rand() >= RAND_MAX / 2;
            
            temp_search.back().push_back(b);
        }
        
    }
//     temp_search = temp;
    
    size_t found_vec = 0;
    size_t found_set = 0;
    
    
    ivy_mike::perf_timer pt;
    
    std::vector<boost::dynamic_bitset<> > svec = temp;
    
    
    std::sort(svec.begin(), svec.end());
    
    pt.add_int();
    for( int i = 0; i < 1000; ++i ) {
    
        for( std::vector< boost::dynamic_bitset< long unsigned int > >::iterator it = temp_search.begin(); it != temp_search.end(); ++it ) {
            bool found = std::binary_search( svec.begin(), svec.end(), *it );
            
            if( found ) {
                found_vec++;
            }
        }
    }
    pt.add_int();
    std::tr1::unordered_set<boost::dynamic_bitset<>, bitset_hash> sset( temp.begin(), temp.end() );
    pt.add_int();
    
    for( int i = 0; i < 1000; ++i ) {
        for( std::vector< boost::dynamic_bitset< long unsigned int > >::iterator it = temp_search.begin(); it != temp_search.end(); ++it ) {
            bool found = sset.find(*it) != sset.end();
            if( found ) {
                found_set++;
                
            }
        }
    }
    pt.add_int();
    pt.print();
    std::cout << svec.size() << " " << sset.size() << "\n";
    std::cout << found_vec << " " << found_set << "\n";
    
}