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



#include <iterator>
#include <stdexcept>
#include <algorithm>

#include <deque>
#include <boost/dynamic_bitset.hpp>

#include <boost/tr1/unordered_set.hpp>
#include <boost/tr1/unordered_map.hpp>

#include <boost/static_assert.hpp>

#include "ivymike/tree_parser.h"
#include "ivymike/smart_ptr.h"

#include <ivymike/time.h>
#include "ivymike/tree_split_utils.h"
#include "tree_utils.h"



using namespace ivy_mike::tree_parser_ms;
using namespace ivy_mike;

#if 0
//
// generate hash value for boost::bitset. this is a bit verbose, but seems to be the only portable solution
//

template<typename Block>
class bitset_hash_iterator : public std::iterator<std::output_iterator_tag,void,void,void,void> {
    Block &hash;
    size_t i;
public:
    bitset_hash_iterator ( Block &out_hash ) : hash(out_hash), i(1) { hash = 1234; }
    
	inline bitset_hash_iterator<Block>& operator= (const bitset_hash_iterator<Block> &other ) {
        hash = other.hash;
        return *this;
    }

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
        #ifndef WIN32
	// TODO: find out why the asser fails under 64bit win32
	BOOST_STATIC_ASSERT( sizeof( size_t ) == sizeof( dynamic_bitset<>::block_type ) );
        #endif
        dynamic_bitset<>::block_type hash = 0;
        
        to_block_range( bs, bitset_hash_iterator<dynamic_bitset<>::block_type>(hash));
        
        
        
        // FIXME: what to do when size_t and Block have different width?
        // would a 'static if' with no overhead work?
        return size_t(hash);
    }
};


static bool comp_tip_name( const lnode* n1, const lnode* n2 ) {
    return n1->m_data->tipName < n2->m_data->tipName;
}

static void sort_tip_serial( std::vector <lnode *> &tips ) {
    std::sort( tips.begin(), tips.end(), comp_tip_name );
    
    int serial = 0;
    
    for( std::vector< lnode * >::const_iterator it = tips.begin(); it != tips.end(); ++it ) {
        (*it)->m_data->setTipSerial( serial );
        
//         std::cout << (*it)->m_data->tipName << " " << (*it)->m_data->tipSerial << "\n";
        
        serial++;
    }
    
}

void get_all_splits( lnode *t, std::vector< std::pair< lnode*, lnode* > > &edges, std::vector<boost::dynamic_bitset<> > &splits, std::vector<lnode *> &sorted_tips ) {
    
    tip_collector_dumb<lnode> tc;
    visit_lnode( t, tc );
    
//     std::cout << tc.m_nodes.size() << "\n";
    
    sort_tip_serial( tc.m_nodes );
    
    const size_t ntips = tc.m_nodes.size();

    
    std::deque<rooted_bifurcation<lnode> > trav_order;
    rooted_traversal_order(t, t->back, trav_order, false );
    
    std::tr1::unordered_map<int, boost::dynamic_bitset<> > res;
    
    
//     std::cout << "start: " << t->m_data->m_serial << " " << t->back->m_data->m_serial << "\n";
    
    
    // add trivial splits
    for( std::vector< lnode* >::const_iterator it = tc.m_nodes.begin(); it != tc.m_nodes.end(); ++it ) {
        assert( res.find( (*it)->m_data->m_serial ) == res.end() );
        boost::dynamic_bitset<> &bs = res[(*it)->m_data->m_serial];
        bs.resize(ntips);
        bs[(*it)->m_data->tipSerial] = true;
    }
    
    for( std::deque< rooted_bifurcation< ivy_mike::tree_parser_ms::lnode > >::const_iterator it = trav_order.begin(); it != trav_order.end(); ++it ) {
        
//         std::cout << *it << "\n";
        switch( it->tc ) {
        case TIP_TIP:    
        {
            const int pser = it->parent->m_data->m_serial;
            boost::dynamic_bitset<> &bs = res[pser];
            bs.resize(ntips);
            
            bs[it->child1->m_data->tipSerial] = true;
            bs[it->child2->m_data->tipSerial] = true;
            
           
            
            break;
        }
        case TIP_INNER:
        {
            const int pser = it->parent->m_data->m_serial;
            boost::dynamic_bitset<> &bs = res[pser];
            
            const int c2ser = it->child2->m_data->m_serial;
            assert( res.find( c2ser ) != res.end() );
            boost::dynamic_bitset<> &bs_c2 = res[c2ser];
            
            bs = bs_c2;
            assert( !bs[it->child1->m_data->tipSerial] );
            bs[it->child1->m_data->tipSerial] = true;
            
            break;
        }
        case INNER_INNER:
        {
            const int pser = it->parent->m_data->m_serial;
            boost::dynamic_bitset<> &bs = res[pser];
            
            const int c1ser = it->child1->m_data->m_serial;
            const int c2ser = it->child2->m_data->m_serial;
            assert( res.find( c1ser ) != res.end() );
            assert( res.find( c2ser ) != res.end() );
            boost::dynamic_bitset<> &bs_c1 = res[c1ser];
            boost::dynamic_bitset<> &bs_c2 = res[c2ser];
            
            bs = bs_c1;
            bs |= bs_c2;
                        
            //std::cout << "inner inner: " << bs_c1.count() << " " << bs_c2.count() << "\n";
            
            break;
        }
        }
        
    }
    edge_collector<lnode> ec;
    visit_edges(t, ec);

    splits.clear();
    splits.reserve( ec.m_edges.size() );
    
    
    for( std::vector< std::pair< lnode*, lnode* > >::iterator it = ec.m_edges.begin(); it != ec.m_edges.end(); ++it ) {
        
//         std::cout << "edge: " << it->first->m_data->m_serial << " " << it->second->m_data->m_serial << "\n";
        
        assert( res.find( it->first->m_data->m_serial ) != res.end() );
        assert( res.find( it->second->m_data->m_serial ) != res.end() );
        
        boost::dynamic_bitset<> &bs1 = res[it->first->m_data->m_serial];
        boost::dynamic_bitset<> &bs2 = res[it->second->m_data->m_serial];
        
//         std::cout << "count: " << bs1.count() << " " << bs2.count() << "\n";
        
        
        const size_t c1 = bs1.count();
        const size_t c2 = bs2.count();
        const size_t c = std::min( c1, c2 );
        const boost::dynamic_bitset<> &smaller_bs = c1 < c2 ? bs1 : bs2;
        
        splits.push_back( smaller_bs );
        

        // if more than half of the bits are set, flip the bitvector (=make it the smaller split set)
        // if _exactly_ half of the bits are set, modify the vector (=flip it or don't flip it) such that the lowest bit is true (=make it deterministic)
        if( c > ntips / 2 || (ntips % 2 == 0 && c == ntips / 2 && splits.back()[0] == false )) {
            splits.back().flip();
        }
        
//         std::cout << "split: " << splits.back().count() << "\n";
    }   
    
    edges.swap( ec.m_edges );
    sorted_tips.swap(tc.m_nodes);
    
}

inline bool equal_tip_names( const lnode * n1, const lnode * n2 ) {
    assert( n1->m_data->isTip && n2->m_data->isTip );
    
    return n1->m_data->tipName == n2->m_data->tipName;
}


bool split_sets_equal( const std::vector<boost::dynamic_bitset<> > &s1, const std::vector<boost::dynamic_bitset<> > &s2 ) {
    if( s1.size() != s2.size() ) {
        throw std::runtime_error( "split sets have different size" );
    }
    
    std::tr1::unordered_set<boost::dynamic_bitset<>, bitset_hash > m1(s1.begin(), s1.end());
    
    for( std::vector< boost::dynamic_bitset< long unsigned int > >::const_iterator it = s2.begin(); it != s2.end(); ++it ) {
        if( m1.find( *it ) == m1.end() ) {
            return false;
        }
    }
    return true;
}

double compare_trees( lnode *t1, lnode *t2, split_set_t &splits2 ) {

	// TODO: is this valid as a sanity check? If an lnode is deallocated, the m_thisptr is in an undefined state,
	// and there should be tree possible behavious:
	// 1: the memory has not been reused or has been reused for some unrelated stuff,
	//    the internal pointer of m_thisptr is most likely != this => assert
	// 2: the memory has been reused for another lnode => we're out of luck (this situation can not be caught)
	// 3: the memory has become unmapped => reproducible segfault

	assert( t1->m_thisptr.get() == t1 );
	assert( t2->m_thisptr.get() == t2 );

    std::vector<lnode *> sorted_tips;
    std::tr1::unordered_set<boost::dynamic_bitset<>, bitset_hash > split_map;
    {
        
        std::vector< std::pair< lnode*, lnode* > > edges; // found no better way to safely get rid of edges/splits than scoping it...
        std::vector<boost::dynamic_bitset<> > splits;
        get_all_splits( t1, edges, splits, sorted_tips );
        
        
        // TODO: change this to move semantics once rval-refs are supported
        split_map.insert(splits.begin(), splits.end()); 
    }
    
    
    std::vector< std::pair< lnode*, lnode* > > edges2;
    //std::vector<boost::dynamic_bitset<> > splits2;
    splits2.clear();
    std::vector<lnode *> sorted_tips2;
    get_all_splits( t2, edges2, splits2, sorted_tips2 );
    
    size_t nfound = 0;
    for( std::vector< dynamic_bitset<> >::const_iterator it = splits2.begin(); it != splits2.end(); ++it ) {
        //std::cout << "contains: " << (split_map.find( *it ) != split_map.end()) << "\n";
        if(split_map.find( *it ) != split_map.end()) {
            nfound++;
        }
    }
    
//     std::cout << "found: " << nfound << " of " << splits2.size() << "\n";
    
    
    
    return 1.0 - (double(nfound) / splits2.size());
    
}

#endif
int main2( int argc, char *argv[] ) {
//     assert(false);
    ln_pool pool;
    
    if( argc < 3 ) {
        std::cout << "argc < 3\n";
        return -1;
    }
    
    parser p( argv[1], pool );
    lnode *t = p.parse();
    
    
    
    std::vector<lnode *> sorted_tips;
    std::tr1::unordered_set<boost::dynamic_bitset<>, bitset_hash > split_map;
    {
        
        std::vector< std::pair< lnode*, lnode* > > edges; // found no better way to safely get rid of edges/splits than scoping it...
        std::vector<boost::dynamic_bitset<> > splits;
        get_all_splits( t, edges, splits, sorted_tips );
        
        
        // TODO: change this to move semantics once rval-refs are supported
        split_map.insert(splits.begin(), splits.end()); 
    }
    
    
    
    
    
    for( int i = 2; i < argc; i++ ) {
        ivy_mike::timer t1;
        
        parser p2( argv[i], pool );
        lnode *t2 = p2.parse();
        
        std::vector< std::pair< lnode*, lnode* > > edges2;
        std::vector<boost::dynamic_bitset<> > splits2;
        std::vector<lnode *> sorted_tips2;
        get_all_splits( t2, edges2, splits2, sorted_tips2 );
        
        if( sorted_tips.size() != sorted_tips2.size() ) {
            throw std::runtime_error( "trees have different size" );
        }
        if( !std::equal( sorted_tips.begin(), sorted_tips.end(), sorted_tips2.begin(), equal_tip_names ) ) {
            throw std::runtime_error( "tipsets differ" );
        }
        
        
        std::cout << "time1: " << t1.elapsed() << "\n";
        
        size_t nfound = 0;
        for( std::vector< dynamic_bitset<> >::const_iterator it = splits2.begin(); it != splits2.end(); ++it ) {
            //std::cout << "contains: " << (split_map.find( *it ) != split_map.end()) << "\n";
            if(split_map.find( *it ) != split_map.end()) {
                nfound++;
            }
        }
        
        std::cout << "found: " << nfound << " of " << splits2.size() << "\n";
        
        pool.clear();
        pool.mark(t);
        pool.sweep();
        std::cout << "time2: " << t1.elapsed() << "\n";
    }
    
    return 0;

}
