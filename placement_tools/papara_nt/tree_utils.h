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


#ifndef __tree_utils_h
#define __tree_utils_h
#if 0
#include <algorithm>
#include <stdexcept>
#include <vector>
#include <set>
#include <cassert>
#include <stdint.h>
#include <boost/tr1/unordered_set.hpp>
#include "ivymike/smart_ptr.h"
#include "ivymike/tree_parser.h"

enum tip_case {
    TIP_TIP,
    TIP_INNER,
    INNER_INNER
    
};



template<class lnode>
struct rooted_bifurcation {
    
    lnode * parent;
    lnode * child1;
    lnode * child2;
    
        
    tip_case tc;
    
    rooted_bifurcation() : parent(0), child1(0), child2(0) {}
    
    rooted_bifurcation( lnode *p, lnode *c1, lnode *c2, tip_case t ) 
        : parent(p),
        child1(c1),
        child2(c2),
        tc(t)
    {}
};


template<typename lnode>
inline bool operator==( const rooted_bifurcation<lnode> &n1, const rooted_bifurcation<lnode> &n2 ) {
	return n1.parent == n2.parent && n1.child1 == n2.child1 && n1.child2 == n2.child2 && n1.tc == n2.tc;
}

template<class lnode>
inline std::ostream &operator<<( std::ostream &os, const rooted_bifurcation<lnode> &rb ) {
    const char *tc = "W.T.F.";
    
    switch( rb.tc ) {
    case TIP_TIP:
        tc = "TIP_TIP";
        break;
        
    case TIP_INNER:
        tc = "TIP_INNER";
        break;
        
	case INNER_INNER:
        tc = "INNER_INNER";
        break;
    }
    
    os << tc << " " << *(rb.parent->m_data) << " " << *(rb.child1->m_data) << " " << *(rb.child2->m_data);
    return os;
}

template <class lnode, class container>
void rooted_traveral_order_rec( lnode *n, container &cont, bool incremental = false ) {
    
    
    
    lnode *n1 = n->next->back;
    lnode *n2 = n->next->next->back;
    
    assert( n->m_data->isTip || n1 != 0 );
    assert( n->m_data->isTip || n2 != 0 );
    
    // FIXME: why did I allow n to be a tip in the assertions above. fix this some time!
    // this function will crash if n is a tip...
    // TODO: should this function silently ignore+return if n is a tip, or
    // is it better to enforce strict handling at a higher level?
    // This function itself does not descent into links to tips, so high level handling
    // makes more sense.

    assert( n1 != 0 && n2 != 0 );

    // STUPID: this just bit me in a non-debug build (assertions disabled), so make it a real runtime error (but keep in assertion, because it is better for debugging)
    if( n1 == 0 || n2 == 0 ) {
        throw std::runtime_error( "n1 == 0 || n2 == 0" );
    }
    
    n->towards_root = true;
    n->next->towards_root = false;
    n->next->next->towards_root = false;
    
    
    
    if( n1->m_data->isTip && n2->m_data->isTip ) {
        cont.push_front( rooted_bifurcation<lnode>( n, n1, n2, TIP_TIP ));
    } else if( n1->m_data->isTip && !n2->m_data->isTip ) {
        cont.push_front( rooted_bifurcation<lnode>( n, n1, n2, TIP_INNER ));
        
        if( !incremental || !n2->towards_root ) {
            rooted_traveral_order_rec( n2, cont );
        }
    } else if( !n1->m_data->isTip && n2->m_data->isTip ) {
        cont.push_front( rooted_bifurcation<lnode>( n, n2, n1, TIP_INNER ));
        
        if( !incremental || !n1->towards_root ) {
            rooted_traveral_order_rec( n1, cont );    
        }
        
    } else {
        cont.push_front( rooted_bifurcation<lnode>( n, n1, n2, INNER_INNER ));
        
        if( !incremental || !n1->towards_root ) {
            rooted_traveral_order_rec( n1, cont );
        }
        if( !incremental || !n2->towards_root ) {
            rooted_traveral_order_rec( n2, cont );
        }
    }
}

//template <class lnode, class container>
//void rooted_traveral_order( lnode *n1, lnode *n2, container &cont, bool incremental ) {
//
//    if( !n1->m_data->isTip ) {
//        rooted_traveral_order_rec<lnode, container>( n1, cont, incremental );
//    }
//    if( !n2->m_data->isTip ) {
//        rooted_traveral_order_rec<lnode, container>( n2, cont, incremental );
//    }
//
//
//    //std::reverse( cont.begin(), cont.end());
//}




template <class lnode, class container>
void rooted_traversal_order( lnode *n1, lnode *n2, lnode *n3, container &cont, bool incremental ) {

    if( !n1->m_data->isTip ) {
        rooted_traveral_order_rec<lnode, container>( n1, cont, incremental );
    }
    if( !n2->m_data->isTip ) {
        rooted_traveral_order_rec<lnode, container>( n2, cont, incremental );
    }
    if( n3 != 0 ) {
    	if( !n3->m_data->isTip ) {
    		rooted_traveral_order_rec<lnode, container>( n3, cont, incremental );
    	}
    }


    //std::reverse( cont.begin(), cont.end());
}


template <class lnode, class container>
void rooted_traversal_order( lnode *n1, lnode *n2, container &cont, bool incremental ) {
	rooted_traversal_order<lnode,container>( n1, n2, 0, cont, incremental );

}


template <class lnode, typename oiter>
void rooted_preorder_traversal( lnode *n, oiter start, bool incremental = false ) {



	if( n == 0 ) {
		return;
	}

	if( n->m_data->isTip ) {
		return;
	}

    lnode *n1 = n->next->back;
    lnode *n2 = n->next->next->back;



    assert( n->m_data->isTip || n1 != 0 );
    assert( n->m_data->isTip || n2 != 0 );

    n->towards_root = true;
    n->next->towards_root = false;
    n->next->next->towards_root = false;



    if( n1->m_data->isTip && n2->m_data->isTip ) {
        *(start++) = rooted_bifurcation<lnode>( n, n1, n2, TIP_TIP );
    } else if( n1->m_data->isTip && !n2->m_data->isTip ) {
    	*(start++) = rooted_bifurcation<lnode>( n, n1, n2, TIP_INNER );

        if( !incremental || !n2->towards_root ) {
            rooted_preorder_traversal( n2, start );
        }
    } else if( !n1->m_data->isTip && n2->m_data->isTip ) {
    	*(start++) = rooted_bifurcation<lnode>( n, n2, n1, TIP_INNER );

        if( !incremental || !n1->towards_root ) {
            rooted_preorder_traversal( n1, start );
        }

    } else {
    	*(start++) = rooted_bifurcation<lnode>( n, n1, n2, INNER_INNER );

        if( !incremental || !n1->towards_root ) {
            rooted_preorder_traversal( n1, start );
        }
        if( !incremental || !n2->towards_root ) {
            rooted_preorder_traversal( n2, start );
        }
    }
}





template <class lnode>
lnode *towards_tree( lnode *n ) {
    int ct = 0;
    
    while( n->back == 0 ) {
        n = n->next;
        
        if( ct > 3 ) {
        
            throw std::runtime_error( "node not connected to tree" );
        }
        
        ct++;
    }
    
    return n;
    
}




template <class visitor>
void visit_lnode( typename visitor::lnode *n, visitor &v, bool go_back = true ) {
    v( n );
    
    if( go_back && n->back != 0 ) {
        visit_lnode( n->back, v, false );
    }
    if( n->next->back != 0 ) {
        visit_lnode( n->next->back, v, false );   
    }

    if( n->next->next->back != 0 ) {
        visit_lnode( n->next->next->back, v, false );
    }
};

template <class visitor>
void visit_lnode_postorder( typename visitor::lnode *n, visitor &v, bool go_back = true ) {

    if( go_back && n->back != 0 ) {
        visit_lnode( n->back, v, false );
    }
    if( n->next->back != 0 ) {
        visit_lnode( n->next->back, v, false );   
    }

    if( n->next->next->back != 0 ) {
        visit_lnode( n->next->next->back, v, false );
    }
    
    v( n );
    
};

template <class LNODE, class CONT = std::vector<std::shared_ptr<LNODE> > >
struct tip_collector {
    typedef LNODE lnode;
    typedef CONT container;
    
    //container<lnode *> m_nodes;
  
    container m_nodes;
    

    void operator()( lnode *n ) {
        if( n->m_data->isTip ) {
            m_nodes.push_back(n->get_smart_ptr().lock());
        }
    }
};


template <class LNODE>
struct tip_collector_dumb {
    typedef LNODE lnode;
    
    
    //container<lnode *> m_nodes;
  
    std::vector<lnode *> m_nodes;
    

    void operator()( lnode *n ) {
        if( n->m_data->isTip ) {
            m_nodes.push_back(n);
        }
    }
};

//template <class LNODE>
//struct node_collector_dumb {
//    typedef LNODE lnode;
//    typedef std::vector<lnode *> container;
//    container m_nodes;
//
//
//    void operator()( lnode *n ) {
//        if( n->m_data->isTip ) {
//            m_nodes.push_back(n);
//        }
//    }
//};



template <typename visitor>
void visit_edges( typename visitor::lnode *n, visitor &v, bool at_root = true ) {
    assert( n->back != 0 );
    
    if( !at_root ) {
        v( n, n->back );
    }
    
    if( at_root && n->back != 0 ) {
        visit_edges( n->back, v, false );
    }
    if( n->next->back != 0 ) {
        visit_edges( n->next->back, v, false );   
    }

    
    if( n->next->next->back != 0 ) {
        visit_edges( n->next->next->back, v, false );
    }
    // at the root, the edge between n and n->back will be visited when recursing to n->back
      
};

template <class LNODE>
struct edge_collector {
public:
	typedef LNODE lnode;
    typedef std::pair<LNODE *, LNODE *> edge;
    typedef std::vector<edge> container;

    void operator()( lnode *n1, lnode *n2 ) {
//         std::cout << "edge: " << n1 << " " << n2 << "\n";
        m_edges.push_back( edge( n1, n2 ) );
        
    }

    container m_edges;
};


class node_level_assignment {
	typedef ivy_mike::tree_parser_ms::lnode lnode;

    std::vector<std::pair<int,lnode *> > m_level_mapping;

    std::tr1::unordered_set<lnode *>m_mix;
    std::tr1::unordered_set<lnode *>m_closed;




    size_t round( int level ) {

//        std::cerr << "round " << level << " " << m_mix.size() << "\n";

        std::tr1::unordered_set<lnode *> cand;

        std::vector<lnode *>rm;


        for( std::tr1::unordered_set<lnode *>::iterator it = m_mix.begin(); it != m_mix.end(); ++it ) {
        	lnode *n = *it;
//        	std::cout << "mix: " << level << " " << n << " " << n->next << "\n";

        	assert( !n->m_data->isTip );

        	if( m_mix.find(n->next) != m_mix.end() ) {
//        		std::cout << "level: " << level << " " << *(n->next->next->m_data) << "\n";

//        		m_level_mapping.push_back( std::pair<int,lnode*>(level, n->next->next ));
        		rm.push_back( n );
//        		rm.push_back( n->next );
//        		if( m_closed.find( n->next->next->back ) == m_closed.end() ) {
//        			m_mix.insert( n->next->next->back );
//        		}
        	}

        }

        for( std::vector<lnode *>::iterator it = rm.begin(); it != rm.end(); ++it ) {
        	lnode *n = *it;
        	m_level_mapping.push_back( std::pair<int,lnode*>(level, n->next->next ));
        	m_closed.insert(n->next->next);

        	m_mix.erase( *it );
        	m_mix.erase( (*it)->next );
        }
        for( std::vector<lnode *>::iterator it = rm.begin(); it != rm.end(); ++it ) {
        	lnode *n = *it;
        	if( m_closed.find( n->next->next->back ) == m_closed.end() ) {
        		m_mix.insert( n->next->next->back );
        	}
        }

    	return m_mix.size();
    }

public:
    node_level_assignment( std::vector<lnode *> tips ) {

        for( std::vector<lnode *>::iterator it = tips.begin(); it != tips.end(); ++it ) {
            m_level_mapping.push_back( std::pair<int,lnode*>( 0, *it ) );
//            m_closed.insert( *it );

            m_closed.insert(*it);
            assert( (*it)->back != 0 );

            if( !(*it)->back->m_data->isTip ) {
            	m_mix.insert( (*it)->back );
            }

        }

        int level = 1;

        while( round( level++ ) != 0 ) {}


//        for( std::vector<std::pair<int,lnode *> >::iterator it = m_level_mapping.begin(); it != m_level_mapping.end(); ++it ) {
//        	std::cout << "level: " << it->first << " " << *(it->second->m_data) << "\n";
//        }

//        std::cout << "end\n";
    }

    std::vector<std::pair<int,lnode *> > &get_level_mapping() {
    	return m_level_mapping;

    }
};


//////////////////////////////////////////////////////////////////////////
// STL inspired replacements to the 'visitor' stuff above.
// (Which actually has nothing to do with the stupid java hack that is
// commonly called visitor pattern)

// STL compatible function to iterate over the lnodes (kind of like a std::copy)
template <typename lnode, typename oiter>
void iterate_lnode( lnode *n, oiter start, bool go_back = true ) {
    *(start++) = n;
    //outer++;

    if( go_back && n->back != 0 ) {
        iterate_lnode( n->back, start, false );
    }
    if( n->next->back != 0 ) {
    	iterate_lnode( n->next->back, start, false );
    }

    if( n->next->next->back != 0 ) {
    	iterate_lnode( n->next->next->back, start, false );
    }
};

template <typename lnode, typename Tfunc>
void apply_lnode( lnode *n, Tfunc func, bool go_back = true ) {
    //*(start++) = n;
    func(n);
    
    //outer++;

    if( go_back && n->back != 0 ) {
        apply_lnode( n->back, func, false );
    }
    if( n->next->back != 0 ) {
        apply_lnode( n->next->back, func, false );
    }

    if( n->next->next->back != 0 ) {
        apply_lnode( n->next->next->back, func, false );
    }
};

// TEST" implement interate_node in terms of the more generic lnode_apply
template <typename lnode, typename oiter>
void iterate_lnode_test( lnode *n, oiter start, bool go_back = true ) {
    // call recurseve apply with a 'iterator inserter' lambda
    apply_lnode(n, [&start](lnode *n) {*(start++)=n;}, go_back );
};



// UNTESTED: back_insert_iterator that only only inserts if predicate is true (TDD is for pussies)
// combined with is_tip and iterate_lnodes it replaces most of the stupid hard-coded collectors from above.

template <class Container, typename Pred>
class back_insert_if_iterator :
    public std::iterator<std::output_iterator_tag,void,void,void,void>
{
protected:
	Container* container_;
	Pred p_;
public:
	typedef Container container_type;
	explicit back_insert_if_iterator (Container& x, Pred p ) : container_(&x), p_(p) {}
	back_insert_if_iterator<Container, Pred>& operator= (typename Container::const_reference value) {


		if( p_( value ) ) {
			container_->push_back(value);
		}
		return *this;
	}
	back_insert_if_iterator<Container,Pred>& operator* () {
		return *this;
	}
	back_insert_if_iterator<Container,Pred>& operator++ () {
		return *this;
	}
	back_insert_if_iterator<Container,Pred> operator++ (int) {
		return *this;
	}
};
template <class Container, typename Pred>
back_insert_if_iterator<Container, Pred> back_insert_ifer(Container &cont, Pred p ) {
	return back_insert_if_iterator<Container, Pred>(cont, p);
}

inline bool is_tip( ivy_mike::tree_parser_ms::lnode * l ) {
	assert( l != 0 );
	assert( l->m_data != 0 );
	return l->m_data->isTip;
}






#endif
#endif
