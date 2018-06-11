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

#ifndef __ivymike__tree_parser_h
#define __ivymike__tree_parser_h


#include "compiler_capabilities.h"

#include <vector>
#include <fstream>
#include <stdexcept>
#include <cstdlib>
#include <sstream>
#include <cctype>
#include <cassert>
#include <iostream>
#include <memory>
#include "demangle.h"


#include <boost/intrusive/slist.hpp>
#include <list>
namespace ivy_mike {






namespace tree_parser_ms {
//using namespace boost;
    
struct adata {

    bool isTip;
    std::string tipName;
    void setTipName( const std::string &name ) {
        tipName = name;
    }

    int tipSerial;
    void setTipSerial( int serial ) {
        tipSerial = serial;
    }

    double support;
    double getSupport() {
        return support;
    }
    void setSupport(double sup) {

        support = sup;
    }

    //
    // CHANGE: 26.10.2011: nodeLabel now stores the string after branch length (the
    // same string is also used for the branch label.
    // Previously the string directly after the closing ')' of a inner node has been used,
    // which only works for inner nodes. This will now go into extra_node_label_
    // CHECK IF THIS BREAKS ANYTHING!
    //
    std::string nodeLabel;
    std::string extra_node_label_;
    void setNodeLabel(const std::string &nl) {
        nodeLabel = nl;
    }

    void set_extra_node_label( const std::string &l ) {
    	extra_node_label_ = l;
    }

    static int s_serial;
    const int m_serial;
    
    adata() : isTip( false ), tipSerial(-1), support(-1), m_serial(s_serial++) {}

    virtual ~adata() ;
    
    virtual void visit() {}

    virtual void print( std::ostream &os ) const {
        if( isTip ) {
        
            os << "(TIP:" << m_serial << ": " << tipName.c_str() << ")";
        } else {
            os << "(INNER:" << m_serial << ")";
        }    
    }
    
    template<typename T>
    T* get_as() {
        if( !ivy_mike::isa<T>(*this) ) {
            std::cerr << "lnode::adata::get_as: trying incompatible dynamic type conversion:\n"
                << ivy_mike::demangle(typeid(*this).name())
                << "\ninto\n"
                << ivy_mike::demangle(typeid(T).name()) << "\n";

            throw std::runtime_error( "dynamic cast failed\n" );
        }
        
        return static_cast<T*>(this);
    }

};
    
inline std::ostream &operator<<(std::ostream &os, const adata &ad ) {
    ad.print( os );
    return os;
}


struct ldata {
    virtual ~ldata() ;
};


class ln_pool;

    
class lnode
    : public boost::intrusive::slist_base_hook<> 

{

    lnode( const lnode &other ) {}
    const lnode &operator=(const lnode &other ) { return *this; }

public:
    static lnode *create( ln_pool &pool );

    
    
    lnode() : next(0), back(0), backLen(-1), backSupport(-1), mark(false), towards_root(false) {
        m_thisptr.reset(this);
    }
    virtual ~lnode() ;
    
    std::weak_ptr<lnode>get_smart_ptr() {
        return std::weak_ptr<lnode>(m_thisptr);
    }
    
    std::shared_ptr<lnode>dealloc() {
        std::shared_ptr<lnode>tmp(m_thisptr);
        m_thisptr.reset();
        return tmp;
    }
    std::shared_ptr<adata> m_data;
    std::shared_ptr<ldata> m_ldata;
    //LN *next;
    lnode *next;
    lnode *back;

    double backLen;
    double backSupport;
    std::string backLabel;

    bool mark;
    
    std::shared_ptr<lnode> m_thisptr;
    
    bool towards_root;
    
    
};



class node_data_factory {

public:
    virtual adata *alloc_adata() {
        return new adata();
    }
    
    
    virtual ldata *alloc_ldata() {
        return 0;
    }

    virtual ~node_data_factory() ;
};

class ln_pool 
{
     
  
public:
#if IVY_MIKE__USE_CPP11 
	typedef std::unique_ptr<node_data_factory> fact_ptr_type;
    // this version takes ownership of fact!
//     ln_pool( fact_ptr_type fact ) : m_ad_fact(std::move(fact)) {}
    ln_pool( fact_ptr_type fact ) : m_ad_fact(std::move(fact)) {}
#else
    typedef std::auto_ptr<node_data_factory> fact_ptr_type;
    // this version takes ownership of fact!
    ln_pool( fact_ptr_type fact ) : m_ad_fact(fact) {}
#endif
    
    
    //ln_pool( std::shared_ptr<node_data_factory> fact ) : m_ad_fact(fact) {}
    ln_pool() : m_ad_fact(new node_data_factory) {}
    
    ~ln_pool() {
        clear();
        sweep();
    }
    node_data_factory &get_adata_factory() {
        return *m_ad_fact;
    }
    
    lnode *alloc() {
        lnode *n = new lnode();
        m_list.push_front( *n );
        return &m_list.front();
        
        //return new LN();
    }
    
    void mark( lnode *n ) ;
    
    void sweep() ;
    void clear() ;
    
    void pin_root( lnode *n );
    void unpin_root( lnode *n );
    
private:
    typedef boost::intrusive::slist<lnode> lt;
    std::vector<lnode *> m_pinned_root;

    lt m_list;

    //std::shared_ptr<node_data_factory> m_ad_fact;
    
    
    
    fact_ptr_type m_ad_fact;
};


class ln_pool_pin {
	lnode * m_pin;
	ln_pool &m_pool;

public:
	ln_pool_pin( lnode *node, ln_pool &pool ) : m_pin( node ), m_pool(pool) {
		assert( m_pin != 0 );
		m_pool.pin_root(m_pin);
	}

	~ln_pool_pin() {
		if( m_pin != 0 ) {
			release();
		}
	}

	void release() {
		if( m_pin == 0 ) {
			throw std::runtime_error( "ln_pool_pin::release called with m_pin == 0" );
		}
		m_pool.unpin_root(m_pin);
		m_pin = 0;
	}


};

class parser {
public:

    //String input;
private:
    // input as char array
    typedef std::vector<char> id_t;
    id_t inputA;

    typedef id_t::iterator idi_t;
    idi_t ptr;


    int nLeafs;
    int nInnerNodes;
    ln_pool &m_pool;

    bool QUIET;
//      public TreeParser(String input) {
//      //      this.input = input;
//              this.inputA = input.toCharArray();
//              ptr = 0;
//      }


    void readFile(const char *f, std::vector<char> &data) ;
    void substring( const idi_t &from, const idi_t &to, std::string &out ) ;

    std::string substring( const idi_t from, const idi_t to ) ;

    void printLocation();
    void print_location( std::ostream &os ) ;
    

    void skipWhitespace() ;



    std::string parseBranchLabel() ;

    lnode *parseNode() ;
    idi_t findNext( idi_t pos, char c) ;

    bool isFloatChar(char c) ;

    idi_t findFloat(idi_t pos) ;

    double parseBranchLength() ;

    lnode *parseLeaf() ;

    idi_t findEndOfBranch(idi_t pos) ;

    bool isOneOf(char c, const char *cs) ;
    lnode *parseInnerNode() ;


public:
     /**
     * create an edge (=double link) between the two nodes and set branch length
     *
     * @param n1
     * @param n2
     * @param branchLen
     */
    static void twiddle( lnode *n1, lnode *n2, double branchLen, std::string branchLabel, double support ) ;

    
    
    parser(const char *f, ln_pool &pool ) : m_pool(pool) {
        QUIET = false;
        readFile(f, inputA);
        //  this.input = new String(inputA);
        ptr = inputA.begin();
    }
    
    template<typename iiter>
    parser( iiter first, iiter last, ln_pool &pool ) : m_pool(pool) {
        QUIET = false;
        
        inputA.assign( first, last );
        ptr = inputA.begin();
    }
    
    /**
     * Call this after object construction to parse the complete tree,
     * @return pseudo root of the tree
     */

    lnode *parse() ;

};

inline void print_newick( lnode *node, std::ostream &os, bool root = true ) {
    if( node->m_data->isTip ) {
        //         assert( !root ); // MEEEEP: unhandled case: tip given as starting node! Can't be bothered now to fix it.
        if( root ) {
            os << "(";
        }


        os << node->m_data->tipName;
        if( root ) {
            os << ",";
            print_newick(node->back, os, false);
            os << ");";
            
        } else {
            os << ":" << node->backLen;
            if( !node->backLabel.empty() ) {
                os << "[" << node->backLabel << "]";
            }
        }
    } else {
    	assert( node->next->back != 0 );
    	assert( node->next->next->back != 0 );
        os << "(";
        print_newick(node->next->back, os, false);
        os << ",";
        print_newick(node->next->next->back, os, false);
        if( root ) {
            os << ",";
            print_newick(node->back, os, false);
            os << ");";
        } else {
            os << "):" << node->backLen;
            
            if( !node->backLabel.empty() ) {
                os << "[" << node->backLabel << "]";
            }
        }
    }
}

inline lnode *towards_tree( lnode *node ) {
    if( !node->m_data->isTip || node->back != 0 ) {
        return node;
    } else if( node->next->back != 0 ) {
        return node->next;
    } else if( node->next->next->back != 0 ) {
        return node->next->next;
    } else {
        throw std::runtime_error( "cannot find valid back pointer." );
    }
}

inline lnode *next_non_tip( lnode *node ) {
    if( !node->m_data->isTip ) {
        return node;
    } else {
        if( node->back->m_data->isTip ) {
            throw std::runtime_error( "cannot find non-tip node: two connected tips." );
        }
        return node->back;
    }
}


// 'transaction like' insert subtree starting at 'c' into the edge between e and e->back
// unless 'commit' is called, 'c' will be removed on destruction, restoring the tree
// (hopefully) to the same state as it was before insertion.
//
// c must be a 'tip' in a tree with at least two nodes, and c->next and c->next->next
// are used to link into the other tree.

class splice_with_rollback {
	double m_len;
	double m_support;
	std::string m_label;

	lnode *m_e1;
	lnode *m_e2;
	lnode *m_c;

	bool m_commit;
public:
	splice_with_rollback( lnode *e, lnode *c ) : m_commit(false) {
		m_len = e->backLen;
		m_support = e->backSupport;
		m_label = e->backLabel;

		m_e1 = e;
		m_e2 = e->back;
		m_c = c;

		m_e1->back = 0;
		m_e2->back = 0;

		assert( c->next->back == 0 );
		assert( c->next->next->back == 0 );

		parser::twiddle(m_e1, m_c->next, m_len / 2, "" /*"tmp1"*/, 0 );
		parser::twiddle(m_e2, m_c->next->next, m_len / 2, "" /*"tmp2"*/, 0 );
	}

	~splice_with_rollback() {

		if( !m_commit ) {
			m_e1->back = 0;
			m_e2->back = 0;
			m_c->next->back = 0;
			m_c->next->next->back = 0;

			parser::twiddle( m_e1, m_e2, m_len, m_label, m_support );
		}
	}

	void commit() {
		m_commit = true;
	}

};

class prune_with_rollback {
public:
#if IVY_MIKE__USE_CPP11

    prune_with_rollback( prune_with_rollback && other ) : 
    serial_( other.serial_ ),
    commit_( other.commit_ ),
    
    len1_(other.len1_),
    len2_(other.len2_),
    
    label1_(other.label1_),
    label2_(other.label2_),
    
    support1_(other.support1_),
    support2_(other.support2_),

    n_(other.n_),
    back1_(other.back1_),
    back2_(other.back2_),
    valid_(other.valid_)
    {
        // invalidate donor
        other.valid_ = false;
    }
#endif
    
    prune_with_rollback( lnode *n ) : serial_(s_serial_++), commit_(false), n_(n) {
//         std::cout << "pwr: " << serial_ << "\n";
        
        if( n->next->back == 0 || n->next->next->back == 0 ) {
            throw std::runtime_error( "trying to prune tip or unlinked node" );
        }
            
        assert( n != 0 );
        back1_ = n->next->back;
        back2_ = n->next->next->back;
        
        len1_ = n->next->backLen;
        len2_ = n->next->next->backLen;
        
        label1_ = n->next->backLabel;
        label2_ = n->next->next->backLabel;
        
        support1_ = n->next->backSupport;
        support2_ = n->next->next->backSupport;
        
       
        n->next->back = 0;
        n->next->next->back = 0;
        
        back1_->back = back2_;
        back2_->back = back1_;
        
        double new_len = back1_->backLen + back2_->backLen;
        std::string back_label = back1_->backLabel + back2_->backLabel;
        back1_->backLen = back2_->backLen = new_len;
        back1_->backLabel = back2_->backLabel = back_label;
        back1_->backSupport = back2_->backSupport = 0.0;
        
        valid_ = true;
    }
    
    ~prune_with_rollback() {
        // do nothing for 'moved-from' pwrs
        if( !valid_ ) {
            return;
        }
        
//         std::cout << "~pwr: " << serial_ << "\n";
        if( !commit_ ) {
            // std::cerr << "WARNING: untested code: ~prune_with_rollback!!!!\n"; // it seems to work in the spr_vis_test code, what else could possibly go wrong?
            
            n_->next->back = back1_;
            n_->next->next->back = back2_;
        
            back1_->back = n_->next;
            back2_->back = n_->next->next;
            
            n_->next->back->backLen = n_->next->backLen = len1_;
            n_->next->next->back->backLen = n_->next->next->backLen = len2_;
        
            n_->next->back->backLabel = n_->next->backLabel = label1_;
            n_->next->next->back->backLabel = n_->next->next->backLabel = label2_;
        
            n_->next->back->backSupport = n_->next->backSupport = support1_;
            n_->next->next->back->backSupport = n_->next->next->backSupport = support2_;
        }
    }
    
    void commit() {
        commit_ = true;
    }
    
    lnode *get_save_node() const {
        // return lnode that is guaranteed to be still part of the tree after the pruning.
        // WARNING: make sure that this is always a node next to the pruned node (=member of the new branch),
        // because other stuff depends on that (e.g., stepwise addition with destiny tree)
        return back1_;
    }
    
private:
    
    prune_with_rollback();
    prune_with_rollback(prune_with_rollback &);
    prune_with_rollback & operator=(prune_with_rollback &);
    
    static size_t s_serial_;
    size_t serial_;
    
    
    bool commit_;
    
    double len1_;
    double len2_;
    std::string label1_;
    std::string label2_; 
    double support1_;
    double support2_;
    
    lnode *n_;
    
    lnode *back1_;
    lnode *back2_;
    
    bool valid_;
};
// forward to stupid static method (why did I put it in a static method? bad java habit?)
inline void twiddle_nodes( lnode *n1, lnode *n2, double branchLen, std::string branchLabel, double support ) {
    parser::twiddle( n1, n2, branchLen, branchLabel, support );

}

} // namespace tree_parser_ms

}
#endif
