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


#include <memory>
#include <algorithm>
#include "ivymike/tree_parser.h"

#include <iostream>

using std::cout;



namespace ivy_mike {
namespace tree_parser_ms {

int adata::s_serial = 0;

lnode *lnode::create( ln_pool &pool ) {
    
    lnode *n = pool.alloc();
    n->next = pool.alloc();
    n->next->next = pool.alloc();
    n->next->next->next = n;
    
    n->m_ldata.reset(pool.get_adata_factory().alloc_ldata());
    n->next->m_ldata.reset(pool.get_adata_factory().alloc_ldata());
    n->next->next->m_ldata.reset(pool.get_adata_factory().alloc_ldata());
    n->m_data.reset( pool.get_adata_factory().alloc_adata() );
    n->next->m_data = n->m_data;
    n->next->next->m_data = n->m_data;
    return n;
    
//     LN *n = new LN();
//     n->next = new LN();
//     n->next->next = new LN();
//     n->next->next->next = n;
//     n->data = n->next->data = n->next->next->data = new AData();
// 
//     return n;

}


void ln_pool::sweep() {

    // TEST: auto-mark pinned roots
    for( std::vector<lnode *>::iterator it = m_pinned_root.begin(); it != m_pinned_root.end(); ++it ) {
        mark(*it);
    }
    
    size_t size1 = m_list.size();

    for ( lt::iterator it = m_list.begin(); it != m_list.end(); ) {
        lt::iterator next = it;
        next++;
        if ( !it->mark ) {
            //lnode *ln = &(*it);
            // TODO: review: is this safe? shouldn't be able to hit any dead lnode (if they unlink themselves properly and the tree is consistent)
            if ( it->back != 0 ) {
                it->back->back = 0;
            }

            m_list.erase(it);
            //delete ln;
            it->dealloc();
        }

        it = next;
    }

    size_t size2 = m_list.size();

   // printf( "sweep: %zd -> %zd\n", size1, size2 );
    std::cerr << "sweep: " << size1 << " -> " << size2 << "\n";
}
void ln_pool::clear() {
    for ( lt::iterator it = m_list.begin(); it != m_list.end(); ++it ) {
        it->mark = false;
    }
}



void ln_pool::mark(lnode* n) {

    n->mark = true;
    n->next->mark = true;
    n->next->next->mark = true;

    if ( n->back != 0 ) {
        if ( !n->back->mark ) {
            mark( n->back );
        }
    }

    if ( n->next->back != 0 ) {
        if ( !n->next->back->mark ) {
            mark( n->next->back );
        }
    }

    if ( n->next->next->back != 0 ) {
        if ( !n->next->next->back->mark ) {
            mark( n->next->next->back );
        }
    }

}


void ln_pool::pin_root( lnode *node ) {
    m_pinned_root.push_back(node);
}

void ln_pool::unpin_root( lnode *node ) {
    const size_t oldsize = m_pinned_root.size();
    
    m_pinned_root.erase( std::remove( m_pinned_root.begin(), m_pinned_root.end(), node ), m_pinned_root.end());
    
    if( oldsize != m_pinned_root.size() + 1 ) {
        std::cerr << oldsize << " -> " << m_pinned_root.size() << std::endl;
        throw std::runtime_error( "ln_pool::unpin_root: root not pinned or duplicate");
    }
}

void parser::readFile(const char* f, std::vector< char >& data) {

    std::ifstream is(f);
    if ( !is.good() ) {
        throw std::runtime_error( "cannot open newick file" );
    }

    is.seekg( 0, std::ios_base::end );
    std::ifstream::off_type size = is.tellg();
    is.seekg( 0, std::ios_base::beg );

    data.resize( size_t(size) );

    is.read( data.data(), size );

}

void parser::substring(const parser::idi_t& from, const parser::idi_t& to, std::string& out) {
//         return new String( Arrays.copyOfRange( inputA, from, to));

    out.clear();

    out.append(from, to );
}



void parser::printLocation() {
    if ( QUIET ) {
        return;
    }

    idi_t pos1 = std::max(inputA.begin(), ptr - 40);
    idi_t pos2 = std::min(inputA.end(), ptr + 40);



    printf( "%s\n", substring(pos1, pos2).c_str());

    for (idi_t i = pos1; i < ptr; ++i) {
        printf(" ");
    }
    printf("^\n");
}

void parser::print_location(std::ostream& os) {
    idi_t pos1 = std::max(inputA.begin(), ptr - 10);
    idi_t pos2 = std::min(inputA.end(), ptr + 10);
    os << "\n";
    os << substring(pos1, pos2).c_str() << "\n";
    for (idi_t i = pos1; i < ptr; ++i) {
        os << " ";
    }
    os << "^\n";
    
}


void parser::skipWhitespace() {
    while ( ptr != inputA.end() && std::isspace(*ptr) ) {
        ++ptr;
    }
    if ( ptr == inputA.end() ) {

        std::stringstream ss;
        
        ss << "hit end of input while skipping white-spaces:\n";
        print_location(ss);
        
        throw std::runtime_error( ss.str() );
    }

}
std::string parser::parseBranchLabel() {
    if ( *ptr == '[' ) {
        idi_t lstart = ptr;
        ++ptr;


        idi_t lend = findNext(ptr, ']' );

        ptr = lend+1;

        // labels have the form [blablabla], so the label content starts at lstart + 1

        if ( lend - (lstart+1) <= 0 ) {
            std::stringstream ss;
            
            ss <<  "bad branch label: " << substring(lstart, ptr) << "\n";
            print_location(ss);
            
            throw std::runtime_error( ss.str() );

        }

        return substring(lstart + 1, lend);


    } else if ( *ptr == '{' ) {
        // alternative label form used in jplace files.
        idi_t lstart = ptr;
        ++ptr;


        idi_t lend = findNext(ptr, '}' );

        ptr = lend+1;

        // alternative labels have the form {blablabla}, so the label content starts at lstart + 1

        if ( lend - (lstart+1) <= 0 ) {
            std::stringstream ss;
            
            ss <<  "bad branch label: " << substring(lstart, ptr) << "\n";
            print_location(ss);
            
            throw std::runtime_error( ss.str() );

        }

        return substring(lstart + 1, lend);


    } else {
        return std::string();
    }


}
lnode* parser::parseNode() {
//         printf( "parseNode\n" );
    skipWhitespace();

    // lookahead: determine node type
    if ( *ptr == '(') {
        return parseInnerNode();
    } else {
        return parseLeaf();
    }
}
parser::idi_t parser::findNext(parser::idi_t pos, char c) {


    while ( pos < inputA.end() && *pos != c) {
        ++pos;
    }

    if ( pos >= inputA.end() ) {
        throw std::runtime_error("reached end of input in find next");

    }
    return pos;
}
bool parser::isFloatChar(char c) {
    return std::isdigit(c) || c == '.' || c == 'e' || c == 'E' || c == '-';
}

std::string parser::substring(const parser::idi_t from, const parser::idi_t to) {
//         return new String( Arrays.copyOfRange( inputA, from, to));



    return std::string(from, to );
}
lnode* parser::parseLeaf() {
//         printf( "parseLeaf\n" );

    skipWhitespace();

    // a leaf consists just of a data string. use the ':' as terminator for now (this is not correct, as there doesn't have to be a branch length (parsr will crash on tree with only one leaf...));
    //int end = findNext(ptr, ':');
    idi_t end = findEndOfBranch(ptr);
    std::string ld = substring(ptr, end);

    ptr = end;


    //      System.out.printf("leaf: %s\n", ld);
    lnode *n = lnode::create( m_pool );
    n->m_data->setTipName(ld);
    n->m_data->setTipSerial(nLeafs);
    //n.data = ld;
    n->m_data->isTip = true; // fake

    nLeafs++;
    return n;
}

double parser::parseBranchLength() {
    skipWhitespace();

    // expect + consume ':'
    if ( *ptr != ':') {
        std::stringstream ss;
        ss <<  "parse error: parseBranchLength expects ':' at " << size_t(ptr - inputA.begin()) << "\n";
        print_location(ss);
        throw std::runtime_error( ss.str() );

    } else {

        ++ptr;

        skipWhitespace();

        idi_t lend = findFloat(ptr);
        if (lend == ptr) {
            std::stringstream ss;
            ss <<  "missing float number at " << size_t(ptr - inputA.begin()) << "\n";
            print_location(ss);
            throw std::runtime_error( ss.str() );
        }

        double l = atof(substring(ptr, lend).c_str());
        ptr = lend;

        return l;
    }
}
parser::idi_t parser::findFloat(parser::idi_t pos) {
    while (isFloatChar(*pos)) {
        ++pos;
    }

    return pos;
}
void parser::twiddle(lnode* n1, lnode* n2, double branchLen, std::string branchLabel, double support) {
    if ( n1->back != 0 ) {
        throw std::runtime_error( "n1.back != null" );
    }

    if ( n2->back != 0 ) {
        throw std::runtime_error("n2.back != null");
    }

    n1->back = n2;
    n2->back = n1;

    n1->backLen = branchLen;
    n2->backLen = branchLen;
    n1->backLabel = branchLabel;
    n2->backLabel = branchLabel;
    n1->backSupport = support;
    n2->backSupport = support;

}

lnode* parser::parseInnerNode() {
//         printf( "parseInnerNode\n" );
    skipWhitespace();


    // expect + consume '('
    if ( *ptr != '(') {
        
        std::stringstream ss;
        ss << "parse error: parseInnerNode expects '(' at " << size_t(ptr - inputA.begin()) << "\n";
        print_location(ss);
        throw std::runtime_error( ss.str() );
        
    }
    ptr++;

    // parse left node + branch length
    lnode *nl = parseNode();
    
    double l1 = 1.0;
    if( *ptr == ':' ) {
        l1 = parseBranchLength();
    }
    
    std::string label1 = parseBranchLabel();

    nl->m_data->setNodeLabel(label1);

    skipWhitespace();


    // expect + consume ','
    if ( *ptr != ',') {
        std::stringstream ss;
        ss << "parse error: parseInnerNode expects ',' at " << size_t(ptr - inputA.begin()) << "\n";
        print_location(ss);
        throw std::runtime_error( ss.str() );
        
    }
    ptr++;


    // parse right node + branch length
    lnode *nr = parseNode();
    double l2 = 1.0;
    if( *ptr == ':' ) {
        l2 = parseBranchLength();
    }
    
    
    std::string label2 = parseBranchLabel();
    nr->m_data->setNodeLabel(label2);

    skipWhitespace();


    nInnerNodes++;
    if ( *ptr == ')') {
        // 'normal' inner node: two childs
        ptr++;

        double support;
        std::string nodeLabel;

        if ( *ptr == ';') {
            // oh my god. a fucking rooted tree
            double sup = std::max( nl->m_data->getSupport(), nr->m_data->getSupport());
            //System.out.printf( "rooted shit: %s %s %f %f %f %f\n", label1, label2, nl.data.getSupport(), nr.data.getSupport(), l1, l2);
            twiddle( nl, nr, l1 + l2, label1, sup );

            return nl;
        }

        if ( *ptr != ':' && *ptr != ',' && *ptr != ')' ) {
            // the stuff between the closing '(' and the ':' of the branch length
            // is interpreted as node-label. If the node label corresponds to a float value
            // it is interpreted as branch support (or node support as a rooted-trees-only-please biologist would say)

            idi_t lend = findEndOfBranch(ptr);

            nodeLabel = substring(ptr, lend);
            ptr = lend;

            bool isDigit = true;
            for ( size_t i = 0; i < nodeLabel.size(); i++ ) {
                isDigit = isDigit && std::isdigit(nodeLabel.at(i));

                if ( i == 0 ) {
                    isDigit = isDigit && (nodeLabel.at(i) != '0');
                }
            }

            if ( isDigit ) {

                support = std::atof(nodeLabel.c_str());
            } else {

                support = -1;
            }


//                int lend = findFloat(ptr);
//                if (lend == ptr) {
//                    printLocation();
//                    throw new RuntimeException("missing float number at " + ptr);
//                }
//
//                support = Double.parseDouble(substring(ptr, lend));
//                ptr = lend;
        } else {
            support = -1.0;
        }

        lnode *n = lnode::create( m_pool );
        n->m_data->setSupport(support);
        n->m_data->setNodeLabel(nodeLabel);

        twiddle( nl, n->next, l1, label1, nl->m_data->getSupport() );
        twiddle( nr, n->next->next, l2, label2, nr->m_data->getSupport() );


        return n;
    } else if ( *ptr == ',') {
        // second comma found: three child nodes == pseudo root
        ptr++;

        lnode *nx = parseNode();

        double l3 = 1.0;
        if( *ptr == ':' ) {
            l3 = parseBranchLength();
        }
        
        std::string label3 = parseBranchLabel();
        //   System.out.printf( "l3: %s\n", nx.data.getTipName() );

        nx->m_data->setNodeLabel(label3);
        skipWhitespace();

        if ( *ptr != ')' ) {
            
            
           std::stringstream ss;
           ss << "parse error: parseInnerNode (at root) expects ') at " << size_t(ptr - inputA.begin()) << "\n";
           print_location(ss);
           throw std::runtime_error( ss.str() );
           
        }
        ptr++;
        skipWhitespace();

        lnode *n = lnode::create( m_pool );

        twiddle( nl, n->next, l1, label1, nl->m_data->getSupport() );
        twiddle( nr, n->next->next, l2, label2, nr->m_data->getSupport() );
        twiddle( nx, n, l3, label3, nx->m_data->getSupport() );

//                      System.out.printf( "root: %f %f %f\n", nl.data.getSupport(), nr.data.getSupport(), nx.data.getSupport() );
//                      System.exit(0);
        return n;
    } else {
         
        std::stringstream ss;
        ss << "parse error: parseInnerNode expects ')'or ',' at " << size_t(ptr - inputA.begin()) << "\n";
        print_location(ss);
        throw std::runtime_error( ss.str() );
        
    }


}
bool parser::isOneOf(char c, const char* cs) {
    while ( *cs != 0 ) {
        if ( c == *cs ) {
            return true;
        }
        ++cs;
    }
    return false;
}
parser::idi_t parser::findEndOfBranch(parser::idi_t pos) {
    const char *termchars = ":,)";


    while ( pos < inputA.end() && !isOneOf( *pos, termchars )) {
        ++pos;
    }
    if ( pos == inputA.end() ) {
        std::stringstream ss;
        ss << "reached end of input while looking for end of branch label\n";
        print_location(ss);
        throw std::runtime_error( ss.str() );
        
    }

    return pos;
}
lnode* parser::parse() {
    nLeafs = 0;
    nInnerNodes = 0;

    skipWhitespace();

    if ( ptr >= inputA.end() ) {
        // seems like we hit the end of the file
        return 0;
    }
    // expect at least one node
    lnode *node = parseNode();

    // expect terminating ';'
    if (ptr >= inputA.end()) {
        throw std::runtime_error("parse error. parse: end of input. missing ';'");
    }

    // branch length might be present for historical reasons
    if ( *ptr == ':' ) {
        parseBranchLength();
    }

    if ( *ptr != ';') {
        std::stringstream ss;
        ss << "parse error. parse expects ';'";
        print_location(ss);
        throw std::runtime_error( ss.str() );
        
    }
    // consume terminating ;
    ++ptr;
    return node;
}

lnode::~lnode() {}
adata::~adata() {}
node_data_factory::~node_data_factory() {}

size_t prune_with_rollback::s_serial_ = 0;

} // namespace tree_parser_ms

} // namespace ivy_mike





// int main() {
// //     getchar();
//     //ivymike::TreeParser tp( "./RAxML_bipartitions.1604.BEST.WITH" );
//     
//     boost::timer t;
//     
//     ivymike::TreeParser tp( "/space/newPNAS/PNAS.ntree" );
//     std::auto_ptr<ivymike::TreeParser::LN> n(tp.parse());
//     
//     //ivymike::LN *n = tp.parse();
//     printf( "n: %f %d\n", n->backLen, n->data->isTip );
// //     getchar();
//     //ivymike::LN::free( n );
// //     delete n;
// //     getchar();
//     
//     cout << t.elapsed() << std::endl;
// }
