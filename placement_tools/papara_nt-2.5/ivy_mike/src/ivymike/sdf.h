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

//
// read and write sdf/mdl (like) files, mostly according to the format described in
// (Arthur Dalby et al., J. Chem. Inf. Comput. Sci, 1992, 32, 244-255) 
// pdf: /http://pubs.acs.org/doi/pdf/10.1021/ci00007a012
//


#ifndef __ivymike_sdf
#define __ivymike_sdf
#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <cstring>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cstdlib>
#include <stdint.h>
#include <cmath>
#include "boost/array.hpp"
#include "boost/io/ios_state.hpp"

namespace ivy_mike {


struct sdf_int_full {
    typedef sdf_int_full sdf_int;
    struct atom {
        float   m_x, m_y, m_z;
//        char    m_ele;

        typedef  boost::array<char,3> ele_t;
        ele_t m_ele;
        int     m_can_ele;
        //std::string     m_extra;

        atom( float x, float y, float z, char *ele, int can_ele ) : m_x(x), m_y(y), m_z(z), m_can_ele(can_ele)
        {
            m_ele = to_ele( ele );

            //      printf( "atom: %f %f %f %c\n", x,y,z,ele);
        }

        atom( float x, float y, float z, ele_t ele, int can_ele ) : m_x(x), m_y(y), m_z(z), m_ele(ele), m_can_ele(can_ele)
        {


            //      printf( "atom: %f %f %f %c\n", x,y,z,ele);
        }

        static inline ele_t to_ele( char *ele ) {

            ele_t r;

            size_t nns = 0;
            for ( size_t i = 0; i < ele_t::size(); i++ ) {

                if ( !isspace(ele[i] ) ) {
                    nns++;

                }
            }

            if ( nns < 1 || nns > 2 ) {
                std::stringstream ss;
                ss << "nns < 1 || nns > 2: '" << ele[0] << ";" << ele[1] << ";" << ele[2] << "'";


                throw std::runtime_error( ss.str() );
            }



            std::copy( ele, ele + ele_t::size(), r.begin() );
            return r;
        }
        void print_mdl( std::ostream &os ) const {
                  
            if( std::fabs( m_x ) > 9999 || std::fabs( m_y ) > 9999 || std::fabs( m_z ) > 9999 ) {
                std::cerr << "xyz: " << m_x << " " << m_y << " " << m_z << "\n";
                throw std::runtime_error( "molecule coordinate out of bounds" );
            }
            boost::io::ios_all_saver ias( os );
            os << std::fixed << std::setprecision(4) <<  std::setw(10) <<  m_x <<  std::setw(10) << m_y << std::setw(10) << m_z <<  m_ele[0] << m_ele[1] << m_ele[2] << "  0  0  0  0  0  0  0  0  0  0  0  0\n";
//             os << "    0.0000    0.0000    0.0000" <<  m_ele[0] << m_ele[1] << m_ele[2] << "  0  0  0  0  0  0  0  0  0  0  0  0\n";
        }

    };
    struct bond {
        std::pair <int,int>     m_atoms;
        int                     m_type;
        int                     m_can_type;
        // std::string             m_extra;


        bond( int first, int second, int type, int can_type ) : m_atoms(first,second), m_type(type), m_can_type(can_type)
        {
            //    printf( "bond: %d %d %d\n", m_atoms.first, m_atoms.second, m_type );
        }
        void print_mdl( std::ostream &os ) const {
            os << std::setw(3) << int(m_atoms.first) << std::setw(3) << int(m_atoms.second) << std::setw(3) << int(m_type) << "  0  0  0  0\n";
        }
    };

    struct molecule {
        std::string     m_header;
        std::string     m_comment;
        std::vector<atom> m_atoms;
        std::vector<bond> m_bonds;
//     molecule( const molecule &other )
//         : m_header(other.m_header),
//         m_comment(other.m_comment),
//         m_atoms( other.m_atoms ),
//         m_bonds( other.m_bonds )
//     {}
//
//     molecule &operator=(const molecule &other ) {
//         m_header = other.m_header;
//         m_comment = other.m_comment;
//         m_atoms = other.m_atoms;
//         m_bonds = other.m_bonds;
//
//         return *this;
//     }

        molecule() {}

        molecule( const char *header, const char *comment ) : m_header( header ), m_comment(comment)
        {}

        inline int size() const {

            return int(m_atoms.size());
        }
        inline bool operator<( const molecule &other ) {
            return size() < other.size();
        }
        //std::map<std::string,std::string> m_extra;
        molecule &swap( molecule &other ) {
            m_header.swap( other.m_header );
            m_comment.swap( other.m_comment );
            m_atoms.swap( other.m_atoms );
            m_bonds.swap( other.m_bonds );

            return *this;
        }
    };


};


struct sdf_int_eco {
    typedef sdf_int_eco sdf_int;
    struct atom {
        typedef  boost::array<char,3> ele_t;
        ele_t m_ele;
        uint8_t     m_can_ele;
        //std::string     m_extra;

        atom( float x, float y, float z, char *ele, int can_ele ) : m_can_ele(can_ele)
        {
            //assert(can_ele >= 0 && can_ele <= 255 );

            if ( can_ele < 0 || can_ele > 255 ) {

                throw std::runtime_error( "can_ele out of range for sdf_eco" );
            }

            m_ele = to_ele( ele );

            //      printf( "atom: %f %f %f %c\n", x,y,z,ele);
        }

        atom( float x, float y, float z, ele_t ele, int can_ele ) : m_ele(ele), m_can_ele(can_ele)
        {
            //assert(can_ele >= 0 && can_ele <= 255 );

            if ( can_ele < 0 || can_ele > 255 ) {

                throw std::runtime_error( "can_ele out of range for sdf_eco" );
            }


            //      printf( "atom: %f %f %f %c\n", x,y,z,ele);
        }

        static inline ele_t to_ele( char *ele ) {

            ele_t r;

            size_t nns = 0;
            for ( size_t i = 0; i < ele_t::size(); i++ ) {

                if ( !isspace(ele[i] ) ) {
                    nns++;

                }
            }

            if ( nns < 1 || nns > 2 ) {
                std::stringstream ss;
                ss << "nns < 1 || nns > 2: '" << ele[0] << ";" << ele[1] << ";" << ele[2] << "'";


                throw std::runtime_error( ss.str() );
            }



            std::copy( ele, ele + ele_t::size(), r.begin() );
            return r;
        }
        void print_mdl( std::ostream &os ) const {
            os << "    0.0000    0.0000    0.0000" <<  m_ele[0] << m_ele[1] << m_ele[2] << "  0  0  0  0  0  0  0  0  0  0  0  0\n";
        }
    };
    struct bond {
        std::pair <uint16_t,uint16_t>     m_atoms;
        uint8_t                     m_type;
        uint8_t                     m_can_type;
        // std::string             m_extra;


        bond( unsigned int first, unsigned int second, unsigned int type, unsigned int can_type ) : m_atoms(first,second), m_type(type), m_can_type(can_type)
        {
            if ( first > 999 || second > 999 || type > 255 || can_type > 255 ) {

                throw std::runtime_error( "bond parameter out of range for sdf_eco" );
            }
            //    printf( "bond: %d %d %d\n", m_atoms.first, m_atoms.second, m_type );
        }
        void print_mdl( std::ostream &os ) const {
            os << std::setw(3) << int(m_atoms.first) << std::setw(3) << int(m_atoms.second) << std::setw(3) << int(m_type) << "  0  0  0  0\n";
        }
    };

    struct molecule {
        std::string     m_header;
        std::string     m_comment; // keep the comment field, but do not fill it
        std::vector<atom> m_atoms;
        std::vector<bond> m_bonds;

//     molecule( const molecule &other )
//         : m_header(other.m_header),
//         m_atoms( other.m_atoms ),
//         m_bonds( other.m_bonds )
//     {}
//
//     molecule &operator=(const molecule &other ) {
//         m_header = other.m_header;
//         m_atoms = other.m_atoms;
//         m_bonds = other.m_bonds;
//
//         return *this;
//     }


        molecule() {}

        molecule( const char *header, const char *comment ) : m_header(header)
        {}

        inline int size() const {

            return int(m_atoms.size());
        }
        inline bool operator<( const molecule &other ) {
            return size() < other.size();
        }
        //std::map<std::string,std::string> m_extra;

        molecule &swap( molecule &other ) {
            m_header.swap( other.m_header );
            m_comment.swap( other.m_comment );
            m_atoms.swap( other.m_atoms );
            m_bonds.swap( other.m_bonds );

            return *this;
        }

    };

};


template<class sdf_int_>
class sdf_impl {
    class sanity_check : public std::runtime_error {
    };

public:
    typedef typename sdf_int_::sdf_int sdf_int;
    typedef typename sdf_int_::atom::ele_t atom_ele_t;
    typedef typename sdf_int_::atom atom;
    typedef typename sdf_int_::bond bond;
    typedef typename sdf_int_::molecule molecule;
private:


    std::map<int,int> m_bondtype_map;
    std::map<atom_ele_t,int> m_atomtype_map;

    const bool m_allow_hydrogen;

    int canonicalize_element( const atom_ele_t &ele ) {


        typename std::map< atom_ele_t, int >::iterator it = m_atomtype_map.find( ele );
        if ( it == m_atomtype_map.end() ) {
            int nid = int(m_atomtype_map.size());
            m_atomtype_map[ele] = nid;
            return nid;
        } else {

            return it->second;
        }

    }

    int canonicalize_element(char *ele) {

        atom_ele_t et = atom::to_ele(ele);
        return canonicalize_element(et);

    }

    int canonicalize_bond_type(int type) {
        std::map< int, int >::iterator it = m_bondtype_map.find(type);
        if ( it == m_bondtype_map.end() ) {
            int nid = int(m_bondtype_map.size());
            m_bondtype_map[type] = nid;
            return nid;
        } else {

            return it->second;

        }
    }

    std::vector<molecule> m_molecules;
    bool parse_molecule( std::vector<molecule> &cont, std::istream &is) {
        // parse one molecule block from the stream.
        // implented as a 'free interpretation' of (Arthur Dalby et al., J. Chem. Inf. Comput. Sci, 1992, 32, 244-255) pdf: http://pubs.acs.org/doi/pdf/10.1021/ci00007a012


        const size_t line_len = 256;
        char line[line_len];
        line[0] = 0;

        while ( line[0] == 0 && !is.eof() ) {
            is.getline(line, line_len);
        }
        if ( is.eof()) {
            return false;
        }


        cont.resize( cont.size() + 1 );
        molecule &mol = cont.back();

        mol.m_header = line;
        mol.m_header.erase( std::remove_if(mol.m_header.begin(), mol.m_header.end(), isspace ), mol.m_header.end() );

        is.getline(line, line_len);
        assert( !is.eof() );

        //mol.m_comment = line;
        is.getline(line, line_len);
        assert( !is.eof() );
        int natoms, nbonds;
        // parse 'count line'

        {
            char nt[5];
            is.getline(line, line_len);
            assert( !is.eof() );

            nt[3] = 0;

            std::copy( line, line + 3, nt );
            natoms = atoi(nt);

            std::copy( line + 3, line + 3 + 3, nt );
            nbonds = atoi(nt);

//             printf( "%d %d\n", natoms, nbonds );

            //assert( memcmp( line + 33, " V2000", 6 ) == 0 );
        }

        mol.m_atoms.reserve( natoms );
        mol.m_bonds.reserve( nbonds );

        // parse 'atom block'
        for ( int i = 0; i < natoms; i++ ) {
            is.getline(line, line_len);
            assert( !is.eof() );

            char nt[11];
            nt[10] = 0;

            std::copy( line, line + 10, nt );
            const float x = float(atof( nt ));

            std::copy( line + 10, line + 10 + 10, nt );
            const float y = float(atof( nt ));

            std::copy( line + 20, line + 20 + 10, nt );
            const float z = float(atof( nt ));

            std::copy( line + 30, line + 30 + 3, nt );
            nt[3] = 0;

            // find first non-space character, which is the element symbol
//             char *ele = nt;
//             while( *ele != 0 && isspace(*ele)) {
//                 ele++;
//             }
//
//             char *ele_end = ele;
//
//             while( *ele_end != 0 && !isspace(*ele_end)) {
//                 ele_end++;
//             }
//             assert( ele_end - ele == 1 || ele_end - ele == 2 );

//             char *ele = line + 30 + 3;

            //  printf( "ele: '%c' '%c' '%c'\n", nt[0], nt[1], nt[2] );

            if ( !m_allow_hydrogen ) {

                int nns = 0;
                bool is_h = false;

                for ( int i = 0; i < 3; i++ ) {
                    if ( !isspace(nt[i]) ) {
                        nns++;
                    }
                    if ( nt[i] == 'h' || nt[i] == 'H' ) {
                        is_h = true;
                    }
                }
                if ( nns == 1 && is_h ) {
                    throw std::runtime_error( "molecule in sdf contains hydrogen, while allow_hydrogen is not set" );
                }

            }


            mol.m_atoms.push_back(atom(x,y,z,nt,canonicalize_element(nt)));
        }

        // parse 'bond block'
        for ( int i = 0; i < nbonds; i++ ) {
            is.getline(line, line_len);
            assert( !is.eof() );
            char nt[4];
            nt[3] = 0;
            std::copy( line, line + 3, nt );
            int first = atoi(nt);

            std::copy( line + 3, line + 3 + 3, nt );
            int second = atoi(nt);

            std::copy( line + 6, line + 6 + 3, nt );
            int type = atoi(nt);

            mol.m_bonds.push_back(bond(first,second,type, canonicalize_bond_type(type)));
        }

        // just ignore the rest for now:
        // consume input up to and including the next occurence of '$$$$' + newline
        {
            // there seems to be some weird stuff in the 'associated data' section of some sdf files, which confuses the ifstream::getline ...
            char c;
            while (!is.eof()) {
                // (1) find the next '$'
                do {
                    is.get(c);
                } while ( c != '$' && !is.eof() );


                // (2) match '$$$' or go back to (1)
                is.get(c);
                if ( c != '$' ) {
                    continue;
                }
                is.get(c);
                if ( c != '$' ) {
                    continue;
                }
                is.get(c);
                if ( c != '$' ) {
                    continue;
                }
                // it should be save to use getline to skip the newline now (TODO: is ifstream::getline actually windows line end safe?)
                is.getline( line, line_len );
//                 printf( "line: '%s'\n", line );

                break;
            }
        }
//         while( memcmp( line, "$$$$", 4 ) != 0 && !is.eof() ) {
//             is.getline(line, line_len);
//         }
        //printf( "mol: %zd %zd\n", mol.m_atoms.size(), mol.m_bonds.size() );
        return true;

    }

public:
    sdf_impl( std::istream &is, bool allow_hydrogen = true ) : m_allow_hydrogen(allow_hydrogen) {
        append(is);
//         getchar();
    }


    void append( std::istream &is ) {
        while ( parse_molecule( m_molecules, is ) ) {
            // printf( "mol: %zd\n", m_molecules.size() );
        }
    }

    const std::vector<molecule> &get_molecules() {

        return m_molecules;
    }

    std::vector<const molecule *> get_molecule_ptrs( size_t minsize = 0 ) {
        std::vector<const molecule *> molptr;//( m_molecules.size() );

        for ( typename std::vector< molecule >::const_iterator it = m_molecules.begin(); it != m_molecules.end(); ++it ) {
//             printf( "put: %p\n", &(*it));
            //if( it->m_header != "20637" && it->m_header != "20715" ) {

            if ( it->m_atoms.size() >= minsize ) {
                molptr.push_back(&(*it));
            }
            //}
        }

        return molptr;
    }

    int num_can_atom_types() {
        return int(m_atomtype_map.size());
    }


    static void pad( char c, const std::string &s, int w, std::ostream &os ) {

        std::stringstream sp;
        size_t len = 0;
        while ( len + s.size() < size_t(w) ) {
            sp << c;
            len++;
        }

        os << sp.str() << s;
    }
    static void pad( char c, int v, int w, std::ostream &os ) {
        std::stringstream ss;
        ss << v;
        pad( c, ss.str(), w, os );
    }

    static void pad( char c, double v, int w, std::ostream &os ) {
        std::stringstream ss;
        ss << v;
        pad( c, ss.str(), w, os );
    }


    static void write_mdl( const molecule &mol, std::ostream &os ) {
        os << mol.m_header << "\n";
        //os << "\n\n";
        os << " " << mol.m_comment << "\n\n";
        os << std::setw(3) << mol.m_atoms.size() << std::setw(3) << mol.m_bonds.size() << "  0  0  0  0  0  0  0  0  1 V2000\n";


        for ( typename std::vector<atom>::const_iterator it = mol.m_atoms.begin(); it != mol.m_atoms.end(); ++it ) {

            //os << it->m_ele[1] << it->m_ele[2] << "\n";
            it->print_mdl( os );
        }

        for ( typename std::vector<bond>::const_iterator it = mol.m_bonds.begin(); it != mol.m_bonds.end(); ++it ) {

            //os << "bond: " << int(it->m_atoms.first) << " " << int(it->m_atoms.second) << " " << int(it->m_type) << "\n";
            it->print_mdl( os );
        }
        os << "M  END\n$$$$\n";

    }


    class mol_builder {
    public:
        typedef typename sdf_int_::sdf_int sdf_int;
        typedef typename sdf_int_::atom::ele_t atom_ele_t;
        typedef typename sdf_int_::atom atom;
        typedef typename sdf_int_::bond bond;
        typedef typename sdf_int_::molecule molecule;
    private:

        molecule m_mol;

    public:
        mol_builder() :m_mol( "default", "built by sdf_mol_builder_impl") {

        }
        void add_atom( unsigned int aid, atom_ele_t &at ) {

            if ( aid > m_mol.m_atoms.size() ) {

                throw std::runtime_error( "aid > m_mol.m_atoms.size()" );
            }

            if ( aid == m_mol.m_atoms.size() ) {

                m_mol.m_atoms.push_back( atom( 0.0, 0.0, 0.0, at, 0 ) );
                //             m_mol.add_atom( atom( 0.0, 0.0, 0.0, at, 0 ) );
            } else {

                if ( m_mol.m_atoms[aid].m_ele != at ) {
                    throw std::runtime_error( "inconsistence: m_mol.m_atoms[aid].m_ele != at" );
                }
            }
        }

        void add_bond( int aid1, int aid2, atom_ele_t at1, atom_ele_t at2, int bt ) {


            if ( aid1 < aid2 ) {
                add_atom( aid1, at1 );
                add_atom( aid2, at2 );
            } else {
				// FIXME: is this supposed to be the same as the 'then' branch?
                add_atom( aid1, at1 );
                add_atom( aid2, at2 );
            }

            m_mol.m_bonds.push_back( bond( aid1 + 1, aid2 + 1, bt, 0 ));
        }

        const molecule &get_mol() {

            return m_mol;
        }

        static void print_mol( const molecule &mol ) {
//             std::cout << "molecule: " << mol.m_atoms.size() << " " << mol.m_bonds.size() << "\n";


        }


    };
};


// template<typename sdf_int_>
// inline bool operator<(const typename sdf_int_::molecule &m1, const typename sdf_int_::molecule &m2 ) {
//     return m1.size() < m2.size();
// }

#ifndef WIN32
extern template class sdf_impl<sdf_int_full>;
extern template class sdf_impl<sdf_int_eco>;
#endif
typedef sdf_impl<sdf_int_full> sdf_full;
typedef sdf_impl<sdf_int_eco> sdf_eco;

} // namespace ivy_mike
#endif

