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

#ifndef __ivy_mike__fasta_h
#define __ivy_mike__fasta_h

#include <stdint.h>

#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <stdexcept>
#include <cassert>

#include <cstdlib>
#include <cstddef>
#include <stdexcept>
#include <sstream>

#include "ivymike/aligned_buffer.h"

namespace ivy_mike {

#if !defined( WIN32 ) && !defined(__native_client__)
#include <sys/types.h>
#include <sys/stat.h>

#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

class mapped_file {
    int m_fd;
    
    void *m_base;
    size_t m_size;
    size_t m_ptr;
    
    mapped_file() {}
    mapped_file( const mapped_file &mf );
    const mapped_file &operator=(const mapped_file &other );
public:
    
    mapped_file( const char *name ) 
        : m_ptr(0)
    {
        m_fd = open( name ,O_RDONLY );
        
        if( m_fd == -1 ) {
         
            throw new std::runtime_error( "cannot open file for mapping" );
        }
        
        m_size = lseek( m_fd, 0, SEEK_END );
        
        
        m_base = mmap( 0, m_size, PROT_READ, MAP_SHARED, m_fd, 0 );
        
        if( m_base == 0 ) {
            throw std::runtime_error( "mmap failed" );   
        }
        
        madvise( m_base, m_size, MADV_SEQUENTIAL );
        
    }
    
    ~mapped_file() {
        munmap( m_base, m_size );
    }


    inline void check_bounds() {
     
        
    }
    int get() {
        if( eof() ) {
            return -1;
        } else {
            return ((char*)m_base)[m_ptr++];
        }
    }
 
    void unget() {
        if( m_ptr > 0 ) {
            m_ptr--;
        }
    }
 
    bool eof() {
        return m_ptr >= m_size;   
    }
    bool good() {
        return !eof();
    }
    
    void seekg( off_t ptr ) {
        m_ptr = ptr;   
    }
    void clear() {
        
    }
};
#endif

template<class input>
static inline void get_line( input &is, std::vector<char> &linebuf ) {
    linebuf.resize(0);
    while(true) {
        char c = is.get();
        if( c == '\n' || is.eof() ) {
            break;
        }
        linebuf.push_back(c);
    }

}

static inline bool xisnewline( int c ) {
    return c == '\n' || c == '\r';
}

static inline bool xisspace( int c ) {
    return c == ' ' || xisnewline(c); //c == '\n';   
}



// template<class input>
// static void read_fasta( input &is, std::vector<std::string> &names, std::vector<std::vector<uint8_t> > &data ) {
//  
//     // WARNING: this whole fasta reader is utter crap. I am very ashamed of it...
//     
// //     std::vector<char> linebuf;//1024 * 1024); // uhm, this is ugly...
//     //std::string linebuf;
//     
//     names.clear();
//     data.clear();
//     std::vector<uint8_t> *data_accum = 0;
//     
//     while( is.good() && !is.eof() ) {
//     
//         //is.getline( linebuf.data(), linebuf.size() );
//             
//             
//         int c = is.get();
// //         get_line( is, linebuf );
// //         std::vector< char >::const_iterator li = linebuf.begin();
//         
//         if( is.eof() ) {
//             break;
//         }
//         
//         if( c == '>' ) {
//             
//             
//             while( isspace(is.get()) ) {}
//             is.unget();
//             
//             names.push_back(std::string());
//             std::string &str = names.back();
//         
//             while( true ) {
//                 c = is.get();
//                 
//                 if( xisspace(c) || is.eof() ) {
//                     break;   
//                 }
//                 str.push_back(c); 
//             }
//             
//         
// //             std::cout << "name: " << names.back() << std::endl;
//             
//             data.push_back( std::vector<uint8_t>() );
//             //data_accum = &data.back();
//             data_accum = &data.back();
//         } else {
//             is.unget();
//             if( data_accum == 0 ) {
//                 throw std::runtime_error( "data_accum == 0. probably bad fasta file\n" );
//             }
//             
//             while( xisspace( is.get() )) {}
//             
// //             std::cout << "xxx: " << int(xxx) << "\n";
//             
//             if( !is.eof() ) {
//                 
//                 is.unget();
//                 
//                 while( true ) {
//                     c = is.get();
// //                     std::cout << "c: " << int(c) << "\n";
//                     if( c == -1 ) {
//                         std::cout << " c:" << int(c) << " " << is.eof() ;
//                     }
//                     if( xisspace(c) || is.eof() ) {
//                         break;   
//                     }
//                     data_accum->push_back(c); 
//                     
//                 }
//             }
//        }
//     }
//     
// //     std::cout << "size: " << names.size() << " " << data.size() << "\n";
//     data.resize( names.size() );
//     
// }

template<class input, class StateMap>
class inc_fasta {
    input &m_input;
    
    const StateMap &m_state_map;
    
public:
    inc_fasta( input &inp, const StateMap &state_map, bool do_reset = true ) : m_input(inp), m_state_map(state_map) {
        
        if( do_reset ) {
            reset();  // why the hell is this in here? completely blows up when reading from a pipe...
        }
     
    }
    
    void reset() {
        m_input.clear();
        m_input.seekg(0);   
    }
    
    bool next_seq( std::string &name, std::vector<uint8_t> &seq ) {
        
      
        
        
        // look for '>'
        while( m_input.good() && m_input.get() != '>' ) {}

        if( !m_input.good() ) {
//             std::cerr << "not good\n";
            return false;
        } else {
//             std::cerr << "good\n";
        }

        int c;
        
        // skip spaces
        while( xisspace(m_input.get() )) {}
        
        m_input.unget();
        
        // read name until first space
        while( true ) {
            c = m_input.get();
            
            if( xisspace(c) || !m_input.good() ) {
                break;   
            }
            name.push_back(c); 
        }
        
        // read remaining name line until first newline
        while( !m_input.eof() && !xisnewline(c) ) {
            c = m_input.get();
        }
        
        
        // gobble up remaining spaces (the newline check above could have been triggered by @#$%* windows newline, aaahhhhrggg!)    
        while( xisspace( m_input.get() )) {}
        
        m_input.unget();
            
        while( true ) {
            c = m_input.get();
            
            if( !m_input.good() ) {
                break;
            }
            
            if( c == '>' ) {
                m_input.unget();
                break;
            }
                
//             std::cout << "c: " << int(c) << "\n";
            if( !isspace(c) && m_state_map.state_valid(c) ) {
                seq.push_back( m_state_map.state_backmap(c)); 
            }
                
        }
        return true;
    }   
    
};

namespace waiting_for_N1427 {
struct null_backmap {
    
    inline uint8_t state_backmap( uint8_t c ) const {
        return c;
    }
    
    inline bool state_valid( int ) const {
        return true; // everything goes...
    }
};
}

template<class input>
static void read_fasta( input &is, std::vector<std::string> &names, std::vector<std::vector<uint8_t> > &data, bool do_reset = true ) {
    
    // optionally read_fasta can take a state-map object (i.e., normally a scoring_matrix object).
    // use 'neutral' mapping if none explicitely given
    waiting_for_N1427::null_backmap nb;
    
    inc_fasta<input,waiting_for_N1427::null_backmap> f(is, nb, do_reset );
    
//     std::cerr << "goodx: " << is.good() << std::endl;
    
    while( true ) {
        names.push_back(std::string());
        data.push_back(std::vector<uint8_t>());
        
        bool success = f.next_seq( names.back(), data.back() );
        


        if( !success ) {
            names.pop_back();
            data.pop_back();
            break;
        }

        if( data.back().empty() ) {
			std::cout << "empty: " << names.back() << "\n";
			names.pop_back();
			data.pop_back();
		}

    }
}

template<typename input, typename StateMap>
static void read_fasta( input &is, const StateMap &sm, std::vector<std::string> &names, std::vector<std::vector<uint8_t> > &data ) {
    inc_fasta<input,StateMap> f(is, sm );
    
    while( true ) {
        names.push_back(std::string());
        data.push_back(std::vector<uint8_t>());
        
        bool success = f.next_seq( names.back(), data.back() );
        
        if( !success ) {
            names.pop_back();
            data.pop_back();
            break;
        }
    }
}



// static void write_fasta( std::ostream &os, std::vector<std::string> &names, std::vector<std::string> &data ) {
//  
//     assert( names.size() == data.size() );
//     
//     for( size_t i = 0; i < names.size(); i++ ) {
//      
//         os << ">" << names[i] << "\n";
//         os << data[i] << "\n";
//     }
//         
//         
//         
// }

class scoring_matrix {
public:
    typedef int8_t score_t;
private:
    
    const static size_t MAX_SIZE = 256;
    
    aligned_buffer<char> m_cmatrix;
    size_t m_cmatrixsize;
    
    //std::vector<score_t> m_cmatrix;
    std::vector<score_t> m_matrix;
    std::vector<int> m_backmap;
    std::string m_alphabet;
    size_t addr( size_t a, size_t b ) const {
        return a + b * MAX_SIZE;
    }
    
    size_t caddr( size_t a, size_t b ) {
        return a + b * m_cmatrixsize;
    }
    
    void compress() {
        
        
        
        // the cmatrix is always at least one column/row larger than the number of states. The last row/column is for a 'zero-state'
        // which returns a score of 0 against all states, to make qprofile padding easier.
        const size_t asize = m_alphabet.size();
        
        // the pysical dimension of the cmatrix are different from the alphabet size, for better cache-line alignment.
        // 32 should give good alignment on current hardware (two cachelines per row?)
        //m_cmatrixsize = asize+1;
        m_cmatrixsize = 32;
        if( m_cmatrixsize < asize + 1 ) {
            // the current standard of 32 bytes is large enough for protein sequences
            throw std::runtime_error( "m_cmatrixsize < asize. alphabet too large." );
            
        }
        
        
        //const size_t asize = m_alphabet.size() + 1;
        m_cmatrix.resize( m_cmatrixsize * m_cmatrixsize );
        std::fill( m_cmatrix.begin(), m_cmatrix.end(), 0 );
//         std::cerr << "cmatrix: " << m_cmatrix.size() << " " << m_cmatrixsize << "\n";
        
        // copy scores from the raw scoring matrix to the compressed matrix
        for( size_t i = 0; i < asize; i++ ) {
            for( size_t j = 0; j < asize; j++ ) {
                
                m_cmatrix[caddr(i,j)] = m_matrix[addr(m_alphabet[i], m_alphabet[j])];
                
            }
        }
    }
public:
    template<typename cont>
    class valid_state_appender {
        cont &m_cont;
        const scoring_matrix &m_sm;
    public:
        valid_state_appender<cont>( const scoring_matrix &sm, cont &container ) : m_cont(container), m_sm(sm) {
            
        }
        
        inline void operator()( int c ) {
            if( c >= 0 && c < int(scoring_matrix::MAX_SIZE) ) {
                
                c = toupper(c);
                
                if( m_sm.m_backmap[c] != -1 ) {
                    m_cont.push_back( m_sm.m_backmap[c] );
//                     std::cout << int(m_cont.back()) << "\n";
                }
            }
        }
        
       
        
    };
    
    inline size_t num_states() const {
        return m_alphabet.size();
        
    }
    
    inline int get_state( int i ) {
        return m_alphabet[i];
    }
    
    inline bool state_valid( int s ) const {
        if( s > 0 && size_t(s) < MAX_SIZE ) {
            return m_backmap[s] != -1;
        } else {
            return false;
        }
    }

    inline int state_backmap( int s ) const {
        if( s < 0 || size_t(s) >= MAX_SIZE ) {
            throw std::runtime_error( "sequence character out of bounds" );
        }
        
        if( m_backmap[s] == -1 ) {
            std::stringstream ss;
            ss << "invalid sequence character " << int(s);
            
            if( std::isprint(s) ) {
                ss << "(" << char(s) << ")"; 
            }
            
            throw std::runtime_error( ss.str() );
        }
        
        return m_backmap[s];
    }

    inline score_t get_score( int a, int b ) const {
//         assert( a >= 0 && a < MAX_SIZE && b >= 0 && b < MAX_SIZE );
        return m_matrix[addr(a,b)];
    }
    
    inline score_t *get_slice( size_t a ) {
        return &m_matrix[a * MAX_SIZE];
    }
    
    inline int get_zero_state() const {
        return int(m_alphabet.size());
    }
    
    inline const char *get_cslice( size_t a ) const {
        //return &m_cmatrix[a * m_alphabet.size()];
        return m_cmatrix(a * m_cmatrixsize);
    }
    
    scoring_matrix( int match, int mismatch )
        : m_matrix( MAX_SIZE * MAX_SIZE ),
            m_backmap(MAX_SIZE)
    {
                
        std::fill(m_backmap.begin(), m_backmap.end(), -1);
        
        const char *acgt = "ACGT";
        
        for( int i = 0; i < 4; i++ ) {
            const char ci = acgt[i];
            m_backmap[ci] = int(m_alphabet.size());
            m_alphabet.push_back(ci);
            
            
            for( int j = 0; j < 4; j++ ) {
                const char cj = acgt[j];
                
                int score;
                if( i == j ) {
                    score = match;
                } else {
                    score = mismatch;
                }
                
                m_matrix[addr( ci, cj )] = score;
                m_matrix[addr( cj, ci )] = score;
                
            }
            
        }
        
        compress();
    
    }
    
    
    scoring_matrix( std::istream &is ) 
        : m_matrix( MAX_SIZE * MAX_SIZE ),
            m_backmap(MAX_SIZE)
        
    {
        std::fill(m_backmap.begin(), m_backmap.end(), -1);
        std::vector<char> linebuf;
        
        bool have_firstline = false;
//         size_t mline = 0;
        while( is.good() && !is.eof() ) {
            
            char c;
            
           
            do {
                c = is.get();
            } while( std::isspace(c) );
        
            if( is.eof() ) {
                break;
            }
            
//            get_line( is, linebuf );
                
            if( c == '#' ) {
                while( !is.eof() && is.get() != '\n' ) {}
                continue;
            }
            
  
            
            if( !have_firstline ) {
                while( c != '\n' ) {
                    c = toupper(c);
                    
                    m_backmap[c] = int(m_alphabet.size());
                    m_alphabet.push_back(c);
                    
                    do {
                        c = is.get();
                        
                        if( c == '\n' ) {
                            break;
                        }
                    } while( std::isspace(c) );
                }
            
//                 std::cout << "alphabet: " << m_alphabet << "\n";
                have_firstline = true;
            } else {
                is.unget();
//                 char c = is.get();
//                 std::cout << "c: " << c << "\n";
                
                std::string lname;
                is >> lname;
                
//                 std::cout << "name: " << lname << "\n";
                
                char lnc = toupper(lname[0]);
                
                for( size_t i = 0; i < m_alphabet.size(); i++ ) {
                 
                    int score;
                    
                    is >> score;
                    
//                     std::cout << "score: " << int(lnc) << " " << int(m_alphabet[i]) << " " << score << "\n";
                    m_matrix[addr( lnc, m_alphabet[i] )] = score;
                    m_matrix[addr( m_alphabet[i], lnc )] = score;
                }
            }
        }
        
        compress();
//         for( size_t i = 0; i < m_alphabet.size(); i++ ) {
//             for( size_t j = 0; j < m_alphabet.size(); j++ ) {
//                 std::cout << int(get_score( m_alphabet[i], m_alphabet[j] )) << " ";
//             }
//             std::cout << "\n";
//             
//             
//         }
     
    }
    
    inline int operator()( int c ) const {
        return state_backmap(c);
    }
};

}
#endif
