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


#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <typeinfo>
#include <cctype>
#include <algorithm>
#include "ivymike/getopt.h"

// OK, the parser is a bit overengineered for cmdline parsing. it's turned out to be some kind of backtracking recursive descent parser...
// maybe I can use the bits for other parsers in the future...

namespace ivy_mike {
namespace getopt {
    
// template <typename pinput_t, int X = 0>
// class pisaver {
//     pinput_t m_pis;
//     pinput_t &m_pi_ref;
//     bool     m_commit;
// 
// public:
//     pisaver( pinput_t &pi ) : m_pis(pi), m_pi_ref(pi), m_commit(false) {}
//     ~pisaver() {
//         if ( !m_commit ) {
//             if( false ) {
//                 std::cout << "pis rollback " << X << " at " << (m_pi_ref.end() - m_pi_ref.begin());
//                 
//                 if ( m_pi_ref.has_next() ) {
//                     std::cout << " at char: " << m_pi_ref.peek() << "\n";
//                 } else {
//                     std::cout << " at end\n";
//                 }
//             }
//             
//             // do the rollback
//             m_pi_ref = m_pis;
// 
//         }
//     }
// 
//     void commit() {
//         m_commit = true;
//     }
// };
//     
// parser::auto_token parser::parse_argument(ivy_mike::getopt::parser::my_pinput& pi) {
//     pisaver<my_pinput> pis(pi);
// 
// 
//     pi.skip_whitespace();
// 
//     if ( pi.peek() == '-' ) {
//         pi.next();
// 
//         unsigned char optid = pi.peek();
// 
//         if ( std::isalnum(optid) && m_options[optid] ) {
//             pi.next();
//             pis.commit();
//             return auto_token(new option( optid ));
//         }
// 
//     }
// 
// 
//     return std::auto_ptr<token>();
// 
// }
// 
// parser::auto_token parser::parse_bare_string(parser::my_pinput& pi) {
//     pisaver<my_pinput,1> pis(pi);
// 
//     pi.skip_whitespace();
// 
// 
//     if ( !pi.has_next()) {
// 
//         return auto_token();
//     }
// 
// 
// 
//     std::string str;
// 
//     while ( pi.has_next() ) {
// 
//         if ( std::isprint( pi.peek() ) && !std::isspace(pi.peek()) && pi.peek() != '"' ) {
//             str.push_back(pi.peek());
// 
//             pi.next();
//         } else {
//             pi.next();
//             break;
//         }
// 
//     }
// 
//     if ( str.empty() ) {
//         return auto_token();
//     }
// 
// 
//     auto_token s(new string( str ));
//     pis.commit();
// 
//     return s;
// 
// }
// parser::auto_token parser::parse_quoted_string(parser::my_pinput& pi) {
//     pisaver<my_pinput,2> pis(pi);
// 
//     pi.skip_whitespace();
// 
// 
//     if ( !pi.has_next()) {
//         return auto_token();
//     }
// 
// 
//     if ( pi.peek() != '"' ) {
//         return auto_token();
//     }
//     pi.next();
// 
//     std::string str;
//     bool closed = false;
//     while ( pi.has_next() ) {
// 
//         if ( pi.peek() == '"' ) {
//             closed = true;
//             pi.next();
//             break;
//         } else {
//             str.push_back(pi.peek());
//             pi.next();
//         }
// 
//     }
// 
//     if ( str.empty() || !closed) {
//         return auto_token();
//     }
// 
// 
//     auto_token s(new string( str ));
//     pis.commit();
// 
//     return s;
// 
// }
// parser::auto_token parser::parse_string(parser::my_pinput& pi) {
//     // the two string parsers do the pi state reconstruction
//     auto_token t(parse_bare_string(pi));
// 
//     if ( t.get() == 0 ) {
//         t = parse_quoted_string(pi);
//     }
// 
//     return t;
// 
// }
// parser::auto_token parser::parse(parser::my_pinput& pi) {
//     if ( !pi.has_next() ) {
//         return auto_token();
//     }
// 
//     auto_token t = parse_argument( pi );
// 
//     if ( t.get() != 0 ) {
//         return t;
//     }
// 
//     t = parse_string( pi );
// 
//     if ( t.get() != 0 ) {
//         return t;
//     }
// 
//     return auto_token();
// 
// }
// void parser::parse_main(parser::my_pinput& pi) {
//     
//     while ( true ) {
//         auto_token t = parse(pi);
//         
//         if( t.get() == 0 ) {
//             break;
//         }
//         
// //         std::cout << "add token: " << t->to_string() << "\n";
//         m_tokenstream.push_back( t.get() );
//         t.release();
//     }
// 
//     if ( pi.has_next() ) {
//         std::cerr << "WARNING: parser not finished\n";
//     }
// 
// }




template <typename pinput_t, int X = 0>
class pisaver {
    pinput_t m_pis;
    pinput_t &m_pi_ref;
    bool     m_commit;

public:
    pisaver( pinput_t &pi ) : m_pis(pi), m_pi_ref(pi), m_commit(false) {}
    ~pisaver() {
        if ( !m_commit ) {
            if( false ) {
//                 std::cout << "pis rollback " << X << " at " << (m_pi_ref.end() - m_pi_ref.begin());
//                 
//                 if ( m_pi_ref.has_next() ) {
//                     std::cout << " at char: " << m_pi_ref.peek() << "\n";
//                 } else {
//                     std::cout << " at end\n";
//                 }
            }
            
            // do the rollback
            m_pi_ref = m_pis;

        }
    }

    void commit() {
        m_commit = true;
    }
};


void parser::parse_option( parser_input &pi ) {
    if( pi.past_end() ) {
        return;
    }
    
     pisaver<parser_input,0> pis(pi);
     
     char *arg = pi.peek();
//      std::cout << "arg: " << arg << "\n";
     if( arg[0] != '-' ) {
        throw std::runtime_error( "WARNING: expected option. bailing out." );
     }
     
     char opt = arg[1];
     if( !std::isalnum(opt) ) {
         throw std::runtime_error( "WARNING: bad option. bailing out." );
     }
     if( !m_options[opt] ) {
         throw std::runtime_error( "WARNING: unknown option. bailing out." );
     }
     
     m_option_count[opt]++;
     pi.next();
     
     if( m_opt_has_argument[opt] ) {
         
         
         if( m_opt_values[opt] != 0 ) {
             if( m_opt_values[opt]->allow_empty() ) {
                 m_opt_values[opt]->set("1"); // HACK
             } else {
                 if( pi.past_end() ) {
                     throw std::runtime_error( "WARNING: option expects argument. bailing out." );
                 }
                 
                 m_opt_values[opt]->set(pi.peek());
                 pi.next();
             }
         }
         
     }
     pis.commit();

    
}

bool parser::parse_new( int argc, char **argv ) {
    
    
    parser_input pi( argc, argv );
    
    
    try {
        while( !pi.past_end() ){
            parse_option(pi);
        } 
    } catch( std::runtime_error x ) {
        std::cerr << x.what() << "\n";
        std::cerr << "exit getopt\n";
        return false;
    }
    
    return true;
}

bool parser::parse(int argc, char** argv) {
//     std::stringstream ss;
// 
//     for ( size_t i = 1; i < size_t(argc); i++ ) {
//         if ( i > 1 ) {
//             ss << " ";
//         }
// 
// 
//         if ( std::count( argv[i], argv[i] + strlen(argv[i]), ' ' ) != 0 ) {
// //             std::cout << "quoting\n";
//             ss << "\"" << argv[i] << "\"";
//         } else {
//             ss << argv[i];
//         }
// 
//     }
//     std::string s = ss.str();
// 
// //     std::cout << "cmdline: '" << s << "'\n" << std::endl;
//     my_pinput pi( s.begin(), s.end() );
// 
//     parse_main( pi );
    
    return parse_new( argc, argv );
//     for ( std::list< token* >::iterator it = m_tokenstream.begin(); it != m_tokenstream.end(); ++it ) {
//         token *t = *it;
// 
//         if ( typeid(*t) == typeid( option ) ) {
//             option *opt = static_cast<option *>(t);
// 
//             std::cout << "option\n";
// 
//             if ( m_options[opt->get_opt()] ) {
//                 m_option_count[opt->get_opt()]++;
// 
//                 if ( m_opt_has_argument[opt->get_opt()] ) {
//                     ++it;
//                     
//                     std::list< token* >::iterator it_next = it;
//                     it_next++;
//                     if ( it_next != m_tokenstream.end() ) {
//                         
//                         t = *it_next;
//    
//                         if ( typeid( *t ) == typeid(string)) {
//                             m_opt_strings[opt->get_opt()] = static_cast<string*>(t)->get_string();
//                             it = it_next;
//                             std::cout << "missing option argument\n";
//                             
//                         } else {
//                             m_opt_strings[opt->get_opt()] = std::string();
//                         }
//                     }
//                     
//                     
//                     if( m_opt_values[opt->get_opt()] != 0 ) {
//                         base_value &opt_value = *m_opt_values[opt->get_opt()];
//                         std::string &opt_string = m_opt_strings[opt->get_opt()];
//                         std::cout << "opt: '" << opt_string << "' " << opt_value.allow_empty() << "\n";
//                        
//                        
//                         if( opt_string.empty() && !opt_value.allow_empty() ) {
//                             std::cout << "option " << opt->get_opt() << " requires an argument\n";
//                             return false;
//                         }
//                         
//                         opt_value.set(opt_string);
//                     }
//                 }
//             }
// 
//         } else {
//             std::cout << "unexpected token: " << typeid(*t).name() << " " << t->to_string() << "\n";
//             return false;
//         }
//     }
    return true;

}
void parser::add_opt(unsigned char c, bool argument) {
    m_options[c] = true;
    if ( argument ) {
        m_opt_has_argument[c] = true;
    }
}

void parser::add_opt_proxy( unsigned char c, std::unique_ptr<base_value> bv ) {
    if( m_opt_values[c] != 0 ) {
        std::stringstream ss;
        ss << "there is already a value proxy registered for option " << c;
        throw std::runtime_error( ss.str() );
    }
    m_options[c] = true;
    m_opt_has_argument[c] = true;
    m_opt_values[c] = bv.get(); // hmm, is vector::operator[] officially allowed to throw exceptions? using sequencial get+release should make this code exception safe anyway!?
    bv.release();
}

parser::~parser() {
    std::for_each( m_tokenstream.begin(), m_tokenstream.end(), delete_object() );
    std::for_each( m_opt_values.begin(), m_opt_values.end(), delete_object() );
}

} // namespace getopt
} // namespace ivy_mike

