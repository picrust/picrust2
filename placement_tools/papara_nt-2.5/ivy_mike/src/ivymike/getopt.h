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



#ifndef __ivy_mike__getopt_h
#define __ivy_mike__getopt_h

#include <cstdlib>
#include <cctype>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <string>
#include <vector>
#include <set>
#include <list>
#include <map>
#include <algorithm>
#include <memory>
#include <functional>
#include <cctype>
namespace ivy_mike {

template<typename T>
T literal_cast( const std::string &str );

template<>
inline int literal_cast( const std::string &str ) {
    return atoi(str.c_str());
}

template<>
inline std::string literal_cast( const std::string &str ) {
    return str;
}

template<>
inline bool literal_cast( const std::string &str ) {
    std::string tmp = str;
    
    std::transform( tmp.begin(), tmp.end(), tmp.begin(), std::ptr_fun<int,int>(std::toupper) );
    
    //std::cerr << "literal_cast<bool>: '" << str << "'\n";
    
    return str != "FALSE" && str != "F";
}



template <typename T>
class literal_cast_fun {
public:
    T operator()( const std::string &str ) {
        return literal_cast<T>(str);
//         std::stringstream ss;
//         ss << "literal cast not implemented for type " << typeid(T).name();
//
//         throw std::runtime_error( ss.str() );
    }

};


    
struct delete_object {

    template<typename T>
    void operator()(const T* ptr ) {
        delete( ptr );
    }
};

namespace getopt {



template <typename iter_t>
class pinput {
    typedef typename iter_t::value_type value_type;

    iter_t m_begin;
    iter_t m_end;


    pinput();

public:
    pinput( const iter_t &beg, const iter_t &end_ ) : m_begin(beg), m_end(end_) {}

    void skip_whitespace() {
//         while( m_begin != m_end && std::isspace(*m_begin) ) { // FIXME: do iterators generally allow to dereference iter.end()?
//             ++m_begin;
//         }

        m_begin = std::find_if( m_begin, m_end, std::not1(std::ptr_fun<int,int>(std::isspace)) );
    }


    value_type peek() {
        if ( m_begin >= m_end ) {
            throw std::runtime_error( "m_begin >= m_end" );
        }

        return *m_begin;
    }

    value_type next() {
        if ( m_begin >= m_end ) {
            throw std::runtime_error( "m_begin >= m_end" );
        }

        return *(m_begin++);
    }

    bool has_next() {
        return m_begin < m_end;
    }

    const iter_t &begin() {
        return m_begin;
    }
    const iter_t &end() {
        return m_end;
    }
};




class base_value {

public:
    virtual ~base_value() {}

    virtual void set( const std::string &str ) = 0;
    virtual bool allow_empty() = 0;
};

template<typename T>
class value : public base_value {
//     literal_cast<T> lc;

    T & m_vref;
    const bool m_allow_empty;
public:
//     value() : m_vref(0) {}
    value( T &ref, bool allow_empty_ = false) : m_vref(ref), m_allow_empty(allow_empty_) {}

    virtual ~value() {}

    virtual void set( const std::string &str ) {
//         if( m_vref == 0 ) {
//             throw std::runtime_error( "meeep: uninitialized value proxy\n" );
//         }

        m_vref = literal_cast<T>(str);
//         std::cerr << "set: " << m_vref << "\n";
    }

    const value<T> &set_default( const T &dv ) {
//         if( m_vref == 0 ) {
//             throw std::runtime_error( "meeep: uninitialized value proxy\n" );
//         }
        m_vref = dv;
        return *this;
    }
    
   virtual bool allow_empty() { return m_allow_empty; }
};




class parser {
    class parser_input {
        int m_argc;
        char **m_argv;
        int m_ptr;
    public:
        parser_input( int argc, char **argv ) : m_argc(argc), m_argv(argv), m_ptr(1) {}
        
        
        bool past_end() {
            return m_ptr >= m_argc;
        }
        
        void next() {
            m_ptr++;
        }
        
        char *peek() {
            if( m_ptr >= m_argc ) {
                throw std::runtime_error( "peek: out of bounds." );
            }
            
            return m_argv[m_ptr];
        }
        
    };
    
    typedef pinput<std::string::iterator> my_pinput;

    std::vector<bool> m_options;

    std::vector<bool> m_opt_has_argument;

    std::vector<int> m_option_count;
    std::map <char,std::string> m_opt_strings;

    std::vector<base_value *> m_opt_values;

    parser( const parser & );
    parser &operator=(const parser &);

    class token {
    public:
        virtual ~token() {}

        virtual std::string to_string() = 0;
    };

    typedef std::unique_ptr<token> auto_token;
    
    const static auto_token empty_token;
    
    class option : public token {
        const unsigned char m_option;
    public:
        option( unsigned char opt ) : m_option(opt) {}

        virtual ~option() {}

        virtual std::string to_string() {
            std::stringstream ss;

            ss << "option: " << m_option;

            return ss.str();
        }

        unsigned char get_opt() {
            return m_option;
        }

    };

    class string : public token {
        std::string m_string;

    public:

        string( std::string &s ) : m_string(s) {}

        virtual ~string() {}

        virtual std::string to_string() {
            std::stringstream ss;

            ss << "string: " << m_string;

            return ss.str();
        }

        const std::string &get_string() {
            return m_string;
        }
    };


    std::list <token *> m_tokenstream;
    void add_opt_proxy( unsigned char c, std::unique_ptr<base_value> bv ); // takes ownership of bv

public:
    parser() : m_options(256), m_opt_has_argument(256),  m_option_count(256), m_opt_strings(), m_opt_values(256) {}

    ~parser() ;


    void add_opt( unsigned char c, bool argument );

    template<typename T>
    void add_opt( unsigned char c, const value<T> &vo ) {
        // wow, finally found a scenario that fits inheritance...
        add_opt_proxy(c, std::unique_ptr<base_value>(new value<T>(vo)));
    }


//     auto_token parse_argument( my_pinput &pi ) ;
//     auto_token parse_bare_string( my_pinput &pi ) ;
//     auto_token parse_quoted_string( my_pinput &pi ) ;
//     auto_token parse_string( my_pinput &pi ) ;
//     auto_token parse( my_pinput &pi ) ;
//     void parse_main( my_pinput &pi ) ;


    
    bool parse_new( int argc, char **argv );
    void parse_option( parser_input &pi );
    
    
    bool parse( int argc, char **argv ) ;

    const std::string &get_string( char opt ) {
        return m_opt_strings[opt];
    }

    int get_int( char opt ) {
        return atoi( m_opt_strings[opt].c_str() );
    }

    void get_int_if_present( char opt, int &i ) {
        if ( opt_count(opt) != 0 ) {
            i = get_int(opt);
        }
    }


    int opt_count( char opt ) {
        return m_option_count[opt];
    }

};

} // namespace getopt
} // namespace ivy_mike



#endif
