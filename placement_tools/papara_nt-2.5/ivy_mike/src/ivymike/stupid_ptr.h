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


#ifndef __ivy_mike__stupid_ptr_h
#define __ivy_mike__stupid_ptr_h
#include <stdexcept>
#include <cstdlib>
#include <cstdio>

namespace ivy_mike {
    
//
// stupid ptr
//
// a smart-ptr which can be 'remote controlled'. 
// Somehow useful for Dependency Injection (especially when injecting pointers to scoped variables using stupid_ptr_guard).
//
template<typename T>
class stupid_ptr {
    T *m_ptr;
    
    stupid_ptr( const stupid_ptr & );
    const stupid_ptr &operator=( const stupid_ptr & );
    
public:
    typedef T *ptr_t;
    
    stupid_ptr() : m_ptr(0) {}
    
    inline void set( T* ptr ) {
        m_ptr = ptr;
    }
    
    inline void clear() {
        m_ptr = 0;
    }
    
    inline bool is_valid_ptr() {
        return m_ptr != 0;
    }
    
    inline T* operator->() {
        if( m_ptr == 0 ) {
            throw std::runtime_error( "accessing null stupid_ptr. Oh noes!" );
        }
        return m_ptr;
    }
    
    inline T& operator*() {
        
        if( m_ptr == 0 ) {
            throw std::runtime_error( "dereferencing null stupid_ptr. Oh noes!" );
        }
        
        return *m_ptr;
    }
    
};


template<typename T>
class stupid_ptr_guard {
    stupid_ptr<T> &m_victim;
    
    stupid_ptr_guard( const stupid_ptr_guard & );
    const stupid_ptr_guard &operator=( const stupid_ptr_guard & );
public:
    stupid_ptr_guard( stupid_ptr<T> &victim, T* ptr ) : m_victim(victim) {
        m_victim.set(ptr);
    }
    
    ~stupid_ptr_guard() {
        m_victim.clear();
    }
    
};

//
// freeer/fcloser: exactly what it's name suggests. For lazy c++ programmers.
//

class freeer {
    void *ptr;

    freeer( const freeer &other ) {}
    const freeer &operator=(const freeer &other ) { return *this; }

public:
    freeer(void *ptr_) : ptr(ptr_) {}
    ~freeer() {
        free(ptr);
    }
};

class fcloser {
    FILE *h_;
    fcloser( FILE *h ) : h_(h) {}
    ~fcloser() {
        if( h_ != 0 ) {
            fclose( h_ );
        }
    }
};
    
}

#endif
