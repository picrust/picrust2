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

#ifndef __ivy_mike__tdmatrix_h
#define __ivy_mike__tdmatrix_h
#include <cstddef>

namespace ivy_mike {
template<typename T>
class odmatrix {
public:
    typedef T value_type;
    typedef T* iterator;
private:
    T* m_base;
    size_t m_size;

    bool m_own;
public:
    odmatrix( T* base, size_t size_ ) : m_base(base), m_size(size_), m_own(false) {}
    
    T &operator[](ptrdiff_t a) {
        return m_base[a];
    }
    
    const T &operator[](ptrdiff_t a) const {
        return m_base[a];
    }
    
    iterator begin() {
        return m_base;
    }
    iterator end() {
        return m_base + m_size;
    }
    
    size_t size() const {
        return m_size;
    }
    
};

template<typename T>
class tdmatrix {
public:
    typedef odmatrix<T> row_type;
    typedef T value_type;
    
    typedef T* iterator;
private:
    
    // WARNING: remember to change swap impl when member vars are changed!

    T* m_base;
    size_t m_asize;
    size_t m_bsize;
    
    bool m_own;

    tdmatrix( const tdmatrix &other ) {}
    const tdmatrix &operator=(const tdmatrix &other ) {return *this;};

public:
    class row_iterator {
        T *m_base;
        size_t m_stride;
        
    public:
        row_iterator( T* base, size_t stride ) : m_base(base), m_stride(stride) {}
        row_iterator &operator++() {
            m_base += m_stride;
            return *this;
        }
        
        odmatrix<T> operator*() {
            return odmatrix<T>(m_base, m_stride);
        }
        
        bool operator==( const row_iterator &other ) {
            return m_base == other.m_base && m_stride == other.m_stride;
        }
        
        bool operator!=( const row_iterator &other ) {
            return !((*this) == other);
        }
        
    };

    
    tdmatrix() : m_base(0), m_asize(0), m_bsize(0), m_own(false) {}
    
    tdmatrix( size_t asize, size_t bsize ) : m_asize(asize), m_bsize(bsize) {
        m_base = new T[num_elements()];
        m_own = true;
        
        
    }
    
    
    
    ~tdmatrix() {
        if( m_own ) {
            delete[] m_base;
        }
    }
    
    const tdmatrix &swap( tdmatrix &other ) {
    	std::swap( m_base, other.m_base );
    	std::swap( m_own, other.m_own );
    	std::swap( m_asize, other.m_asize );
    	std::swap( m_bsize, other.m_bsize );

    	return *this;
	}

    void init_size( size_t asize, size_t bsize ) {
        // uuhm, is this the thing to do (TM) in the name of exception safety?
        
        tdmatrix<T> n(asize, bsize);
        
        //std::swap( *this, n );
        swap(n);
    }
    
    
    size_t num_elements() {
        return m_asize * m_bsize;
    }
    

    
    row_type operator[](ptrdiff_t a) {
        return row_type(m_base + (a * m_bsize), m_bsize );
    }
    
    const odmatrix<T> operator[](ptrdiff_t a) const {
        return row_type(m_base + (a * m_bsize), m_bsize );
    }
    
    row_iterator row_begin() {
        return row_iterator(m_base, m_bsize);
    }
    
    row_iterator row_end() {
        return row_iterator(m_base + num_elements(), m_bsize);
    }
    
    // WARNING: the default iterators work on single elements, not rows! This is inconsistend with operator[]
    iterator begin() {
        return m_base;
    }
    iterator end() {
        return m_base + num_elements();
    }
    
    size_t size() const {
        return m_asize;
    }
    


};
} // namespace ivy_mike


#endif
