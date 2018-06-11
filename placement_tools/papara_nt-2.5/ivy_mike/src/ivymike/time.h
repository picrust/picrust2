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


#ifndef __IvyMike_time_h
#define __IvyMike_time_h

#include <iostream>
#include <numeric>
#include <stdexcept>

namespace ivy_mike {

double gettime(void);


class timer {
    double m_st;
    
public:
    timer() : m_st(gettime()) {};
    
    double elapsed() const {
        return gettime() - m_st;
    }
    
    
};


class perf_timer {
    const static size_t n_int = 16;
    size_t cur_int;
    
    double base_ticks;
    double int_ticks[n_int];
    
    static double my_getticks();
    bool m_quiet;
    
public:  
    perf_timer() : cur_int(0), m_quiet(false) {
        base_ticks = my_getticks();
        add_int();
    }
    
    perf_timer(bool quiet) : cur_int(0), m_quiet(quiet) {
        base_ticks = my_getticks();
        add_int();
    }
    
    void add_int() {
        if( cur_int >= n_int - 1 ) {
            return;
        }
        int_ticks[cur_int] = my_getticks() - base_ticks;
        ++cur_int;
    }
    
    void print() {
        
        if( !m_quiet ) {
            
                       
            double sum_nticks = 0.0;
            
            for( size_t i = 1; i < cur_int; ++i ) {
                sum_nticks += int_ticks[i] - int_ticks[i-1];
            }
            for( size_t i = 1; i < cur_int; ++i ) {
                double d = int_ticks[i] - int_ticks[i-1];
                
                std::cout << "int " << i - 1 << ": " << d << " (" << d / sum_nticks * 100 <<  "%)\n";
            }
        }
    }
    
    perf_timer &operator+=(const perf_timer &other ) {
        if( cur_int == 1 ) {
            // WARNING; brainfuck ahead
            bool quiet = m_quiet;
            *this = other;
            m_quiet = quiet;
        } else {
        
            if( n_int != other.n_int ) {
                throw std::runtime_error( "n_int != other.n_int\n" );
                return *this;
            }
            
            for( size_t i = 0; i < n_int; ++i ) {
                int_ticks[i] += other.int_ticks[i];
            }
        }
        return *this;
    }
    
};

}
#endif
