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

#ifndef __LARGEPHYLIP_H
#define __LARGEPHYLIP_H


#include <cstdio>
#include <cstdlib>
//#include <unistd.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>
//#include <sys/mman.h>
#include <fcntl.h>
//#include <unistd.h>
#include <cassert>

#include <stdint.h>

#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <algorithm>
#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>


typedef unsigned char u1_t;
// namespace gnu = __gnu_cxx;

struct Rec {
    off_t name;
    int nameLen;
    int nameMax;

    off_t data;
    int dataLen;
    int dataMax;

    inline std::string getName( u1_t *base ) {
//         char *tmp = (char*)alloca(nameLen+1);
//         memcpy( tmp, base + name, nameLen );
//         tmp[nameLen] = 0;
// 
//         return std::string(tmp);

	return std::string( base + name, base + name + nameLen );
    }

    inline std::string getData( u1_t *base ) {
//         char *tmp = (char*)alloca(dataLen+1);
//         memcpy( tmp, base + data, dataLen );
//         tmp[dataLen] = 0;

        
	return std::string( base + data, base + data + dataLen );
    }

};

class LargePhylip {
//    int m_fd;
    boost::interprocess::file_mapping m_fm;
    boost::interprocess::mapped_region m_mapping;

    off_t      m_fileSize;

    u1_t *m_buf;

    int m_nTaxa;
    int m_seqLen;
    int m_maxNameLen;

    std::vector<Rec> m_recs;

    std::map<std::string,size_t> m_nameMap;


    void interpret( off_t line, off_t lineLen, Rec& rec, size_t seq_len = -1 ) ;



public:
    LargePhylip( const char *filename ) ;

    ~LargePhylip() ;
    void print() ;
    void map() ;

    void unmap() ;

    int getIdx( const char *name ) ;
    
    inline std::string getName( int i ) {
	assert( m_buf != 0 );
	return m_recs.at(i).getName(m_buf);
    }
    
    inline size_t getNameLen( int i ) {
        return m_recs.at(i).nameLen;
    }
    
    inline std::string getSequence( int i ) {
	assert( m_buf != 0 );
	return m_recs.at(i).getData(m_buf);
    }
    
    inline u1_t * getSequenceBegin( int i ) {
        const Rec &rec = m_recs.at(i);
        return m_buf + rec.data;
    }
    
    inline u1_t * getSequenceEnd( int i ) {
        const Rec &rec = m_recs.at(i);
        return m_buf + rec.data + rec.dataLen;
    }
    
    inline int size() {
	return m_recs.size();
    }
    
    inline int getSeqLen() {
	return m_seqLen;
    }
    void getSequenceStart ( int i );
    
    inline size_t max_name_len() const {
        return m_maxNameLen;
    }
};

#endif