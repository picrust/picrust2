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

// boost mmap stuff does not compile on android, yet
#if !defined (__ANDROID__) && !defined(__native_client__)

#include "ivymike/LargePhylip.h"
#include <fstream>


LargePhylip::LargePhylip(const char* filename) 
    : m_fm(filename, boost::interprocess::read_only), 
    m_fileSize(0), 
    m_buf(0), 
    m_nTaxa(0), 
    m_seqLen(0), 
    m_maxNameLen(0) 
{
    //m_fd = open( filename, O_RDONLY );
    
    {
        std::ifstream ist(filename);
        ist.seekg(0, std::ios::end);
        m_fileSize = ist.tellg();
    }
    

//     printf( "size: %zd\n", m_fileSize );


    map();
	
    off_t ptr = 0;

    bool haveHeader = false;
    std::vector<Rec>::iterator currec; // When all you have is a hammer, everything looks like a nail
    
    size_t seq_len = -1;
    while ( ptr < m_fileSize ) {
        off_t spos = ptr;

        if( seq_len != size_t(-1)) {
            ptr+=seq_len;
        }
        // FIXME: this assumes left to right execution of the expression!
        while ( ptr < m_fileSize && m_buf[ptr] != '\n' ) {
            ptr++;
        }
        // ptr currently points to the newline (or the first position past the EOF).
        // in the non-EOF case this is one position before the current mbuf.position()

//             u1_t *line = &m_buf[spos];
        
        size_t lineLen = ptr - spos;
        // advance ptr past newline (=beginning of next line)
        ptr++;
    
        if ( !haveHeader ) {
            Rec rec;
//                 printf( "ll: %d %d %d\n", lineLen, spos, ptr );
            interpret( spos, lineLen, rec );

            std::string name = rec.getName(m_buf);
            std::string data = rec.getData(m_buf);

//                 printf( "header: %s %s\n", name.c_str(), data.c_str() );
            haveHeader = true;

            m_nTaxa = atoi( name.c_str() );
            m_seqLen = atoi( data.c_str() );
            seq_len = m_seqLen;
            m_recs.resize( m_nTaxa );
            currec = m_recs.begin();
        } else {
            if ( lineLen > 0 ) {
                interpret(spos, lineLen, *currec, seq_len );
                std::string n = (*currec).getName(m_buf);

                size_t idx = currec - m_recs.begin();
                if ( m_nameMap.find( n ) != m_nameMap.end() ) {
                    std::stringstream out;
                    out << idx;
                    n = out.str();
                }
                m_nameMap[n] = idx;
                m_maxNameLen = std::max( m_maxNameLen, int(n.length()) );
                ++currec;
                
                
            }
        }
    }

    assert( currec == m_recs.end() );
//         print();

}
LargePhylip::~LargePhylip() {
    if ( m_buf != 0 ) {
        unmap();
    }
}
void LargePhylip::print() {
    for ( std::vector< Rec >::iterator it = m_recs.begin(); it != m_recs.end(); ++it ) {
        printf( "name: %s %d\n", (*it).getName(m_buf).c_str(), (*it).dataLen );
    }
}
void LargePhylip::map() {
	
    assert( m_buf == 0 );
    assert(m_fileSize > 0 );
#if 0
    m_buf = (u1_t*)mmap( 0, m_fileSize, PROT_READ, MAP_PRIVATE, m_fd, 0 );
#endif
    m_mapping = boost::interprocess::mapped_region( m_fm, boost::interprocess::read_only );
//     assert(0);
	m_buf = (u1_t*) m_mapping.get_address();
	//madvise(m_buf, m_fileSize, MADV_RANDOM );
    assert( m_buf != 0 );
	
}
void LargePhylip::unmap() {
    assert(m_buf != 0);
    assert(m_fileSize > 0 );
//    munmap(m_buf, m_fileSize);
    m_mapping = boost::interprocess::mapped_region();
    m_buf = 0;
}
int LargePhylip::getIdx(const char* name) {
    std::map< std::string, size_t >::iterator it = m_nameMap.find(std::string(name));
    if ( it != m_nameMap.end() ) {
        return (*it).second;
    } else {
        return -1;
    }
}
void LargePhylip::interpret(off_t line, off_t lineLen, Rec& rec, size_t seq_len) {
    rec.name = line;

    off_t ptr = 0;

        
    while ( ptr < lineLen && !isspace(*(m_buf + line + ptr) )) {

        ptr++;
    }

    assert( ptr > 0 );

    rec.nameLen = ptr;
    rec.nameMax = ptr;


    int nspace = 0;
    
//     if( seq_len != size_t(-1)) {
//         ptr+=seq_len;
//     }
//     assert( pre < lineLen );
    while ( ptr < lineLen && isspace(*(m_buf + line + ptr) )) {
        ptr++;
        nspace++;
    }

    assert( nspace > 0 );

    rec.data = line + ptr;
    rec.dataLen = lineLen - ptr;
    rec.dataMax = rec.dataLen;

}
#endif
