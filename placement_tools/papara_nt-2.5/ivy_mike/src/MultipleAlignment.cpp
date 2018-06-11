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

#include "ivymike/MultipleAlignment.h"


MultipleAlignment* MultipleAlignment::loadPhylip(const char* name) {
    std::ifstream is(name);

    size_t nTaxon;
    size_t seqLen;

    is >> nTaxon;
    is >> seqLen;

    printf( "header: %zd %zd\n", nTaxon, seqLen );

    size_t n = 0;

    MultipleAlignment *ma = new MultipleAlignment();//nTaxon, seqLen);
    while ( !is.eof() ) {
        std::string name;
        std::string seq;

        is >> name;
        is >> seq;

        ma->m_names.push_back(name);
        ma->m_data.push_back(seq);
        //printf( "name: %s\n", name.c_str() );
        n++;
    }

    assert( n == nTaxon );
    printf( "n: %zd\n", n );

    return ma;
}
