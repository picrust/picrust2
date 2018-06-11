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


#ifndef __multiple_alignment_h
#define __multiple_alignment_h

#include <string>
#include <vector>
#include <iosfwd>
#include <stdint.h>

namespace ivy_mike {
	
struct multiple_alignment {
	std::vector <std::string > names;
	std::vector <std::vector<uint8_t> > data;
    
	
	bool load_phylip( std::istream &is );
	bool load_phylip( const char *name );
	
};
	


}
#endif
