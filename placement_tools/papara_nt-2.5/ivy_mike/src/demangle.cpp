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


#include "ivymike/demangle.h"
#include "ivymike/stupid_ptr.h"


#ifndef WIN32
#include <cxxabi.h>
#else
//#include <Windows.h>
//#include <DbgHelp.h>
#endif

namespace ivy_mike {


    
#ifndef WIN32


std::string demangle( const char *tname ) {
    int status;
    char *realname = abi::__cxa_demangle(tname, 0, 0, &status);
    freeer f(realname);
    
    return std::string( realname );
}
#else
std::string demangle( const char *tname ) {
#if 0
	// doing this acually requires linking an extra library with the
	// stupid name 'Dbghelp.dll'. It's not worth it...

	const DWORD len = 256;
	char name[len];
	
	UnDecorateSymbolName( tname, name, len, 0 );

	return std::string( name );
#else

	return std::string( tname );
#endif
}
#endif
}
