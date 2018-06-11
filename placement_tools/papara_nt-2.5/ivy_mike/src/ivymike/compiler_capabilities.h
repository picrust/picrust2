/*
 * Copyright (C) 2012 Simon A. Berger
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


#ifndef ivy_mike__compiler_capabilities_h__
#define ivy_mike__compiler_capabilities_h__


// WARNING: assuming that the 'November 2012 Compiler CTP' is used when on visual studio 2012
#if __cplusplus >= 201103L || _MSC_VER >= 1700
  #define IVY_MIKE__USE_CPP11 (1)
#endif 


#endif