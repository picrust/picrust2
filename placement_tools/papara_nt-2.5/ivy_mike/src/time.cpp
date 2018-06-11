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



#ifndef WIN32
#include <sys/time.h>
#else 
#include "ivymike/disable_shit.h"

#include <windows.h>
static bool g_pc_valid = false;
static LARGE_INTEGER g_pc_freq;
#endif
#include "ivymike/time.h"
#include "ivymike/cycle.h"


double ivy_mike::gettime(void )
{
#ifdef WIN32
    static bool g_pc_valid = false;
    static LARGE_INTEGER g_pc_freq;

	if( !g_pc_valid ) {
		QueryPerformanceFrequency( &g_pc_freq );
		g_pc_valid = true;
	} 
	LARGE_INTEGER pcv;
	QueryPerformanceCounter( &pcv );
	return pcv.QuadPart / double(g_pc_freq.QuadPart);
	
#else
 struct timeval ttime;
 gettimeofday(&ttime , 0);
 return ttime.tv_sec + ttime.tv_usec * 0.000001;
#endif
}


double ivy_mike::perf_timer::my_getticks() {

#ifdef HAVE_TICK_COUNTER
#if defined(WIN32) && !defined(_M_X64)
	return double(getticks().QuadPart);
#else
	return double(getticks());
#endif
#else
    // fall back to msecs (kind of)
    return ivy_mike::gettime() * 1000;
#endif

}

