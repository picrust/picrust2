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

#ifndef __ivy_mike__concurrent_h
#define __ivy_mike__concurrent_h


//
// the 'lock-free spin' barrier from
// S.A. Berger, A. Stamatakis: "Assessment of Barrier Implementions for
// Fine-Grain Parallel Regions on Current Multi-core Architectures"
// IEEE Cluster workshop on Parallel Programming and Applications on Accelerator Clusters,
// Heraklion, Greece, September 2010.
//

#include <cstring>

namespace ivy_mike {

template<int PAD>
class lockfree_spin_barrier {

	lockfree_spin_barrier( const lockfree_spin_barrier &other ) {}
	const lockfree_spin_barrier &operator=(const lockfree_spin_barrier &other ) { return *this; }

public:
	int m_num;
	volatile char *m_flags;

	void clear() {
		for (int i = 0; i < m_num; i++) {
			m_flags[i * PAD] = 0;
		}
	}

	lockfree_spin_barrier(int num) :
			m_num(num), m_flags(0) {
		alloc();

	}

	~lockfree_spin_barrier() {

		free( (void*)m_flags );
	}

	void alloc() {
#if 0
		m_flags = ((char *)memalign(64, num * PAD + 60)) + 60;
		//m_flags = ((char *)memalign(64, num * PAD + 40)) + 40;
#else

		//m_flags = (volatile char *)memalign(64, m_num * PAD );
		m_flags = (volatile char *) malloc(m_num * PAD);
#endif
		clear();
	}

	inline void wait(const int slot) {

		if (slot == 0) {

			int sum = 0;
			do {
				sum = 0;

				for (int i = 0; i < m_num * PAD; i += PAD) {
					sum += m_flags[i];
				}
// 	    printf( "sum: %d %d\n", slot, sum );
			} while (sum != m_num - 1);
			clear();
		} else {
			m_flags[slot * PAD] = 1;
		}
	}

	const char *id() {
		return "lockfree_spin_barrier";
	}
};

}
#endif
