/*
 * Copyright (C) 2009-2012 Simon A. Berger
 * 
 * This file is part of papara.
 * 
 *  papara is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  papara is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with papara.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <climits>
#include <cmath>
#include <time.h>
#include <cstdlib>
#include <cstdio>
#include <ctype.h>
#include <cassert>
#include <cstring>
#include <stdint.h>
#include <algorithm>
#include <vector>
#include <limits>
#include <stdexcept>

#include "pars_align_seq.h"
#include "sequence_model.h"

//pars_align_seq::pars_align_seq(const int* seqA, const unsigned char* seqB,
//		size_t n_a, size_t n_b, size_t aStride, const unsigned int* aAux,
//		size_t aAuxStride, pars_align_seq::arrays &arr,
//		const unsigned int *bvtrans, pars_align_seq::score_t gapOpen,
//		score_t gapExtend, pars_align_seq::score_t mismatch,
//		pars_align_seq::score_t matchCGap) :
//		GAP_OPEN(gapOpen), GAP_EXTEND(gapExtend), GAP_OPEN_EXTEND(
//				GAP_OPEN + GAP_EXTEND), MISMATCH_PENALTY(mismatch), MATCH_CGAP(
//				matchCGap), g_dump(false), m_ncups(0), m_arr(arr), m_tbStartA(
//				-1), m_tbStartB(-1), m_bvtrans(bvtrans)
//
//{
//	if (n_a <= n_b) {
//		throw std::runtime_error("n_a <= n_b\n");
//	}
//
//	m_na = n_a;
//	m_nb = n_b;
//
//	m_ma = n_a + 1;
//	m_mb = n_b + 1;
//
//	m_a = seqA;
//	m_b = seqB;
//
//	m_aAux = aAux;
//	m_aStride = aStride;
//	m_aAuxStride = aAuxStride;
//
//	m_msize = m_ma * m_mb;
//
//	m_arr.size(m_msize);
//
//}
//pars_align_seq::~pars_align_seq() {
//
//}
//void pars_align_seq::alignFreeshiftS11() {
//
//	if (m_na > m_nb) {
//		for (size_t ia = 0; ia <= m_na - m_nb - 1; ia++) {
//			//  std::cout << "init: " << ia << " " << m_na << " " << m_nb << "\n";
//
//			//size_t xxx = saddr ( ia, -1 );
//			//         if( xxx >= m_arr.size() ) {
//			//             std::cout << "meeeep: " << xxx << " " << ia << "\n";
//			//         }
//			m_arr.scoreL[saddr(ia, -1)] = 0;
//			m_arr.score[saddr(ia, -1)] = 0;
//		}
//	}
//
//	m_arr.score[saddr(-1, -1)] = 0;
//	m_arr.scoreL[saddr(-1, -1)] = 0;
//
//	for (size_t i = 0; i < m_nb; i++) {
//		m_arr.score[saddr(i - 1, i)] = LARGE_VALUE;
//		m_arr.scoreL[saddr(i - 1, i)] = LARGE_VALUE;
//		m_arr.scoreL[saddr(i, i)] = LARGE_VALUE;
//		m_arr.score[saddr((m_na - m_nb) + i, i - 1)] = LARGE_VALUE;
//	}
//
//	for (size_t i = 0; i < m_na; i++) {
//		m_arr.score[saddr(i, -1)] = 0;
//	}
//
//	if (m_arr.dir != 0) {
//
//		if (m_na > m_nb) {
//			for (size_t ia = 0; ia <= m_na - m_nb - 1; ia++) {
//				m_arr.dir[saddr(ia, -1)] = BT_STAY_L;
//			}
//		}
//	}
//}
//
//void pars_align_seq::tracebackCompressed(std::vector<uint8_t>& bvec) {
////     int idx_b = m_na-1;
//	ptrdiff_t ba = m_tbStartA;
//	ptrdiff_t bb = m_tbStartB;
//	bool inL = false;
//
//	if (m_tbStartA < ptrdiff_t(m_na) - 1) {
//		for (int ba = m_na - 1; ba > m_tbStartA; ba--) {
//			bvec.push_back(1);
//
//		}
//	}
//
//	// traceback
//
//	while ((ba >= 0 || bb >= 0)) {
//		if (!inL) {
//
//			// we are currently not in a gap. see if we should match ('go diagonal') or open a gap ('jump to gap-matrix')
//			if ((m_arr.dir[saddr(ba, bb)] & BT_STAY) != 0) {
//				bvec.push_back(0);
////                 idx_b--;
//				ba--;
//				bb--;
//			} else if ((m_arr.dir[saddr(ba, bb)] & BT_UP) != 0) {
//				bb--;
//				bvec.push_back(2);
//
//			} else {
//				// open gap. keep position and switch to gap-matrix
//				inL = true;
//			}
//		} else {
//			// we are in a gap
//
//			// output gap
//			bvec.push_back(1);
//
////             idx_b--;
//			// check if we should stay in the gap-matrix
//			inL = (m_arr.dir[saddr(ba, bb)] & BT_STAY_L) != 0;
//			ba--;
//		}
//	}
//
////     for( int i = 0; i <  bvec.size(); i++ ) {
////         printf( "tb %d: %d\n", i, bvec[i] );
////     }
//
//}
