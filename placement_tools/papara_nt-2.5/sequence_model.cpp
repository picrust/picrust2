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


#include "sequence_model.h"


using sequence_model::model;
using sequence_model::tag_aa;
using sequence_model::tag_dna;
using sequence_model::tag_dna4;

namespace raxml_aa_meaning {
// is it officially legal to initialize const static members in the header? I guess c++ removes redundant definitions during linking...
const static char inverse[23] = {'A','R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', '-'};
const static unsigned int bitVector[23] = {1, 2, 4, 8, 16, 32, 64, 128,
                                      256, 512, 1024, 2048, 4096,
                                      8192, 16384, 32768, 65536, 131072, 262144,
                                      524288, 12 /* N | D */, 96 /*Q | E*/, 1048575 /* - */ };
}

const std::vector<char> model<tag_aa>::inverse_meaning(raxml_aa_meaning::inverse, raxml_aa_meaning::inverse + ivy_mike::arrlen(raxml_aa_meaning::inverse));
const std::vector<unsigned int> model<tag_aa>::bit_vector(raxml_aa_meaning::bitVector, raxml_aa_meaning::bitVector + ivy_mike::arrlen(raxml_aa_meaning::bitVector));



namespace raxml_dna_meaning {
const char inverse[16]   = {'_', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', '-' };
//const char bitvector[17]   = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 15 };
}





const std::vector<char> model<tag_dna>::inverse_meaning(raxml_dna_meaning::inverse, raxml_dna_meaning::inverse + ivy_mike::arrlen(raxml_dna_meaning::inverse));
//const std::vector<uint8_t> model<tag_dna>::bit_vector(raxml_dna_meaning::bitvector, raxml_dna_meaning::bitvector + ivy_mike::arrlen(raxml_dna_meaning::bitvector));


namespace raxml_dna4_meaning {
const char inverse[5]   = {'A', 'C', 'G', 'T', '-' };
//const char bitvector[17]   = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 15 };
}


const std::vector<char> model<tag_dna4>::inverse_meaning(raxml_dna4_meaning::inverse, raxml_dna4_meaning::inverse + ivy_mike::arrlen(raxml_dna4_meaning::inverse));
