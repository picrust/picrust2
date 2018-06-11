#pragma once

#include "util/template_magic.hpp"

// map of all IUPAC DNA one letter codes.
// additionally, position is the encoding in 4bit format:
// ACGT
// 0001 <- for example for T, index = 1
constexpr unsigned char NT_MAP[] = 
  {
 '-', // 0000 = 0
 'T', // 0001 = 1
 'G', // 0010 = 2
 'K', // 0011 = 3
 'C', // 0100 = 4
 'Y', // 0101 = 5
 'S', // 0110 = 6
 'B', // 0111 = 7
 'A', // 1000 = 8
 'W', // 1001 = 9
 'R', // 1010 = 10
 'D', // 1011 = 11
 'M', // 1100 = 12
 'H', // 1101 = 13
 'V', // 1110 = 14
 'N'  // 1111 = 15
};
constexpr unsigned char AA_MAP[] = 
  { 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 
    'T', 'V', 'W', 'Y', '-', 'X', 'B', 'Z'};

constexpr size_t NT_MAP_SIZE = array_size(NT_MAP);
constexpr size_t AA_MAP_SIZE = array_size(AA_MAP);
