#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>

#include "util/Matrix.hpp"
#include "util/maps.hpp"



class FourBit
{
  // shorthand
  using uchar = unsigned char;
private:
  static constexpr uchar NONE_CHAR = '-';

  uchar pack_(const uchar lhs, const uchar rhs)
  {
    assert(lhs < 16 and rhs < 16);
    return (lhs << 4) | rhs;
  }

  static inline std::pair<uchar, uchar> unpack_(const uchar c)
  {
    static const uchar upper_mask = 0b11110000;
    static const uchar lower_mask = 0b00001111;
    return {(c & upper_mask) >> 4, c & lower_mask};
  }

public:
  FourBit()
    : to_fourbit_(128, 128)
  {
    static_assert(NT_MAP_SIZE == 16, "Weird NT map size, go adjust encoder code!");

    // both chars are valid characters:
    const auto map_size = static_cast<uchar>(NT_MAP_SIZE);
    for (uchar i = 0; i < map_size; ++i) {
      assert(NT_MAP[i] == std::toupper(NT_MAP[i]));
      for (uchar j = 0; j < map_size; ++j) {
        uchar packed_char = pack_(i,j);

        // uppercase
        auto row = NT_MAP[i];
        auto col = NT_MAP[j];
        to_fourbit_.at(row, col) = packed_char;
        // lowercase
        row = std::tolower(NT_MAP[i]);
        to_fourbit_.at(row, col) = packed_char;

        col = std::tolower(NT_MAP[j]);
        to_fourbit_.at(row, col) = packed_char;

        row = NT_MAP[i];
        to_fourbit_.at(row, col) = packed_char;
      }
    }

    // valid + padding
    for (uchar i = 0; i < map_size; ++i) {
      uchar packed_char = pack_(i,0);

      // uppercase
      auto row = NT_MAP[i];
      auto col = NONE_CHAR;
      to_fourbit_.at(row, col) = packed_char;

      row = std::tolower(NT_MAP[i]);
      to_fourbit_.at(row, col) = packed_char;
    }

    // make the reverse lookup
    for (size_t i = 0; i < from_fourbit_.max_size() * 2; i+=2) {
      auto pair = unpack_(i/2);
      reinterpret_cast<uchar*>(from_fourbit_.data())[i] = 
        NT_MAP[pair.first];
      reinterpret_cast<uchar*>(from_fourbit_.data())[i+1u] = 
        NT_MAP[pair.second];
    }
  }
  ~FourBit() = default;

  inline size_t packed_size(const size_t size)
  {
    return std::ceil(size / 2.0);
  }

  // conversion functions
  inline std::basic_string<char> to_fourbit(const std::string& s)
  {
    const size_t p_size = packed_size(s.size());
    std::basic_string<char> res;
    res.reserve(p_size);

    size_t i = 0;
    for (; i + 1 < s.size(); i += 2) {
      res.push_back(to_fourbit_.at(s[i], s[i+1u]));
    }

    // original string size not divisible by 2: trailing padding
    if (i < s.size()) {
        res.push_back(to_fourbit_.at(s[i], NONE_CHAR));
    }

    return res;
  }

  std::string from_fourbit(const std::basic_string<char>& s, const size_t n)
  {
    assert(s.size() > 0);
    assert(n > 0);

    // prepare the result string
    std::string res;
    res.resize(n);

    // determine whether the packed string has padding
    const bool padded = (s.size() * 2 < n);

    // unpack
    size_t i = 0;
    for (; i < s.size() - 1u; ++i) {
      reinterpret_cast<char16_t*>(&res[0])[i] 
        = from_fourbit_[static_cast<uchar>(s[i])];
    }

    // last element is special
    auto char_pair = unpack_(static_cast<uchar>(s[i]));
    res[i*2] = NT_MAP[char_pair.first];

    if (not padded) {
      res[i*2+1] = NT_MAP[char_pair.second];
    } else {
      assert(NT_MAP[char_pair.second] == NONE_CHAR);
    }

    assert(res.size() == n);

    return res;
  }
  
private:
  Matrix<char> to_fourbit_;
  std::array<char16_t, 256> from_fourbit_;
};
