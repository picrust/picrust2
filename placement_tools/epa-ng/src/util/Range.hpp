#pragma once

#include <string>

class Range {
public:
  Range (const unsigned int begin, const unsigned int span)
    : begin(begin), span(span) {};
  Range() = default;
  ~Range () = default;

  unsigned int begin;
  unsigned int span;

private:

};

inline std::ostream& operator << (std::ostream& out, Range const& rhs)
{
  out << " begin " + std::to_string(rhs.begin);
  out << " span " + std::to_string(rhs.span);
  return out;
}

/*  Returns the range of a sequence outside of which there are ONLY indel characters.
 *  Range starts at the first valid position and ends after <span> characters, where
 *  begin + span is the first element not included in the valid range.
 *  Example:
 *  -  -  -  A  T  A  G  C  T  -  -
 *  0  1  2  3  4  5  6  7  8  9 10
 *  Output: (3,6)
 */
inline Range get_valid_range(const std::string& sequence)
{
  unsigned int lower = 0;
  unsigned int upper = sequence.length();

  while(sequence.c_str()[lower] == '-')
    lower++;

  while(sequence.c_str()[upper - 1] == '-')
    upper--;

  return Range(lower, upper - lower);
}


