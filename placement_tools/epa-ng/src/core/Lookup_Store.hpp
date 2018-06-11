#pragma once

#include <mutex>
#include <vector>
#include <map>
#include <utility>
#include <limits>
#include <cassert>

#include "util/Matrix.hpp"
#include "util/maps.hpp"
#include "util/Range.hpp"

constexpr size_t INVALID = std::numeric_limits<size_t>::max();

class Lookup_Store
{
/**
 * NOTE TO FUTURE DEVS:
 * This class has gotten a bit convoluted, so where is a brief overview of the various maps and tables:
 *
 * store_: vector of matrices, one per branch in the ref tree
 * <matrix in store>-> lookup_matrix: stores one CLV per character suitable for the model (ACGTVH- etc.)
 * char_map_: set of chars for which a lookup_matrix is done (see util/maps.hpp)
 * char_to_posish: maps ascii char to a column in a lookup_matrix. This also normalizes the input!
 *                 meaning: map upper and lowercase to the same CLV site, different variants of
 *                 GAP (-?Xx etc.) and ANY (N), U into T (RNA support) and defines invalid chars
 */
public:
  using lookup_type = Matrix<double>;

  Lookup_Store(const size_t num_branches, const size_t num_states) 
    : branch_(num_branches)
    , store_(num_branches)
    , char_map_size_((num_states == 4) ? NT_MAP_SIZE : AA_MAP_SIZE)
    , char_map_((num_states == 4) ? NT_MAP : AA_MAP)
  {
    const bool dna = (num_states == 4);

    for (size_t i = 0; i < 128; ++i) {
      char_to_posish_[i] = INVALID;
    }

    // build reverse map from char to ID in the charmap
    for (size_t i = 0; i < char_map_size_; ++i) {
      char_to_posish_[char_map_[i]] = i;
      char_to_posish_[std::tolower(char_map_[i])] = i;
    }

    if (dna) {
      // allow for RNA
      char_to_posish_['U'] = char_to_posish_['T'];
      char_to_posish_['u'] = char_to_posish_['T'];
    }

    // gap/any chars
    if (dna) {
      char_to_posish_['X'] = char_to_posish_['-'];
      char_to_posish_['x'] = char_to_posish_['-'];
      char_to_posish_['O'] = char_to_posish_['-'];
      char_to_posish_['o'] = char_to_posish_['-'];
      char_to_posish_['.'] = char_to_posish_['-'];
    } else {
      char_to_posish_['X'] = char_to_posish_['N'];
      char_to_posish_['x'] = char_to_posish_['N'];
    }
    char_to_posish_['?'] = char_to_posish_['-'];
  }

  Lookup_Store()  = delete;
  ~Lookup_Store() = default;

  void init_branch(const size_t branch_id, std::vector<std::vector<double>> precomps)
  {
    store_[branch_id] = Matrix<double>(precomps[0].size(), char_map_size_);

    for(size_t ch = 0; ch < precomps.size(); ++ch) {
      for(size_t site = 0; site < precomps[ch].size(); ++site) {
        store_[branch_id](site, ch) = precomps[ch][site];
      }
    }
  }

  std::mutex& get_mutex(const size_t branch_id)
  {
    return branch_[branch_id];
  }

  bool has_branch(const size_t branch_id) 
  {
    return store_[branch_id].size() != 0; 
  }

  lookup_type& operator[](const size_t branch_id)
  {
    return store_[branch_id];
  }

  unsigned char char_map(const size_t i)
  {
    if (i >= char_map_size_) {
      throw std::runtime_error{
        std::string("char_map access out of bounds! i =") + std::to_string(i)
      };
    }

    return char_map_[i];
  }

  size_t char_map_size()
  {
    return char_map_size_;
  }

  size_t char_position(unsigned char c) const
  {
    auto pos = char_to_posish_[c];

    if (pos == INVALID) {
      throw std::runtime_error{std::string("char is invalid! char = ") + std::to_string(c)};
    }

    return pos;
  }

  double sum_precomputed_sitelk(const size_t branch_id, const std::string& seq, const Range& range) const
  {
    assert(seq.length() == store_[branch_id].rows());
    
    double sum = 0;
    const auto& lookup_matrix = store_[branch_id];
    const auto& lookup = lookup_matrix.get_array();

    // unrolled loop
    size_t site = range.begin;
    const size_t end = range.begin + range.span;

    const size_t stride = 4;
    for (; site + stride-1u < end; site+=stride) {
      double sum_one =
      lookup[lookup_matrix.coord(site, char_to_posish_[seq[site]])]
      + lookup[lookup_matrix.coord(site+1u, char_to_posish_[seq[site+1u]])];
      
      double sum_two =
      lookup[lookup_matrix.coord(site+2u, char_to_posish_[seq[site+2u]])]
      + lookup[lookup_matrix.coord(site+3u, char_to_posish_[seq[site+3u]])];

      sum_one += sum_two;

      sum += sum_one;
    }

    // rest of the horizontal add
    while (site < end) {
      sum += lookup[lookup_matrix.coord(site, char_to_posish_[seq[site]])];
      ++site;
    }
    return sum;
  }

private:
  std::vector<std::mutex> branch_;
  std::vector<lookup_type> store_;
  const size_t char_map_size_;
  const unsigned char * char_map_;
  std::array<size_t, 128> char_to_posish_;
};
