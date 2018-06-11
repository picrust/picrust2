#pragma once

#include <string>
#include <fstream>
#include <streambuf>
#include <unordered_map>

#include "util/logging.hpp"

static std::string parse(const std::string& full, std::string qry, size_t& pos )
{
  std::string result;
  pos = full.find(qry, pos);
  if (pos != std::string::npos) {
    pos += qry.length();
    auto end = full.find('\n', pos);
    if (end == std::string::npos) {
      throw std::runtime_error{"couldnt find terminating newline?!"};
    } else {
      result = full.substr(pos, end - pos);
    }
    pos = end;
  }
  return result;
}

// using ssmap_t = std::unordered_map<std::string, std::string>;

// static const ssmap_t raxml8_to_raxmlng = {
//   {"GTR", "GTR"}
// };

// static std::string translate(const std::string& str, const ssmap_t& map)
// {
//   std::string result;

//   auto token = map.find(str);

//   if (token != map.end()) {
//     result = token->second;
//   } else {
//     throw std::runtime_error{std::string()+"couldnt find translation for token: "+str};
//   }
//   return result;
// }

static std::string from_raxml_8(const std::string& file)
{
  std::string model_desc;
  std::ifstream t(file);
  std::string full((std::istreambuf_iterator<char>(t)),
                    std::istreambuf_iterator<char>());
  size_t pos = 0;

  // parse data type
  auto type = parse(full, "DataType: ", pos);
  const bool dna = (type == "DNA");

  std::string sub_mat;

  // parse  sub matrix identifier
  sub_mat += parse(full, "Substitution Matrix: ", pos);

  if (not dna and sub_mat == "GTR") {
    sub_mat = "PROTGTR";
  }

  model_desc.append(sub_mat);

  // parse alpha
  auto alpha = parse(full, "alpha: ", pos);

  const std::string chars = dna ? std::string("ACGT") : std::string("ARNDCQEGHILKMFPSTWYV");

  // parse rates
  std::string rates;

  rates += "{";

  for (size_t i = 0; i < chars.length() - 1u; ++i) {
    for (size_t k = i + 1u; k < chars.length(); ++k) {
      if (k > 1) {
        rates += "/";
      }
      rates += parse(full,
        std::string("rate ") + chars[i] + " <-> " + chars[k] + ": ",
        pos
      );
    }
  }
  rates += "}";

  model_desc.append(rates);

  // parse stationary frequencies
  std::string freqs;

  freqs += "+FU{";
  for (size_t i = 0; i < chars.length(); ++i) {
    if (i > 0) {
        freqs += "/";
      }
    freqs += parse(full,
      std::string("freq pi(") + chars[i] +"): ",
      pos);
  }
  freqs += "}";

  model_desc.append(freqs);

  // append alpha
  alpha = "+G4{" + alpha + "}";
  model_desc.append(alpha);

  return model_desc;
}

inline std::string parse_model(const std::string& file)
{
  return from_raxml_8(file);
}