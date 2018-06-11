#include "util/stringify.hpp"

std::string stringify(raxml::Model& model)
{
  std::ostringstream output;

  output << "Substitution Matrix Symmetries: " << NEWL;
  output << stringify(model.submodel(0).rate_sym()) << NEWL;

  output << "Base Frequencies: " << NEWL;
  output << stringify(model.base_freqs(0)) << NEWL;

  output << "Substitution Rates: " << NEWL;
  output << stringify(model.subst_rates(0)) << NEWL;

  output << "Alpha: " << model.alpha() << NEWL;

  return output.str();
}

std::vector<std::string> split_by_delimiter(const std::string & text, const std::string delim)
{
  std::vector<std::string> parts;
  size_t start = 0;
  size_t end = 0;

  do
  {
    end = text.find(delim, start);
    end = (end != std::string::npos) ? end : text.length();
    parts.emplace_back(text.substr(start, end - start));
    start = end + delim.length();
  } while (end != std::string::npos and start <= text.length());

  return parts;
}
