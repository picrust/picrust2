#include "seq/MSA_Info.hpp"

#include "io/Binary_Fasta.hpp"

MSA_Info make_msa_info(const std::string& file_path)
{
  MSA_Info info;
  try {
    info = Binary_Fasta::get_info(file_path);
  } catch(const std::exception&) {
    info = MSA_Info(file_path);
  }
  return info;
}