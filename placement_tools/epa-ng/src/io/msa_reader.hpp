#pragma once

#include <memory>

#include "seq/MSA_Stream.hpp"
#include "seq/MSA_Info.hpp"
#include "io/Binary_Fasta.hpp"
#include "io/file_io.hpp"
#include "util/stringify.hpp"
#include "util/logging.hpp"
#include "util/Options.hpp"
#include "io/msa_reader_interface.hpp"

inline auto make_msa_reader(const std::string& file_name,
                            const MSA_Info& info,
                            const bool premasking = true)
{
  std::unique_ptr<msa_reader> result(nullptr);

  try {
    result = std::make_unique<Binary_Fasta_Reader>(file_name, info, premasking);
  } catch(const std::exception&) {
    LOG_DBG << "Failed to parse input as binary fasta (bfast), trying `fasta` instead.";
    result = std::make_unique<MSA_Stream>(file_name, info, premasking);
  }

  return result;
}
