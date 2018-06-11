#pragma once

#include "seq/MSA_Stream.hpp"
#include "seq/MSA_Info.hpp"
#include "util/Options.hpp"
#include "tree/Tree.hpp"
#include "core/raxml/Model.hpp"

#include <string>

void simple_mpi(Tree& tree,
                const std::string& query_file,
                const MSA_Info& msa_info,
                const std::string& outdir,
                const Options& options,
                const std::string& invocation);

