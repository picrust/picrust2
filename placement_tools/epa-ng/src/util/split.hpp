#include "genesis/genesis.hpp"

#include "util/logging.hpp"

void split(std::string ref_msa, std::vector<std::string> query_files, std::string outdir="")
{
    auto outfile = outdir + "query.fasta";

    if (genesis::utils::file_exists(outfile)) {
      throw std::runtime_error{outfile + " already exists!"};
    }

    genesis::sequence::SequenceSet ref_set;
  
    auto reader = genesis::sequence::PhylipReader();
    reader.from_file( ref_msa, ref_set );

    auto ref_labels = labels( ref_set );

    genesis::sequence::SequenceSet qry_set;
    for (const auto& f : query_files) {
      LOG_INFO << "File: " << f;
      reader.from_file( f, qry_set );
    }

    genesis::sequence::filter_by_label_list( qry_set, ref_labels );

    auto writer = genesis::sequence::FastaWriter();
    LOG_INFO << "Writing output: " << outfile;
    writer.to_file( qry_set, outfile );
}