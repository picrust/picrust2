/*
    gappa - Genesis Applications for Phylogenetic Placement Analysis
    Copyright (C) 2017-2018 Lucas Czech and HITS gGmbH

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact:
    Lucas Czech <lucas.czech@h-its.org>
    Exelixis Lab, Heidelberg Institute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "commands/prepare/chunkify.hpp"

#include "options/global.hpp"
#include "tools/version.hpp"

#include "CLI/CLI.hpp"

#include "genesis/sequence/sequence.hpp"
#include "genesis/sequence/sequence_set.hpp"
#include "genesis/sequence/formats/fasta_input_iterator.hpp"
#include "genesis/sequence/formats/fasta_writer.hpp"
#include "genesis/sequence/functions/labels.hpp"

#include "genesis/utils/core/std.hpp"
#include "genesis/utils/io/input_source.hpp"
#include "genesis/utils/io/input_stream.hpp"
#include "genesis/utils/io/output_stream.hpp"
#include "genesis/utils/text/string.hpp"
#include "genesis/utils/tools/md5.hpp"
#include "genesis/utils/tools/sha1.hpp"
#include "genesis/utils/tools/sha256.hpp"

#ifdef GENESIS_OPENMP
#   include <omp.h>
#endif

#include <sparsepp/spp.h>

#include <cstdio>
#include <fstream>
#include <unordered_map>
#include <unordered_set>

// =================================================================================================
//      Typedefs
// =================================================================================================

/**
 * @brief Store the data needed to write one abundace file.
 *
 * That is, for one sequence, we need the chunk it is in, and all abundances of the different
 * labels that this sequence has appeared with.
 */
struct SequenceInfo
{
    // In which chunk was this sequence stored?
    size_t chunk_num;

    // Which label has which abundance?
    std::unordered_map< std::string, size_t > abundances;
};

/**
 * @brief Map from hash to Sequence Info for storing per input file abundances and chunk nums.
 */
using AbundancesHashMap = std::unordered_map< std::string, SequenceInfo >;

// =================================================================================================
//      Setup
// =================================================================================================

void setup_chunkify( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto opt = std::make_shared<ChunkifyOptions>();
    auto sub = app.add_subcommand(
        "chunkify",
        "Chunkify a set of fasta files and create abundance maps."
    );

    // -----------------------------------------------------------
    //     Input options
    // -----------------------------------------------------------

    opt->sequence_input.add_fasta_input_opt_to_app( sub );

    // -----------------------------------------------------------
    //     Fill in custom options
    // -----------------------------------------------------------

    // Chunk Size
    sub->add_option(
        "--chunk-size",
        opt->chunk_size,
        "Number of sequences per chunk file.",
        true
    )->group( "Settings" );

    // Minimum Abundance
    sub->add_option(
        "--min-abundance",
        opt->min_abundance,
        "Minimum abundance of a single sequence. Sequences below are filtered out.",
        true
    )->group( "Settings" );

    // Hash Function
    sub->add_set_ignore_case(
        "--hash-function",
        opt->hash_function,
        { "SHA1", "SHA256", "MD5" },
        "Hash function for re-naming and identifying sequences.",
        true
    )->group( "Settings" );

    // -----------------------------------------------------------
    //     Output options
    // -----------------------------------------------------------

    opt->chunk_output.add_output_dir_opt_to_app( sub, "chunks" );
    opt->chunk_output.add_file_prefix_opt_to_app( sub, "chunk", "chunk_" );

    opt->abundance_output.add_output_dir_opt_to_app( sub, "abundances" );
    opt->abundance_output.add_file_prefix_opt_to_app( sub, "abundance", "abundances_" );

    // -----------------------------------------------------------
    //     Callback
    // -----------------------------------------------------------

    // Set the run function as callback to be called when this subcommand is issued.
    // Hand over the options by copy, so that their shared ptr stays alive in the lambda.
    sub->set_callback( [ opt ]() {
        run_chunkify( *opt );
    });
}

// =================================================================================================
//      Helpers
// =================================================================================================

void write_chunk_file(
    ChunkifyOptions const& options,
    genesis::sequence::SequenceSet const& chunk,
    size_t chunk_count
) {
    using namespace genesis::sequence;

    // Do not write a file if there is no content.
    if( chunk.empty() ) {
        return;
    }

    // Prepare fata writer.
    auto writer = FastaWriter();
    writer.line_length( 0 );

    // Generate output file name.
    auto const ofn
        = options.chunk_output.out_dir()
        + options.chunk_output.file_prefix()
        + std::to_string( chunk_count )
        + ".fasta"
    ;

    // Write
    writer.to_file( chunk, ofn );
}

void write_abundance_map_file(
    ChunkifyOptions const& options,
    AbundancesHashMap const& seq_abundances,
    size_t input_file_counter
) {
    using namespace genesis::utils;

    // Base name of the current input file
    auto const base_fn = options.sequence_input.base_file_name( input_file_counter );

    // Pprepare a new abundance file
    auto const fn
        = options.abundance_output.out_dir()
        + options.abundance_output.file_prefix()
        + base_fn
        + ".json"
    ;
    std::ofstream ofs;
    file_output_stream( fn, ofs );
    ofs << "{\n";

    // Write file metadata.
    ofs << "  \"sample\": \"" << base_fn << "\",\n";
    ofs << "  \"gappa\": \"" << gappa_version() << "\",\n";
    ofs << "  \"invocation\": \"" << global_options.command_line() << "\",\n";
    ofs << "  \"hash\": \"" << options.hash_function << "\",\n";

    // Write name of the input file for later identification.
    ofs << "  \"abundances\": [";

    // Write abundance information for this file.
    bool is_first_seq = true;
    for( auto seq_it = seq_abundances.begin(); seq_it != seq_abundances.end(); ++seq_it ) {

        // Print comma for all but the first entry.
        if( ! is_first_seq ) {
            ofs << ",";
        }
        is_first_seq = false;
        ofs << "\n";

        // Print sequence data.
        ofs << "    [ \"" << seq_it->first << "\", ";
        ofs << seq_it->second.chunk_num << ", {";

        // Write per label abundances.
        bool is_first_abun = true;
        for( auto const& label_abun : seq_it->second.abundances ) {

            // Print comma for all but the first entry.
            if( ! is_first_abun ) {
                ofs << ",";
            }
            is_first_abun = false;
            ofs << "\n";

            ofs << "      \"" << label_abun.first << "\": " << label_abun.second;
        }

        ofs << "\n    }]";
    }

    // Finish the file.
    ofs << "\n  ]\n";
    ofs << "}\n";
    ofs.close();
}

// =================================================================================================
//      Main Work Function
// =================================================================================================

template< class HashFunction >
void run_chunkify_with_hash( ChunkifyOptions const& options )
{
    using namespace genesis::utils;
    using namespace genesis::sequence;

    // using ChunkHashMap = std::unordered_map< typename HashFunction::DigestType, size_t >;
    using ChunkHashMap = spp::sparse_hash_map< typename HashFunction::DigestType, size_t >;

    // Sequences hashes, mapping to the chunk number where they are stored,
    // i.e. where they first occured.
    ChunkHashMap hash_to_chunk;

    // -----------------------------------------------------------
    //     Iterate Input Files
    // -----------------------------------------------------------

    // Collect sequences for the current chunk here.
    SequenceSet current_chunk;
    size_t file_count = 0;
    size_t chunk_count = 0;
    size_t total_seqs_count = 0;
    size_t min_abun_count = 0;
    auto const set_size = options.sequence_input.file_count();

    // Iterate fasta files
    #pragma omp parallel for schedule(dynamic)
    for( size_t fi = 0; fi < set_size; ++fi ) {
        auto const& fasta_filename = options.sequence_input.file_path( fi );

        // User output
        if( global_options.verbosity() >= 2 ) {
            #pragma omp critical(GAPPA_CHUNKIFY_PRINT_PROGRESS)
            {
                ++file_count;
                std::cout << "Processing file " << file_count << " of " << set_size;
                std::cout << ": " << fasta_filename << "\n";
            }
        }

        // Count identical sequences of this fasta file, accessed via their hash.
        AbundancesHashMap seq_abundances;

        // Iterate sequences
        auto it = FastaInputIterator( options.sequence_input.fasta_reader() );
        for( it.from_file( fasta_filename ); it; ++it ) {
            #pragma omp atomic
            ++total_seqs_count;

            // Check for min abundance.
            auto const abundance = guess_sequence_abundance( *it );
            if( abundance < options.min_abundance ) {
                continue;
            }
            #pragma omp atomic
            ++min_abun_count;

            // Calculate (relatively expensive) hashes.
            auto const hash_digest = HashFunction::from_string_digest( it->sites() );
            auto const hash_hex = HashFunction::digest_to_hex( hash_digest );

            // Increment seq abundance for this file and label.
            auto& seq_abun = seq_abundances[ hash_hex ];
            seq_abun.abundances[ it->label() ] += abundance;

            // The hash calculation above is the main work of this loop.
            // The rest is "just" setting some values (and the occasional chunk flush),
            // but we need a fully blown critical section for them.
            #pragma omp critical(GAPPA_CHUNKIFY_UPDATE_MAPS)
            {
                auto const hash_it = hash_to_chunk.find( hash_digest );
                if( hash_it != hash_to_chunk.end() ) {

                    // We saw that sequence before. Don't need to add it to the chunk,
                    // just use its chunk count for the current file.
                    seq_abun.chunk_num = hash_it->second;

                } else {

                    // New sequence: never saw that hash before. Add it to the chunk, store chunk num.
                    current_chunk.add( Sequence( hash_hex, it->sites() ));
                    hash_to_chunk[ hash_digest ] = chunk_count;
                    seq_abun.chunk_num = chunk_count;

                    // If a chunk is full, flush it.
                    if( current_chunk.size() >= options.chunk_size ) {
                        write_chunk_file( options, current_chunk, chunk_count );
                        ++chunk_count;
                        current_chunk.clear();
                    }
                }
            }
        }

        // Finished a fasta file. Write its abundances.
        write_abundance_map_file( options, seq_abundances, fi );
    }

    // -----------------------------------------------------------
    //     Finish
    // -----------------------------------------------------------

    // Write the remaining chunk.
    write_chunk_file( options, current_chunk, chunk_count );

    if( global_options.verbosity() >= 1 ) {
        std::cout << "Processed " << total_seqs_count << " sequences, thereof ";
        std::cout << (total_seqs_count - min_abun_count) << " (";
        std::cout << ( 100 * (total_seqs_count - min_abun_count) / total_seqs_count );
        std::cout << "%) filtered due to low abundance.\n";
        std::cout << "Wrote " << hash_to_chunk.size() << " unique sequences ";
        std::cout << "in " << ( chunk_count + 1 ) << " fasta chunk files.\n";
    }
}

// =================================================================================================
//      Run
// =================================================================================================

void run_chunkify( ChunkifyOptions const& options )
{
    using namespace genesis::utils;

    // -----------------------------------------------------------
    //     Input File Preparations
    // -----------------------------------------------------------

    // Check if any of the files we are going to produce already exists. If so, fail early.
    options.abundance_output.check_nonexistent_output_files(
        { options.abundance_output.file_prefix() + ".*\\.json" }
    );
    options.chunk_output.check_nonexistent_output_files(
        { options.chunk_output.file_prefix() + "[0-9]+\\.fasta" }
    );

    // Print some user output.
    options.sequence_input.print();

    // -----------------------------------------------------------
    //     Run
    // -----------------------------------------------------------

    if( options.hash_function == "SHA1" ) {
        run_chunkify_with_hash<SHA1>( options );
    } else if( options.hash_function == "SHA256" ) {
        run_chunkify_with_hash<SHA256>( options );
    } else if( options.hash_function == "MD5" ) {
        run_chunkify_with_hash<MD5>( options );
    } else {
        throw CLI::ConversionError( "Unknown hash function: " + options.hash_function );
    }

}
