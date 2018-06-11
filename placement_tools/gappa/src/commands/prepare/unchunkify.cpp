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

#include "commands/prepare/unchunkify.hpp"

#include "options/global.hpp"

#include "CLI/CLI.hpp"

#include "genesis/placement/formats/jplace_reader.hpp"
#include "genesis/placement/formats/jplace_writer.hpp"
#include "genesis/placement/sample.hpp"

#include "genesis/utils/containers/mru_cache.hpp"
#include "genesis/utils/core/fs.hpp"
#include "genesis/utils/formats/json/document.hpp"
#include "genesis/utils/formats/json/iterator.hpp"
#include "genesis/utils/formats/json/reader.hpp"
#include "genesis/utils/text/string.hpp"
#include "genesis/utils/tools/md5.hpp"
#include "genesis/utils/tools/sha1.hpp"
#include "genesis/utils/tools/sha256.hpp"

#ifdef GENESIS_OPENMP
#   include <omp.h>
#endif

#include <algorithm>
#include <cassert>
#include <limits>
#include <memory>
#include <utility>

// =================================================================================================
//      Typedefs
// =================================================================================================

/**
 * @brief We offer three modes, depending on whaty type of input is given.
 */
enum class UnchunkifyMode
{
    kNone,
    kChunkListFile,
    kChunkFileExpression,
    kJplaceInput
};

/**
 * @brief Store a sample, along with a map from sequence hash to the pquery index in the sample.
 *
 * Not all modes of the command use the map, it can thus be empty if not needed.
 */
template< class HashFunction >
struct MappedSample
{
    genesis::placement::Sample sample;
    std::unordered_map<typename HashFunction::DigestType, size_t> hash_to_index;
};

/**
 * @brief Store a sample index and a pquery index that tells where a particular hash can be found.
 */
struct SamplePqueryIndices
{
    size_t sample_index;
    size_t pquery_index;
};

/**
 * @brief Map of all sequences hashes to their sample and pquery indices.
 *
 * Needed for the jplace files input mode.
 */
template< class HashFunction >
using HashToIndexMap = std::unordered_map<typename HashFunction::DigestType, SamplePqueryIndices>;

/**
 * @brief Cache for chunk jplace files, mapping frmo file path to sample.
 */
template< class HashFunction >
using ChunkCache = genesis::utils::MruCache<std::string, std::shared_ptr<MappedSample<HashFunction>>>;

// =================================================================================================
//      Setup
// =================================================================================================

void setup_unchunkify( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto opt = std::make_shared<UnchunkifyOptions>();
    auto sub = app.add_subcommand(
        "unchunkify",
        "Unchunkify a set of jplace files using abundace map files and create per-sample jplace files."
    );

    // -----------------------------------------------------------
    //     Input options
    // -----------------------------------------------------------

    opt->abundance_map_input.add_multi_file_input_opt_to_app( sub, "abundances", "json" );
    opt->jplace_input.add_jplace_input_opt_to_app( sub, false )->group( "Input" );

    // -----------------------------------------------------------
    //     Fill in custom options
    // -----------------------------------------------------------

    // Chunk List file.
    auto chunk_list_file_opt = sub->add_option(
        "--chunk-list-file",
        opt->chunk_list_file,
        "If provided, needs to contain a list of chunk file paths in the numerical order that was "
        "produced by the chunkify command."
    )->group( "Input" );

    // Chunk List file.
    auto chunk_file_expression_opt = sub->add_option(
        "--chunk-file-expression",
        opt->chunk_file_expression,
        "If provided, the expression is used to load jplace files by replacing any '@' character "
        "with the chunk number."
    )->group( "Input" );

    // Cache size
    sub->add_option(
        "--jplace-cache-size",
        opt->jplace_cache_size,
        "Cache size to determine how many jplace files are kept in memory. Default (0) means all. "
        "Use this if the command runs out of memory. It however comes at the cost of longer runtime. "
        "In order to check how large the cache size can be, you can run the command with more verbosity "
        "(-vv), which will report the used cache size until it crashes. Then, set the cache size to "
        "something below that.",
        true
    )->group( "Settings" );

    // Hash Function
    sub->add_set_ignore_case(
        "--hash-function",
        opt->hash_function,
        { "SHA1", "SHA256", "MD5" },
        "Hash function that was used for re-naming and identifying sequences in the chunkify command.",
        true
    )->group( "Settings" );

    // Make the three input modes mutually exclusive.
    chunk_list_file_opt->excludes( opt->jplace_input.option() );
    chunk_list_file_opt->excludes( chunk_file_expression_opt );
    chunk_file_expression_opt->excludes( opt->jplace_input.option() );
    chunk_file_expression_opt->excludes( chunk_list_file_opt );
    opt->jplace_input.option()->excludes( chunk_list_file_opt );
    opt->jplace_input.option()->excludes( chunk_file_expression_opt );

    // -----------------------------------------------------------
    //     Output options
    // -----------------------------------------------------------

    opt->file_output.add_output_dir_opt_to_app( sub );

    // -----------------------------------------------------------
    //     Callback
    // -----------------------------------------------------------

    // Set the run function as callback to be called when this subcommand is issued.
    // Hand over the options by copy, so that their shared ptr stays alive in the lambda.
    sub->set_callback( [ opt ]() {
        run_unchunkify( *opt );
    });
}

// =================================================================================================
//      Helpers
// =================================================================================================

/**
 * @brief Check which of the three modes was selected by the user, and return it.
 */
UnchunkifyMode get_unchunkify_mode( UnchunkifyOptions const& options )
{
    UnchunkifyMode mode;
    size_t mode_cnt = 0;

    if( *( options.jplace_input.option() ) && options.jplace_input.file_count() > 0 ) {
        mode = UnchunkifyMode::kJplaceInput;
        ++mode_cnt;

        if( global_options.verbosity() >= 1 ) {
            std::cout << "Selected mode: Jplace Input.\n";
        }
    }
    if( ! options.chunk_list_file.empty() ) {
        mode = UnchunkifyMode::kChunkListFile;
        ++mode_cnt;

        if( global_options.verbosity() >= 1 ) {
            std::cout << "Selected mode: Chunk List File.\n";
        }
    }
    if( ! options.chunk_file_expression.empty() ) {
        mode = UnchunkifyMode::kChunkFileExpression;
        ++mode_cnt;

        if( global_options.verbosity() >= 1 ) {
            std::cout << "Selected mode: Chunk File Expression.\n";
        }
    }
    if( mode_cnt != 1 ) {
        throw CLI::ValidationError(
            "--jplace-path, --chunk-list-file, --chunk-file-expression",
            "Exactly one of the three input modes has to be provided."
        );
    }

    return mode;
}

/**
 * @brief If Jplace Files mode was selected, build the hash map. If not, return an empty map.
 */
template< class HashFunction >
HashToIndexMap<HashFunction> get_hash_to_indices_map(
    UnchunkifyOptions const&  options,
    ChunkCache<HashFunction>& chunk_cache,
    UnchunkifyMode            mode
) {
    HashToIndexMap<HashFunction> hash_map;

    // If we are not in that mode, just return empty.
    if( mode != UnchunkifyMode::kJplaceInput ) {
        return hash_map;
    }

    options.jplace_input.print();

    // Print user output.
    if( global_options.verbosity() >= 2 ) {
        std::cout << "Preparing chunk hash list.\n";
    }

    // Load all (!) chunk files once (possibly removing the earlier ones from the cache while
    // doing so), and for each pquery, store its names and indices for later lookup.
    // As the actual storage is critical, this loop does not scale well.
    // But at least, the file loading is done in parallel... Could optimize if needed.
    #pragma omp parallel for schedule(dynamic)
    for( size_t sample_idx = 0; sample_idx < options.jplace_input.file_count(); ++sample_idx ) {

        auto const file_path = options.jplace_input.file_path( sample_idx );
        auto const chunk = chunk_cache.fetch_copy( file_path );

        for( size_t pquery_idx = 0; pquery_idx < chunk->sample.size(); ++pquery_idx ) {
            auto const& pquery = chunk->sample.at( pquery_idx );

            for( auto const& name : pquery.names() ) {
                auto const digest = HashFunction::hex_to_digest( name.name );
                if( hash_map.count( digest ) > 0 ) {
                    auto const file2 = options.jplace_input.file_path(
                        hash_map[ digest ].sample_index
                    );
                    throw std::runtime_error(
                        "Pquery with hash name '" + name.name + "' exists in multiple files: " +
                        file_path + " and " + file2
                    );
                }

                #pragma omp critical(GAPPA_UNCHUNKIFY_FILL_HASH_INDICES_MAP)
                {
                    hash_map[ digest ] = { sample_idx, pquery_idx };
                }
            }
        }
    }

    // Print user output.
    if( global_options.verbosity() >= 2 ) {
        std::cout << "Prepared chunk hash list.\n";
    }

    return hash_map;
}

/**
 * @brief If Chunk List File mode was selected, read the list file.
 */
std::vector<std::string> get_chunk_list_file(
    UnchunkifyOptions const&  options,
    UnchunkifyMode            mode
) {
    using namespace genesis::utils;
    std::vector<std::string> list;

    // If we are not in that mode, just return empty.
    if( mode != UnchunkifyMode::kChunkListFile ) {
        return list;
    }

    // Read list file.
    auto const cont = split( file_read( options.chunk_list_file ), "\r\n" );
    for( auto const& line : cont ) {
        auto const file_path = trim( line );
        list.push_back( file_path );

        // Fail early.
        if( ! file_exists( file_path )) {
            throw std::runtime_error(
                "In line " + std::to_string( list.size() ) +
                ", chunk file list contains non-existing file: " + file_path
            );
        }
    }

    if( global_options.verbosity() >= 1 ) {
        std::cout << "Found " << list.size() << " jplace files in chunk file list.\n";
    }

    return list;
}

/**
 * @brief Function to load a sample from file, given its file path.
 */
template< class HashFunction >
std::shared_ptr<MappedSample<HashFunction>> load_sample(
    UnchunkifyMode const&           mode,
    ChunkCache<HashFunction> const& chunk_cache,
    std::string const&              file_path
) {
    // Report cache size on every load, that is, whenever the cache actually changes.
    if( global_options.verbosity() >= 3 ) {
        #pragma omp critical(GAPPA_UNCHUNKIFY_PRINT)
        {
            std::cout << "Current jplace cache size: " << chunk_cache.size() << "\n";
        }
    }

    // Create the result and load the sample.
    auto mapped_sample = std::make_shared<MappedSample<HashFunction>>();
    mapped_sample->sample = genesis::placement::JplaceReader().from_file( file_path );

    // If we are in a mode that needs per-sample indicies, create the map from hashes to indices.
    if( mode == UnchunkifyMode::kChunkFileExpression || mode == UnchunkifyMode::kChunkListFile ) {
        for( size_t pquery_idx = 0; pquery_idx < mapped_sample->sample.size(); ++pquery_idx ) {
            auto const& pquery = mapped_sample->sample.at( pquery_idx );

            for( auto const& name : pquery.names() ) {
                auto const digest = HashFunction::hex_to_digest( name.name );
                if( mapped_sample->hash_to_index.count( digest ) > 0 ) {
                    throw std::runtime_error(
                        "Pquery with hash name '" + name.name +
                        "' exists in multiple times in file: " + file_path
                    );
                }

                mapped_sample->hash_to_index[ digest ] = pquery_idx;
            }
        }
    }

    return mapped_sample;
}

/**
 * @brief Sort the abundance by chunk id, in order to minimize loading.
 */
void sort_abundances_by_chunk_id(
    genesis::utils::JsonDocument::ArrayType& arr,
    std::string const& map_filename
) {
    using namespace genesis::utils;

    auto sort_by_chunk_id = [&]( JsonDocument const& lhs, JsonDocument const& rhs )
    {
        // Caution.
        if(
            ! lhs.is_array() || lhs.size() != 3 || ! lhs[1].is_number_unsigned() ||
            ! rhs.is_array() || rhs.size() != 3 || ! rhs[1].is_number_unsigned()
        ) {
            throw std::runtime_error( "Invalid abundance map: " + map_filename );
        }

        return lhs.get_array()[1].get_number_unsigned() < rhs.get_array()[1].get_number_unsigned();
    };
    std::sort( arr.begin(), arr.end(), sort_by_chunk_id );
}

// =================================================================================================
//      Main Work Functions
// =================================================================================================

/**
 * @brief Inner loop function for proessing an abundance map file.
 *
 * Its main input is a sequence entry from the json abundance map.
 * It then uses all kind of other inputs (unfortunately...), and returns the jplace chunk from
 * the cache along with the index of the pquery that was specified in the sequence entry.
 *
 * If no pquery could be found that fits the sequence entry (that is, the hash in there),
 * an empty object is returned, which is checked by the main work function.
 */
template< class HashFunction >
std::pair<std::shared_ptr<MappedSample<HashFunction>>, size_t> get_chunk_and_pquery(
    genesis::utils::JsonDocument const& seq_entry,
    HashToIndexMap<HashFunction> const& hash_to_indices,
    std::vector<std::string>     const& chunk_list,
    UnchunkifyOptions const&            options,
    UnchunkifyMode const&               mode,
    std::string const&                  map_filename,
    ChunkCache<HashFunction>&           chunk_cache
) {
    using namespace genesis::utils;

    // Test the json object:
    //  [0]: the hash hex value
    //  [1]: the chunk number
    //  [2]: the object that maps sequence labels to abundances
    if(
        ! seq_entry.is_array() ||
        seq_entry.size() != 3  ||
        ! seq_entry[0].is_string() ||
        ! seq_entry[1].is_number_unsigned() ||
        ! seq_entry[2].is_object()
    ) {
        throw std::runtime_error( "Invalid abundance map: " + map_filename );
    }

    // Get the hash that we are trying to find the the chunks.
    auto const digest = HashFunction::hex_to_digest( seq_entry[0].get_string() );

    // Get sample path depending on mode.
    std::string sample_file_path;
    size_t pquery_idx = std::numeric_limits<size_t>::max();
    if( mode == UnchunkifyMode::kJplaceInput ) {

        // In jplace input mode, get both sample path and pquery index from the big map.
        if( hash_to_indices.count( digest ) == 0 ) {
            // The hash is not there. Return empty.
            return {};
        }
        auto const& indices = hash_to_indices.at( digest );
        sample_file_path = options.jplace_input.file_path( indices.sample_index );
        pquery_idx = indices.pquery_index;

    } else if( mode == UnchunkifyMode::kChunkFileExpression ) {

        // In expression mode, get the sample path by replacing in the expression.
        auto const chunk_num = seq_entry[1].get_number_unsigned();
        sample_file_path = replace_all(
            options.chunk_file_expression, "@", std::to_string( chunk_num )
        );
        if( ! file_exists( sample_file_path ) ) {
            throw std::runtime_error(
                "Chunk file expression yields non-existing file: " + sample_file_path
            );
        }

    } else if( mode == UnchunkifyMode::kChunkListFile ) {

        // In list file mode, get the path from the list file.
        auto const chunk_num = seq_entry[1].get_number_unsigned();
        if( chunk_num >= chunk_list.size() ) {
            throw std::runtime_error(
                "Chunk index " + std::to_string( chunk_num ) + " does not exist in chunk list file, " +
                "but is required in abundance map file " + map_filename
            );
        }
        sample_file_path = chunk_list[ chunk_num ];

    } else {
        throw std::domain_error( "Invalid unchunkify mode." );
    }

    // Load the chunk.
    auto const chunk = chunk_cache.fetch_copy( sample_file_path );

    // For two modes, we need to get the pquery index from the sample.
    if(
        mode == UnchunkifyMode::kChunkFileExpression ||
        mode == UnchunkifyMode::kChunkListFile
    ) {
        if( chunk->hash_to_index.count( digest ) == 0 ) {
            // The hash is not there. Return empty.
            return {};
        }
        pquery_idx = chunk->hash_to_index.at( digest );
    }

    return { chunk, pquery_idx };
}

/**
 * @brief After a pquery was added to a sample, we need to replace the hash name by
 * the actual sequence labels and abundances fromt he map file.
 */
void add_sequence_names_and_abundances(
    genesis::utils::JsonDocument const& seq_entry,
    genesis::placement::Pquery&         pquery,
    std::string const&                  map_filename
) {
    // That was checked before already, but why not. It's cheap.
    if(
        ! seq_entry.is_array() ||
        seq_entry.size() != 3  ||
        ! seq_entry[0].is_string() ||
        ! seq_entry[1].is_number_unsigned() ||
        ! seq_entry[2].is_object()
    ) {
        throw std::runtime_error( "Invalid abundance map: " + map_filename );
    }

    // Remove the hash name, and add the actual sequence names and abundances.
    pquery.clear_names();
    auto const& mult_arr = seq_entry[2].get_object();
    for( auto const& mult_obj : mult_arr ) {
        auto const& label = mult_obj.first;
        if( ! mult_obj.second.is_number_unsigned() ) {
            throw std::runtime_error( "Invalid abundance map: " + map_filename );
        }
        auto const mult = mult_obj.second.get_number_unsigned();

        pquery.add_name( label, mult );
    }
}

/**
 * @brief Main work function. Loops over all abundance map files and writes a per-sample jplace file
 * for each of them.
 */
template< class HashFunction >
void run_unchunkify_with_hash( UnchunkifyOptions const& options )
{
    using namespace genesis::placement;
    using namespace genesis::utils;

    // Get run mode
    auto const mode = get_unchunkify_mode( options );

    // -----------------------------------------------------------
    //     Prepare Helper Data
    // -----------------------------------------------------------

    // Writer for samples
    auto const jplace_writer = JplaceWriter();

    // Make a cache for storing the jplace chunk files.
    // We load a file given its path. This makes it flexible for the different
    // modes to decide how they get the path from their input.
    ChunkCache<HashFunction> chunk_cache( options.jplace_cache_size );
    chunk_cache.load_function = [&]( std::string const& sample_file_path ){
        return load_sample<HashFunction>( mode, chunk_cache, sample_file_path );
    };

    // Depending on the mode, we need different helper data structures. Prepare them.

    // Mode Jplace Input. We don't have any chunk nums, so instead we just prepare a full
    // lookup from hashes to sample index. This can get big...
    // It is only filled if the mode is actually jplace input.
    auto const hash_to_indices = get_hash_to_indices_map<HashFunction>( options, chunk_cache, mode );

    // Mode Chunk List file. Read the list file.
    // It is only filled if the mode is actually chunk list file.
    auto const chunk_list = get_chunk_list_file( options, mode );

    // -----------------------------------------------------------
    //     Run
    // -----------------------------------------------------------

    // Some statistics for user output.
    size_t file_count = 0;
    size_t total_seqs_count = 0;
    size_t not_found_count = 0;

    // Iterate map files
    #pragma omp parallel for schedule(dynamic)
    for( size_t fi = 0; fi < options.abundance_map_input.file_count(); ++fi ) {
        auto const& map_filename = options.abundance_map_input.file_path( fi );

        // User output
        if( global_options.verbosity() >= 2 ) {
            #pragma omp critical(GAPPA_UNCHUNKIFY_PRINT)
            {
                ++file_count;
                std::cout << "Processing file " << file_count << " of ";
                std::cout << options.abundance_map_input.file_count();
                std::cout << ": " << map_filename << "\n";
            }
        }

        // Read map file and do some checks.
        auto doc = JsonReader().from_file( map_filename );
        if( ! doc.is_object() ) {
            throw std::runtime_error( "Invalid abundance map: " + map_filename );
        }
        auto hash_it = doc.find( "hash" );
        if( hash_it == doc.end() || ! hash_it->is_string() ) {
            throw std::runtime_error( "Invalid abundance map: " + map_filename );
        }
        if( ! equals_ci( hash_it->get_string(), options.hash_function )) {
            throw std::runtime_error(
                "Command was called with hash function " + options.hash_function +
                ", but abundance map file specifies hash function " + hash_it->get_string() +
                ": " + map_filename
            );
        }
        auto abun_it = doc.find( "abundances" );
        if( abun_it == doc.end() || ! abun_it->is_array() ) {
            throw std::runtime_error( "Invalid abundance map: " + map_filename );
        }

        // Sort mapped sequences by chunk id, in order to minimize loading.
        sort_abundances_by_chunk_id( abun_it->get_array(), map_filename );

        // Create empty sample.
        Sample sample;

        // Loop over mapped sequences and add them to the sample.
        for( auto seq_entry_it = abun_it->begin(); seq_entry_it != abun_it->end(); ++seq_entry_it ) {
            #pragma omp atomic
            ++total_seqs_count;

            // Get the chunk and pquery index.
            // If not found, there is no pquery for the current sequence.
            auto const chunk_and_pquery = get_chunk_and_pquery(
                *seq_entry_it, hash_to_indices, chunk_list, options, mode, map_filename, chunk_cache
            );
            auto const chunk = chunk_and_pquery.first;
            auto const pquery_idx = chunk_and_pquery.second;
            if( ! chunk ) {
                #pragma omp atomic
                ++not_found_count;
                continue;
            }

            // New sample: give it a tree!
            if( sample.empty() ) {
                sample = Sample( chunk->sample.tree() );
            }

            // Fill in the sequence, with labels and abundances.
            auto& pquery = sample.add( chunk->sample.at( pquery_idx ));
            add_sequence_names_and_abundances( *seq_entry_it, pquery, map_filename );
        }

        // Get sample name.
        auto sample_name_it = doc.find( "sample" );
        if( sample_name_it == doc.end() || ! sample_name_it->is_string() ) {
            throw std::runtime_error( "Invalid abundance map: " + map_filename );
        }
        auto const& sample_name = sample_name_it->get_string();

        // We are done with the map/sample. Write it.
        jplace_writer.to_file( sample, options.file_output.out_dir() + sample_name + ".jplace" );
    }

    if( global_options.verbosity() >= 1 ) {
        std::cout << "Wrote " << total_seqs_count << " sequences to sample files.\n";
        std::cout << "Could not find " << not_found_count << " sequence hashes.\n";
    }
}

// =================================================================================================
//      Run
// =================================================================================================

void run_unchunkify( UnchunkifyOptions const& options )
{
    using namespace genesis::utils;

    // -----------------------------------------------------------
    //     Input Output File Preparations
    // -----------------------------------------------------------

    // Check if any of the files we are going to produce already exists. If so, fail early.
    options.file_output.check_nonexistent_output_files({ ".*\\.jplace" });

    // Print some user output.
    options.abundance_map_input.print();

    // -----------------------------------------------------------
    //     Run
    // -----------------------------------------------------------

    if( options.hash_function == "SHA1" ) {
        run_unchunkify_with_hash<SHA1>( options );
    } else if( options.hash_function == "SHA256" ) {
        run_unchunkify_with_hash<SHA256>( options );
    } else if( options.hash_function == "MD5" ) {
        run_unchunkify_with_hash<MD5>( options );
    } else {
        throw CLI::ConversionError( "Unknown hash function: " + options.hash_function );
    }
}
