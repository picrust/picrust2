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

#include "commands/prepare/extract.hpp"

#include "options/global.hpp"

#include "CLI/CLI.hpp"

#include "genesis/placement/formats/jplace_reader.hpp"
#include "genesis/placement/formats/jplace_writer.hpp"
#include "genesis/placement/function/functions.hpp"
#include "genesis/placement/function/helper.hpp"
#include "genesis/placement/function/operators.hpp"
#include "genesis/placement/function/sample_set.hpp"
#include "genesis/placement/function/tree.hpp"

#include "genesis/sequence/sequence.hpp"
#include "genesis/sequence/sequence_set.hpp"
#include "genesis/sequence/formats/fasta_input_iterator.hpp"
#include "genesis/sequence/formats/fasta_writer.hpp"
#include "genesis/sequence/functions/labels.hpp"

#include "genesis/tree/tree.hpp"
#include "genesis/tree/bipartition/bipartition.hpp"
#include "genesis/tree/bipartition/functions.hpp"
#include "genesis/tree/default/functions.hpp"
#include "genesis/tree/drawing/functions.hpp"
#include "genesis/tree/drawing/circular_layout.hpp"
#include "genesis/tree/drawing/layout_tree.hpp"
#include "genesis/utils/io/output_stream.hpp"

#include "genesis/utils/core/fs.hpp"
#include "genesis/utils/formats/csv/reader.hpp"
#include "genesis/utils/text/string.hpp"
#include "genesis/utils/tools/color/qualitative_lists.hpp"

#include <algorithm>
#include <cassert>
#include <string>
#include <stdexcept>
#include <vector>

#ifdef GENESIS_OPENMP
#   include <omp.h>
#endif

// =================================================================================================
//     Typedefs
// =================================================================================================

/**
 * @brief Contains a list of clades, each itself containing a list of taxa belonging to that clade.
 */
using CladeTaxaList = std::unordered_map<std::string, std::vector<std::string>>;

/**
 * @brief Contains a list of clades, each itself containing a list of edge indices belonging to
 * that clade. We use a vector to maintain the order of the clades.
 * This is re-computed again for every input file, to make sure that the indices are correct.
 */
using CladeEdgeList = std::vector<std::pair<std::string, std::unordered_set<size_t>>>;

/**
 * @brief Contains a list of clades, each itself contianing all pquery names in that clade.
 */
using PqueryNamesPerCladeList = std::unordered_map<std::string, std::unordered_set<std::string>>;

// =================================================================================================
//      Setup
// =================================================================================================

void setup_extract( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto options = std::make_shared<ExtractOptions>();
    auto sub = app.add_subcommand(
        "extract",
        "Extract placements from clades of the tree and write per-clade jplace files."
    );

    // Jplace input.
    options->jplace_input.add_jplace_input_opt_to_app( sub );

    // Clade list file
    auto clade_list_file_opt = sub->add_option(
        "--clade-list-file",
        options->clade_list_file,
        "File containing a tab-separated list of taxon to clade mapping."
    );
    clade_list_file_opt->required();
    clade_list_file_opt->check( CLI::ExistingFile );
    clade_list_file_opt->group( "Input" );

    // Other Input files.
    options->sequence_input.add_fasta_input_opt_to_app( sub, false );

    // Other options
    auto threshold_opt = sub->add_option(
        "--threshold",
        options->threshold,
        "Threshold of how much placement mass needs to be in a clade for extracting a pquery.",
        true
    );
    threshold_opt->check( CLI::Range( 0.5, 1.0 ));
    threshold_opt->group( "Settings" );

    // Make this a settings option.
    options->jplace_input.add_point_mass_opt_to_app( sub );

    auto color_tree_file_opt = sub->add_option(
        "--color-tree-file",
        options->color_tree_file,
        "If a path is provided, an svg file with a tree colored by clades is written."
    );
    color_tree_file_opt->check( CLI::NonexistentPath );
    color_tree_file_opt->group( "Output" );

    // Output files.
    options->jplace_output.add_output_dir_opt_to_app( sub, "samples", "samples" );
    // options->jplace_output.add_file_prefix_opt_to_app( sub, "sample", "" );
    options->sequence_output.add_output_dir_opt_to_app( sub, "sequences", "sequences" );
    // options->sequence_output.add_file_prefix_opt_to_app( sub, "sequence", "" );

    // Set the run function as callback to be called when this subcommand is issued.
    // Hand over the options by copy, so that their shared ptr stays alive in the lambda.
    sub->set_callback( [options]() {
        run_extract( *options );
    });
}

// =================================================================================================
//     Get Clade Taxa Lists
// =================================================================================================

/**
 * @brief Return a list of clades, each containing a list of taxa.
 *
 * The function takes the clade file path from the options as input. Each line of that file
 * contains a tab-separated entry that maps from a taxon of the tree to the clade name that this
 * taxon belongs to:
 *
 *     Taxon_1, Clade_a
 *
 * White spaces around the taxa and clades are stripped.
 * The return value of this function is a map from clade names to a vector of taxa names.
 */
CladeTaxaList get_clade_taxa_lists( ExtractOptions const& options )
{
    using namespace ::genesis::utils;

    auto const & clade_filename = options.clade_list_file;
    auto csv_reader = CsvReader();
    csv_reader.separator_chars( "\t" );

    // Create a list of all clades and fill each clade with its taxa.
    CladeTaxaList clades;
    auto table = csv_reader.from_file( clade_filename );
    for( size_t i = 0; i < table.size(); ++i ) {
        auto const& line = table[i];
        if( line.size() != 2 ) {
            throw std::runtime_error(
                "Expecting two tab-separated fields in clade file. "
                "Invalid line " + std::to_string(i) + " in clade file: " + clade_filename
            );
        }

        auto const taxon = trim( line[0] );
        auto const clade = trim( line[1] );

        // Overly pedentic check to protect from mistakes:
        // Check if the taxon was used before.
        for( auto const& cl : clades ) {
            for( auto const& tx : cl.second ) {
                if( tx == taxon ) {
                    throw std::runtime_error(
                        "Taxon " + taxon + " occurs multiple times in clade file: " + clade_filename
                    );
                }
            }
        }

        // Add the taxon to its clade.
        clades[ clade ].push_back( taxon );
    }

    // We will use two speciled clades later. Check here that they are not in there yet.
    if( clades.count( options.basal_clade_name ) > 0 ) {
        throw std::runtime_error(
            "Clade file contains reserved clade name \"" + options.basal_clade_name + "\": " + clade_filename
        );
    }
    if( clades.count( options.uncertain_clade_name ) > 0 ) {
        throw std::runtime_error(
            "Clade file contains reserved clade name \"" + options.uncertain_clade_name + "\": " + clade_filename
        );
    }

    // Some user output.
    if( global_options.verbosity() >= 1 ) {
        std::cout << "Clade list file contains " << table.size();
        std::cout << " taxa in " << clades.size() << " clades.\n";
    }

    return clades;
}

// =================================================================================================
//     Get Clade Edges
// =================================================================================================

/**
 * @brief Return a list clades, each itself containing a list of edge indices of that clade.
 *
 * The function takes a list of clades with their taxa as input, and a reference tree.
 * It then inspects all clades and findes the edges of the tree that belong to a clade.
 * Furthermore, a clade "basal_branches" is added for those edges of the tree that do not
 * belong to any clade.
 */
CladeEdgeList get_clade_edges(
    ExtractOptions const&      options,
    CladeTaxaList const&       clades,
    genesis::tree::Tree const& tree,
    std::string const&         sample_name
) {
    using namespace ::genesis::tree;

    // Prepare the result map.
    CladeEdgeList clade_edges;

    // Make a set of all edges that do not belong to any clade (the basal branches of the tree).
    // We first fill it with all edge indices, then remove the clade-edges later,
    // so that only the wanted ones remain. This also serves as a check that we do not use any
    // edges in more than one clade.
    std::unordered_set<size_t> basal_branches;
    for( auto it = tree.begin_edges(); it != tree.end_edges(); ++it ) {
        basal_branches.insert( (*it)->index() );
    }

    // Prepare.
    auto const bipartitions = bipartition_set( tree );

    // Process all clades.
    for( auto const& clade : clades ) {

        // Find the nodes that belong to the taxa of this clade.
        std::vector< TreeNode const* > node_list;
        for( auto const& taxon : clade.second ) {
            auto node = find_node( tree, taxon );
            if( node == nullptr ) {
                std::cout << "Cannot find taxon " << taxon << " in tree of sample " << sample_name << "\n";
                continue;
            }
            node_list.push_back( node );
        }

        // Find the edges that are part of the monophyletic subtrees of this clade.
        auto const subedges = find_monophyletic_subtree_edges( tree, bipartitions, node_list );

        // Add to the list.
        // For now, we convert to an unordered map here by hand.
        // This is a bit messy and could be cleaned up in the future.
        clade_edges.push_back( std::make_pair( clade.first, std::unordered_set<size_t>(
            subedges.begin(), subedges.end()
        )));

        // Remove the edge indices of this clade from the basal branches (non-clade) edges list.
        for( auto const edge : subedges ) {

            // Test whether the edge was already removed. If so, clades overlap.
            // Actually, I'm pretty sure that this cannot happen, because of the way that
            // find_monophyletic_subtree_edges() works: If an edge was in two clades, they would
            // not be monophyletic. Thus, this edge is never used by the function...
            // Anyway, we check it to be sure. Not expensive.
            if( basal_branches.count( edge ) == 0 ) {

                // Find the clade that already has that edge, and use it for the error message.
                std::string colliding_clade;
                for( auto const& coll_clade : clade_edges ) {
                    if( coll_clade.second.count( edge ) > 0 ) {
                        colliding_clade = coll_clade.first;
                        break;
                    }
                }
                throw std::runtime_error(
                    "Clades " + clade.first + " and " + colliding_clade + " collidate! "
                    "That is, there are edges that are in both clades. "
                    "The clade file needs to be cleaned up to fix this."
                );
            }

            // Now, remove the edge from basal branches, because it is in a clade now.
            basal_branches.erase( edge );
        }
    }

    // Now that we have processed all clades, also add the non-clade edges (basal branches)
    // to the list as a special clade "basal_branches". This way, all edges of the reference tree
    // are used by exaclty one clade.
    clade_edges.push_back( std::make_pair( options.basal_clade_name, std::move( basal_branches )));

    return clade_edges;
}

// =================================================================================================
//     Extract Pqueries
// =================================================================================================

/**
 * @brief Take a list of edges per clade and a Sample and fill a SampleSet with single samples
 * for all given clades, where each sample contains those pqueries that fell into the clade.
 *
 * This is the main extraaction method. The SampleSet also contains an additional Sample
 * "uncertain", where all pqueries of the provided sample end up which do not have more than
 * `threshold` of their placement mass in a certain clade.
 */
void extract_pqueries(
    ExtractOptions const&             options,
    CladeEdgeList const&              clade_edges,
    genesis::placement::Sample const& sample,
    genesis::placement::SampleSet&    sample_set
) {
    using namespace ::genesis::placement;

    // Process all pqueries of the given sample.
    // We again use openmp here, so that even if only a single file is given, we make use of threads.
    #pragma omp parallel for
    for( size_t qi = 0; qi < sample.size(); ++qi ) {
        auto const& pquery = sample.at( qi );

        // Prepare an accumulator that collects the mass per clade for this pquery.
        // The indices in the vector are the same as the ones in the clade_edge vector.
        std::vector<double> mass_per_clade( clade_edges.size(), 0.0 );

        // For each placement, find its edge and the clade that this edge belongs to.
        // For this clade, accumulate the placement's like weight ratio.
        for( auto const& placement : pquery.placements() ) {

            // Find the correct clade for the edge of this placement by scanning them all
            // (this is a bit inefficient, but for now, it works).
            for( size_t i = 0; i < clade_edges.size(); ++i ) {
                if( clade_edges[i].second.count( placement.edge().index() ) > 0 ) {

                    // If we found the correct clade, accumulate to its mass.
                    mass_per_clade[i] += placement.like_weight_ratio;
                }
            }
        }

        // Now check whether there is a clade that has equal or more than `threshold` percent
        // of the placement's weight ratio. If so, this is the one we assign the pquery to.
        bool found_clade = false;
        for( size_t i = 0; i < mass_per_clade.size(); ++i ) {
            // Only execute the main body of this loop when we found a fitting clade.
            if( mass_per_clade[i] < options.threshold ) {
                continue;
            }
            found_clade = true;

            // Get the sample from the set that has the current clade's name.
            auto sample_ptr = find_sample( sample_set, clade_edges[i].first );
            if( sample_ptr == nullptr ) {
                throw std::runtime_error( "Internal error: Lost sample " + clade_edges[i].first );
            }

            // Add a copy of the pquery to the sample.
            #pragma omp critical(GAPPA_EXTRACT_ADD_SAMPLE)
            {
                sample_ptr->add( pquery );
            }
        }

        // If there is no sure assignment ( < threshold ) for this pquery, we copy it into the
        // special `uncertain` sample.
        if( ! found_clade ) {
            auto sample_ptr = find_sample( sample_set, options.uncertain_clade_name );
            if( sample_ptr == nullptr ) {
                throw std::runtime_error( "Internal error: Lost sample " + options.uncertain_clade_name );
            }

            #pragma omp critical(GAPPA_EXTRACT_ADD_SAMPLE)
            {
                sample_ptr->add( pquery );
            }
        }
    }
}

// =================================================================================================
//     Write Sample Set
// =================================================================================================

/**
 * @brief Take a SampleSet and a directory and write all Samples in the set to jplace files
 * in that directory, named after the sample name in the set.
 */
void write_sample_set(
    genesis::placement::SampleSet const& sample_set,
    ExtractOptions const&                options
) {
    using namespace ::genesis::placement;

    // User output.
    if( global_options.verbosity() >= 1 ) {
        for( size_t si = 0; si < sample_set.size(); ++si ) {
            auto const& named_sample = sample_set.at( si );

            std::cout << "Collected " << named_sample.sample.size() << " pqueries in clade ";
            std::cout << named_sample.name << "\n";
        }

        std::cout << "Writing " << sample_set.size() << " clade sample files.\n";
    }

    // Write files.
    auto writer = JplaceWriter();
    #pragma omp parallel for schedule(dynamic)
    for( size_t si = 0; si < sample_set.size(); ++si ) {
        auto const& named_sample = sample_set.at( si );

        auto const fn = options.jplace_output.file_prefix() + named_sample.name + ".jplace";
        writer.to_file( named_sample.sample, options.jplace_output.out_dir() + fn );
    }
}

// =================================================================================================
//      Write Color Tree
// =================================================================================================

void write_color_tree(
    ExtractOptions const& options,
    CladeEdgeList const&  clade_edges,
    genesis::tree::Tree const& tree
) {
    using namespace ::genesis;
    using namespace ::genesis::tree;
    using namespace ::genesis::utils;

    // Only write a tree file if user specified a path.
    if( options.color_tree_file.empty() ) {
        return;
    }

    // Prepare. We use a grey base color for basal branches.
    auto const base_color = Color( 0.7, 0.7, 0.7 );
    std::vector<utils::Color> color_vector( tree.edge_count(), base_color );
    std::vector<std::string> names;

    // Make a long list of colors.
    auto colors = color_list_set1();
    colors.insert( colors.end(), color_list_dark2().begin(), color_list_dark2().end() );
    colors.insert( colors.end(), color_list_paired().begin(), color_list_paired().end() );
    auto color_map = ColorMap( colors );

    // Colorize the edges and collect the names.
    for( size_t ci = 0; ci < clade_edges.size(); ++ci ) {

        // Do not colorize basal branches.
        if( clade_edges[ci].first == options.basal_clade_name ) {
            continue;
        }

        names.push_back( clade_edges[ci].first );
        for( auto ei : clade_edges[ci].second ) {
            if( color_vector[ ei ] != base_color ) {
                throw std::runtime_error(
                    "Internal error: Overlapping branches!"
                );
            }
            color_vector[ ei ] = color_map.color( ci );
        }
    }

    // Make a layout tree.
    LayoutParameters params;
    params.stroke.width = 8.0;
    auto layout = CircularLayout( tree, params.type, params.ladderize );

    // Set edge colors.
    std::vector<SvgStroke> strokes;
    for( auto const& color : color_vector ) {
        auto stroke = params.stroke;
        stroke.color = color;
        stroke.line_cap = SvgStroke::LineCap::kRound;
        strokes.push_back( std::move( stroke ));
    }
    layout.set_edge_strokes( strokes );

    // Prepare svg doc.
    auto svg_doc = layout.to_svg_document();

    // Add color list
    auto svg_color_list = make_svg_color_list( color_map, names );
    svg_color_list.transform.append( SvgTransform::Translate(
        svg_doc.bounding_box().width() / 2.0 + 200, 0.0
    ));
    svg_color_list.transform.append( SvgTransform::Scale( 5.0, 5.0 ) );
    svg_doc << svg_color_list;

    // Write to file.
    svg_doc.margin = SvgMargin( 200.0 );
    svg_doc.margin.right += 600;
    std::ofstream ofs;
    std::string const ext = ends_with( options.color_tree_file, ".svg" ) ? "" : ".svg";
    auto const fn = options.color_tree_file + ext;
    dir_create( file_path( fn ), true );
    utils::file_output_stream( fn, ofs );
    svg_doc.write( ofs );
}

// =================================================================================================
//      Clade Names List
// =================================================================================================

/**
 * @brief Get the names of all pqueries per clade.
 */
PqueryNamesPerCladeList get_pqueries_per_clade(
    genesis::placement::SampleSet const& sample_set
) {
    PqueryNamesPerCladeList list;
    size_t duplicate_names = 0;

    for( auto const& named_sample : sample_set ) {
        for( auto const& pquery : named_sample.sample ) {
            for( auto const& pquery_name : pquery.names() ) {

                // Check if it is not already there.
                for( auto const& clade : list ) {
                    if( clade.second.count( pquery_name.name ) > 0 ) {
                        ++duplicate_names;
                    }
                }

                // Add it to the list.
                list[ named_sample.name ].insert( pquery_name.name );
            }
        }
    }

    if( duplicate_names > 0 ) {
        std::cout << "Found " << duplicate_names << " pqueries that have the same name. ";
        std::cout << "This will cause the extraction of sequences with that name to be ";
        std::cout << "randomly assigned to one of the clades that have a pquery with that name. ";
        std::cout << "Thus, this should better be fixed first!\n";
    }

    return list;
}

// =================================================================================================
//      Clade Names List
// =================================================================================================

void extract_sequences(
    ExtractOptions const& options,
    PqueryNamesPerCladeList const& list
) {
    using namespace ::genesis;
    using namespace ::genesis::sequence;
    using namespace ::genesis::utils;

    // User output.
    options.sequence_input.print();

    // Helper: Given a clade name, get the fasta file to write to.
    auto clade_filename = [&]( std::string const& cladename ){
        auto path = options.sequence_output.out_dir();
        path += options.sequence_output.file_prefix() + cladename + ".fasta";
        return path;
    };

    // Lazy prep: We write empty files first, which makes sure (again) that they do not exist,
    // and that we can later append to them.
    for( auto const& clade : list ) {
        file_write( "", clade_filename( clade.first ));
    }

    // Helpers.
    auto const set_size = options.sequence_input.file_count();
    size_t file_count = 0;
    auto writer = FastaWriter();

    // Count for user output.
    size_t total_seqs_count = 0;
    size_t missing_seqs_count = 0;

    // Cache sequences per clade while extracting.
    std::unordered_map< std::string, SequenceSet > clade_sequences;

    // Process the input files.
    #pragma omp parallel for schedule(dynamic)
    for( size_t fi = 0; fi < set_size; ++fi ) {
        auto const& fasta_filename = options.sequence_input.file_path( fi );

        // User output.
        if( global_options.verbosity() >= 2 ) {
            #pragma omp critical(GAPPA_EXTRACT_PRINT_PROGRESS)
            {
                ++file_count;
                std::cout << "Processing file " << file_count << " of " << set_size;
                std::cout << ": " << options.sequence_input.file_path( fi ) << "\n";
            }
        }

        auto it = FastaInputIterator( options.sequence_input.fasta_reader() );
        for( it.from_file( fasta_filename ); it; ++it ) {
            #pragma omp atomic
            ++total_seqs_count;

            // Try to find the fasta sequence name in a clade.
            PqueryNamesPerCladeList::const_iterator clade_it = list.end();
            for( auto lit = list.begin(); lit != list.end(); ++lit ) {
                if( lit->second.count( it->label() ) > 0 ) {
                    clade_it = lit;
                    break;
                }
            }

            // If not, skip this sequence.
            if( clade_it == list.end() ) {
                #pragma omp atomic
                ++missing_seqs_count;
                continue;
            }

            // Add the sequence to the clade set, and write if cache too big.
            #pragma omp critical(GAPPA_CHUNKIFY_UPDATE_SEQUENCES)
            {
                // Get the sequence cache for the clade and add ther sequence to it.
                auto& clade_seqs = clade_sequences[ clade_it->first ];
                clade_seqs.add( *it );

                // If the cache is full, flush it.
                if( clade_seqs.size() >= 1000 ) {

                    // Prepare file for appending.
                    auto const fn = clade_filename( clade_it->first );
                    std::ofstream out_stream( fn, std::ofstream::app );
                    if( out_stream.fail() ) {
                        throw std::runtime_error( "Cannot append sequences to file " + fn + "." );
                    }

                    // Write and clear.
                    writer.to_stream( clade_seqs, out_stream );
                    clade_seqs.clear();
                }
            }
        }
    }

    // Flush clade caches.
    for( auto const& clade_seqs : clade_sequences ) {

        // Prepare file for appending.
        auto const fn = clade_filename( clade_seqs.first );
        std::ofstream out_stream( fn, std::ofstream::app );
        if( out_stream.fail() ) {
            throw std::runtime_error( "Cannot append sequences to file '" + fn + "'." );
        }

        // Write.
        writer.to_stream( clade_seqs.second, out_stream );
    }

    std::cout << "Collected " << total_seqs_count << " sequences in " << list.size() << " clades.\n";
    if( missing_seqs_count > 0 ) {
        std::cout << "Thereof, " << missing_seqs_count << " sequences could not be assigned to any ";
        std::cout << "clade, because their name does not appear in any jplace file.\n";
    }
}

// =================================================================================================
//      Run
// =================================================================================================

void run_extract( ExtractOptions const& options )
{
    using namespace ::genesis;
    using namespace ::genesis::placement;
    using namespace ::genesis::tree;

    // Prepare output file names and check if any of them already exists. If so, fail early.
    options.jplace_output.check_nonexistent_output_files({
        options.jplace_output.file_prefix() + ".*\\.jplace"
    });
    if( options.sequence_input.file_count() > 0 ) {
        options.sequence_output.check_nonexistent_output_files({
            options.sequence_output.file_prefix() + ".*\\.fasta"
        });
    }

    // User output.
    options.jplace_input.print();

    // Preparations.
    size_t file_count = 0;
    auto const set_size = options.jplace_input.file_count();

    // Clade taxa list.
    auto const clade_taxa_list = get_clade_taxa_lists( options );

    // We store one tree for the colour output and for checking that all samples have the same one.
    Tree tree;

    // Resulting sample set, gets filled with the extracted pqueries for each clade.
    SampleSet sample_set;

    #pragma omp parallel for schedule(dynamic)
    for( size_t fi = 0; fi < set_size; ++fi ) {

        // User output.
        if( global_options.verbosity() >= 2 ) {
            #pragma omp critical(GAPPA_EXTRACT_PRINT_PROGRESS)
            {
                ++file_count;
                std::cout << "Processing file " << file_count << " of " << set_size;
                std::cout << ": " << options.jplace_input.file_path( fi ) << "\n";
            }
        }

        // Read the sample.
        auto sample = options.jplace_input.sample( fi );
        auto const fn = options.jplace_input.base_file_name( fi );
        auto const clade_edges = get_clade_edges( options, clade_taxa_list, sample.tree(), fn );

        // Prepare tree and sample set. Add samples for every clade that we are going to use.
        // Slighly too conservative in terms of locking, but okay for now.
        #pragma omp critical(GAPPA_EXTRACT_REF_TREE_SAMPLE_SET)
        {
            if( tree.empty() ) {
                tree = sample.tree();

                // Write a tree with clade colors for error checking.
                write_color_tree( options, clade_edges, tree );

                for( auto const& cl : clade_taxa_list ) {
                    sample_set.add( Sample( tree ), cl.first );
                }
                sample_set.add( Sample( tree ), options.basal_clade_name );
                sample_set.add( Sample( tree ), options.uncertain_clade_name );

            } else if( ! genesis::placement::compatible_trees( tree, sample.tree() )) {
                throw std::runtime_error( "Input jplace files have differing reference trees." );
            }
        }

        // Normalize the like_weight_ratios. This step makes sure that missing placement weights do not
        // lead to a pquery being placed in the uncertain clade. That means, we only use the provided
        // placement masses as given in the jplace files, and scale them so that they sum up to 1.0.
        // In turn, this means that uncertainties resulting from the placement algorithm are ignored.
        normalize_weight_ratios( sample );

        // Do the work!
        extract_pqueries( options, clade_edges, sample, sample_set );
    }

    // Write everything to jplace files.
    write_sample_set( sample_set, options );

    // If there were sequences given as input as well, extract them!
    // We can also delete the samples to save some mem. Not needed any more.
    if( options.sequence_input.file_count() > 0 ) {
        auto const list = get_pqueries_per_clade( sample_set );
        sample_set.clear();
        extract_sequences( options, list );
    }
}
