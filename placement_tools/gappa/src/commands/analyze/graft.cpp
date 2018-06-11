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

#include "commands/analyze/graft.hpp"

#include "options/global.hpp"

#include "CLI/CLI.hpp"

#include "genesis/placement/formats/jplace_reader.hpp"
#include "genesis/placement/function/tree.hpp"
#include "genesis/tree/default/newick_writer.hpp"
#include "genesis/utils/core/fs.hpp"

#ifdef GENESIS_OPENMP
#   include <omp.h>
#endif

// =================================================================================================
//      Setup
// =================================================================================================

void setup_graft( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto opt = std::make_shared<GraftOptions>();
    auto sub = app.add_subcommand(
        "graft",
        "Make a tree with each of the query sequences represented as a pendant edge."
    );

    // Add input options.
    opt->jplace_input.add_jplace_input_opt_to_app( sub );

    // Fill in custom options.
    sub->add_flag(
        "--fully-resolve", opt->fully_resolve,
        "If set, branches that contain multiple pqueries are resolved by creating a new branch "
        "for each of the pqueries individually, placed according to their distal/proximal lengths. "
        "If not set (default), all pqueries at one branch are collected in a subtree "
        "that branches off from the branch."
    )->group( "Settings" );
    sub->add_option(
        "--name-prefix", opt->name_prefix,
        "Specify a prefix to be added to all new leaf nodes, i.e., to the query sequence names.",
        true
    )->group( "Settings" );

    // Add output options.
    opt->file_output.add_output_dir_opt_to_app( sub );

    // Set the run function as callback to be called when this subcommand is issued.
    // Hand over the options by copy, so that their shared ptr stays alive in the lambda.
    sub->set_callback( [opt]() {
        run_graft( *opt );
    });
}

// =================================================================================================
//      Run
// =================================================================================================

void run_graft( GraftOptions const& options )
{
    using namespace genesis;
    using namespace genesis::placement;

    // Prepare output file names and check if any of them already exists. If so, fail early.
    std::vector<std::string> out_tree_files;
    for( auto const& bfn : options.jplace_input.base_file_names() ) {
        out_tree_files.push_back( bfn + ".newick" );
    }
    options.file_output.check_nonexistent_output_files( out_tree_files );

    // Print some user output.
    options.jplace_input.print();

    // TODO add support for external trees, e.g., bootstrap trees.
    // for this, make sure that the attribute tree is used so that all values can be captured.
    // that probably means there has to be an option to specify which values go where
    // (edges or nodes). or, if only newick is used, not, because we do not reroot,
    // so the assignment does not change.

    size_t fc = 0;
    #pragma omp parallel for schedule(dynamic)
    for( size_t i = 0; i < options.jplace_input.file_count(); ++i ) {

        // User output.
        if( global_options.verbosity() >= 2 ) {
            #pragma omp critical(GAPPA_JPLACE_INPUT_PROGRESS)
            {
                ++fc;
                std::cout << "Reading file " << fc << " of " << options.jplace_input.file_count();
                std::cout << ": " << options.jplace_input.file_path( i ) << "\n";
            }
        }

        // Read the sample and make the tree.
        auto const sample = options.jplace_input.sample( i );
        auto const tog    = labelled_tree( sample, options.fully_resolve, options.name_prefix );

        // Write output to file.
        tree::DefaultTreeNewickWriter().to_file( tog, options.file_output.out_dir() + out_tree_files[i] );
    }
}
