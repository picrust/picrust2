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

#include "commands/analyze/krd.hpp"

#include "options/global.hpp"

#include "CLI/CLI.hpp"

#include "genesis/placement/function/functions.hpp"
#include "genesis/placement/function/emd.hpp"
#include "genesis/placement/function/operators.hpp"

#include "genesis/tree/default/distances.hpp"
#include "genesis/tree/function/distances.hpp"
#include "genesis/tree/function/functions.hpp"
#include "genesis/tree/mass_tree/emd.hpp"
#include "genesis/tree/mass_tree/functions.hpp"
#include "genesis/tree/mass_tree/tree.hpp"

#include "genesis/utils/containers/matrix.hpp"
#include "genesis/utils/containers/matrix/operators.hpp"
#include "genesis/utils/io/output_stream.hpp"

#include <fstream>

#ifdef GENESIS_OPENMP
#   include <omp.h>
#endif

// =================================================================================================
//      Setup
// =================================================================================================

void setup_krd( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto opt = std::make_shared<KrdOptions>();
    auto sub = app.add_subcommand(
        "krd",
        "Calcualte the pairwise Kantorovich-Rubinstein (KR) Distance between samples."
    );

    // Add common options.
    opt->jplace_input.add_jplace_input_opt_to_app( sub );
    opt->jplace_input.add_point_mass_opt_to_app( sub );
    opt->matrix_output.add_matrix_output_opts_to_app( sub, "distance" );

    // Add custom options.
    sub->add_option(
        "--exponent",
        opt->exponent,
        "Exponent for KR integration.",
        true
    );

    // Set the run function as callback to be called when this subcommand is issued.
    // Hand over the options by copy, so that their shared ptr stays alive in the lambda.
    sub->set_callback( [ opt ]() {
        run_krd( *opt );
    });
}

// =================================================================================================
//      Run
// =================================================================================================

void run_krd( KrdOptions const& options )
{
    using namespace genesis;
    using namespace genesis::placement;
    using namespace genesis::tree;
    using namespace genesis::utils;

    // Check if any of the files we are going to produce already exists. If so, fail early.
    // TODO this is ugly.
    options.matrix_output.check_nonexistent_output_files({ options.matrix_output.output_filename() });

    // Print some user output.
    options.jplace_input.print();
    if( global_options.verbosity() >= 1 ) {
        std::cout << "Reading samples.\n";
    }

    // Prepare storage.
    auto const set_size = options.jplace_input.file_count();
    auto mass_trees = std::vector<MassTree>( set_size );
    size_t file_count = 0;

    // TODO branch length and compatibility checks!

    // Load files.
    #pragma omp parallel for schedule(dynamic)
    for( size_t fi = 0; fi < set_size; ++fi ) {

        // User output.
        if( global_options.verbosity() >= 2 ) {
            #pragma omp critical(GAPPA_NHD_PRINT_PROGRESS)
            {
                ++file_count;
                std::cout << "Processing file " << file_count << " of " << set_size;
                std::cout << ": " << options.jplace_input.file_path( fi ) << "\n";
            }
        }

        // Read in file.
        auto const sample = options.jplace_input.sample( fi );

        // Turn it into a mass tree.
        mass_trees[fi] = convert_sample_to_mass_tree( sample ).first;
    }

    if( global_options.verbosity() >= 1 ) {
        std::cout << "Calculating pairwise KR distances.\n";
    }

    // Calculate result matrix.
    auto const krd_matrix = earth_movers_distance( mass_trees, options.exponent );

    if( global_options.verbosity() >= 1 ) {
        std::cout << "Writing distance matrix.\n";
    }
    options.matrix_output.write_matrix( krd_matrix );

    if( global_options.verbosity() >= 1 ) {
        std::cout << "Finished.\n";
    }
}
