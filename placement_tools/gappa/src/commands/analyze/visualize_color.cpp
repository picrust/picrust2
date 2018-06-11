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

#include "commands/analyze/visualize_color.hpp"

#include "options/global.hpp"

#include "CLI/CLI.hpp"

#include "genesis/placement/formats/jplace_reader.hpp"
#include "genesis/placement/function/helper.hpp"
#include "genesis/placement/function/masses.hpp"
#include "genesis/placement/function/operators.hpp"
#include "genesis/placement/function/tree.hpp"
#include "genesis/tree/default/newick_writer.hpp"
#include "genesis/utils/core/fs.hpp"

#include <algorithm>
#include <stdexcept>

#ifdef GENESIS_OPENMP
#   include <omp.h>
#endif

// =================================================================================================
//      Setup
// =================================================================================================

void setup_visualize_color( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto options = std::make_shared<VisualizeColorOptions>();
    auto sub = app.add_subcommand(
        "visualize-color",
        "Make a tree with edges colored according to the placement mass of the samples."
    );

    // Input files.
    options->jplace_input.add_jplace_input_opt_to_app( sub );
    options->jplace_input.add_mass_norm_opt_to_app( sub, false );
    options->jplace_input.add_point_mass_opt_to_app( sub );
    options->jplace_input.add_ignore_multiplicities_opt_to_app( sub );

    // Color. We allow max, but not min, as this is always 0.
    options->color_map.add_color_list_opt_to_app( sub, "BuPuBk" );
    options->color_norm.add_log_scaling_opt_to_app( sub );
    options->color_norm.add_max_value_opt_to_app( sub );
    options->color_map.add_over_color_opt_to_app( sub );
    options->color_norm.add_min_value_opt_to_app( sub );
    options->color_map.add_under_color_opt_to_app( sub );
    options->color_norm.add_mask_value_opt_to_app( sub );
    options->color_map.add_mask_color_opt_to_app( sub );

    // Output files.
    options->tree_output.add_tree_output_opts_to_app( sub );
    options->file_output.add_output_dir_opt_to_app( sub );
    options->file_output.add_file_prefix_opt_to_app( sub, "tree", "tree" );

    // Set the run function as callback to be called when this subcommand is issued.
    // Hand over the options by copy, so that their shared ptr stays alive in the lambda.
    sub->set_callback( [options]() {
        run_visualize_color( *options );
    });
}

// =================================================================================================
//      Run
// =================================================================================================

void run_visualize_color( VisualizeColorOptions const& options )
{
    using namespace genesis;
    using namespace genesis::placement;
    using namespace genesis::tree;

    // Prepare output file names and check if any of them already exists. If so, fail early.
    std::vector<std::string> files_to_check;
    for( auto const& e : options.tree_output.get_extensions() ) {
        files_to_check.push_back( options.file_output.file_prefix() + "\\." + e );
    }
    options.file_output.check_nonexistent_output_files( files_to_check );

    // User output.
    options.jplace_input.print();

    // Prepare results.
    Tree tree;
    std::vector<double> total_masses;
    size_t file_count = 0;

    // Read all jplace files and accumulate their masses.
    #pragma omp parallel for schedule(dynamic)
    for( size_t fi = 0; fi < options.jplace_input.file_count(); ++fi ) {

        // User output
        if( global_options.verbosity() >= 2 ) {
            #pragma omp critical(GAPPA_VISUALIZE_PRINT_PROGRESS)
            {
                ++file_count;
                std::cout << "Processing file " << file_count << " of " << options.jplace_input.file_count();
                std::cout << ": " << options.jplace_input.file_path( fi ) << "\n";
            }
        }

        // Read in file. This also already applies all normalizations.
        auto const sample = options.jplace_input.sample( fi );

        // Get masses per edge.
        auto const masses = placement_mass_per_edges_with_multiplicities( sample );

        // The main accumulation is single threaded.
        // We could optimize more, but seriously, it is fast enough already.
        #pragma omp critical(GAPPA_VISUALIZE_ACCUMULATE)
        {
            // Tree
            if( tree.empty() ) {
                tree = sample.tree();
            } else if( ! genesis::placement::compatible_trees( tree, sample.tree() ) ) {
                throw std::runtime_error( "Input jplace files have differing reference trees." );
            }

            // Masses
            if( total_masses.empty() ) {
                total_masses = masses;
            } else if( total_masses.size() != masses.size() ) {
                throw std::runtime_error( "Input jplace files have differing reference trees." );
            } else {
                for( size_t i = 0; i < masses.size(); ++i ) {
                    total_masses[i] += masses[i];
                }
            }
        }
    }

    // If we use relative masses, we normalize the whole mass set once more, so that the sum is 1.
    if( options.jplace_input.mass_norm_relative() ) {
        auto const sum = std::accumulate( total_masses.begin(), total_masses.end(), 0.0 );
        for( size_t i = 0; i < total_masses.size(); ++i ) {
            total_masses[i] /= sum;
        }
    }

    // Get color map and norm.
    auto color_map  = options.color_map.color_map();
    auto color_norm = options.color_norm.get_sequential_norm();

    // First, autoscale to get the max.
    // Finally, apply the user settings that might have been provided.
    color_norm->autoscale( total_masses );
    if( options.color_norm.log_scaling() ) {

        // Some user friendly safety. Min of 0 does not work with log scaling.
        // Instead, if we have a max > 1, we set min to 1, which is a good case for absolute abundances.
        // For relative abundances (normalized samples), the max is < 1, so we set the min to some
        // value below that that spans some orders of magnitude. This is all used as default anyway,
        // as users can overwrite this via --min-value.
        if( color_norm->min_value() <= 0.0 ) {
            if( color_norm->max_value() > 1.0 ) {
                color_norm->min_value( 1.0 );
            } else {
                color_norm->min_value( color_norm->max_value() / 10e4 );
            }
            color_map.clip_under( true );

            // Can't really issue a warning here that easily. Users might have set the --min-value options,
            // which will overwrite this. So, we'd need some more complex checks... not now!
            // std::cout << "Warning: Some branches have mass 0, which cannot be shown using --log-scaling. ";
            // std::cout << "Hence, the minimum was set to " << color_norm->min_value() << " instead.\n";
        }

    } else {
        color_norm->min_value( 0.0 );
    }

    // Now overwrite the above "default" settings with what the user specified
    // (in case that they actually did specify certain values).
    options.color_norm.apply_options( *color_norm );

    // Now, make a color vector and write to files.
    auto const colors = color_map( *color_norm, total_masses );
    options.tree_output.write_tree_to_files(
        tree,
        colors,
        color_map,
        *color_norm,
        options.file_output.out_dir() + options.file_output.file_prefix()
    );
}
