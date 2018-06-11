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

#include "commands/analyze/kmeans_phylogenetic.hpp"

#include "commands/analyze/kmeans.hpp"
#include "options/global.hpp"

#include "CLI/CLI.hpp"

#include "genesis/placement/function/epca.hpp"
#include "genesis/placement/function/operators.hpp"
#include "genesis/tree/function/tree_set.hpp"
#include "genesis/tree/mass_tree/emd.hpp"
#include "genesis/tree/mass_tree/functions.hpp"
#include "genesis/tree/mass_tree/kmeans.hpp"
#include "genesis/tree/mass_tree/tree.hpp"
#include "genesis/utils/io/output_stream.hpp"
#include "genesis/utils/text/string.hpp"

#include <fstream>

#ifdef GENESIS_OPENMP
#   include <omp.h>
#endif

// =================================================================================================
//      Setup
// =================================================================================================

void setup_pkmeans( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto opt = std::make_shared<PkmeansOptions>();
    auto sub = app.add_subcommand(
        "phylogenetic-kmeans",
        "Run Phylogenetic k-means clustering on a set of samples."
    );

    // Setup common kmeans options.
    setup_kmeans( opt.get(), sub, "pkmeans_" );

    // Binning.
    auto bins_opt = sub->add_option(
        "--bins",
        opt->bins,
        "Bin the masses per-branch in order to save time and memory. "
        "Default is 0, that is, no binning. If set, we recommend to use 50 bins or more.",
        true
    );
    bins_opt->group( "Settings" );

    // Set the run function as callback to be called when this subcommand is issued.
    // Hand over the options by copy, so that their shared ptr stays alive in the lambda.
    sub->set_callback( [ opt ]() {
        run_pkmeans( *opt );
    });
}

// =================================================================================================
//      Helper Functions
// =================================================================================================

void write_pkmeans_cluster_trees(
    PkmeansOptions const& options,
    std::vector<genesis::tree::MassTree> const& centroids,
    size_t k
) {

    if( centroids.size() != k ) {
        throw std::runtime_error(
            "Internal Error: Differing number of centroids (" + std::to_string( centroids.size() ) +
            ") and  k (" + std::to_string( k ) + ")."
        );
    }

    // Get color map and norm.
    auto color_map  = options.color_map.color_map();
    auto color_norm = options.color_norm.get_sequential_norm();

    // Out base file name
    auto const base_fn = options.file_output.out_dir() + cluster_tree_basepath( options, k );

    // Write all centroid trees
    for( size_t ci = 0; ci < centroids.size(); ++ci ) {
        auto const& centroid = centroids[ci];

        // Prepare colors
        auto const masses = mass_tree_mass_per_edge( centroid );
        color_norm->autoscale_max( masses );

        // Now, make a color vector and write to files.
        auto const colors = color_map( *color_norm, masses );
        auto const cntr_fn = base_fn + std::to_string( ci );
        options.tree_output.write_tree_to_files(
            centroid,
            colors,
            color_map,
            *color_norm,
            cntr_fn
        );
    }
}

// =================================================================================================
//      Main Run Function
// =================================================================================================

void run_pkmeans( PkmeansOptions const& options )
{
    using namespace genesis;
    using namespace genesis::placement;
    using namespace genesis::tree;
    using namespace genesis::utils;

    // Print some user output.
    options.jplace_input.print();

    // Base check
    if( options.jplace_input.file_count() < 2 ) {
        throw std::runtime_error( "Cannot run k-means with fewer than 2 samples." );
    }

    // Get the values of k to run.
    auto const ks = get_k_values( options );

    // Check for existing files.
    check_kmeans_output_files( options );

    if( global_options.verbosity() >= 1 ) {
        std::cout << "Reading samples.\n";
    }

    // Read in the trees and immediately convert them to mass trees to save storage.
    auto mass_trees = options.jplace_input.mass_tree_set();

    // Binning.
    if( options.bins > 0 ) {
        for( auto& mt : mass_trees ) {
            mass_tree_binify_masses( mt, options.bins );
        }
    }

    // Set up kmeans.
    auto mkmeans = MassTreeKmeans();
    mkmeans.report_iteration = [&]( size_t iteration ){
        if( global_options.verbosity() >= 2 ) {
            std::cout << " - Iteration " << iteration << "\n";
        }
    };
    if( options.bins > 0 ) {
        mkmeans.accumulate_centroid_masses( options.bins );
    }

    // Run kmeans for every specified k.
    std::vector<KmeansClusterOverview> overview;
    for( auto const& k : ks ) {

        // Run it.
        std::cout << "\nRunning Phylogenetic Kmeans with k=" << k << "\n";
        auto const iterations = mkmeans.run( mass_trees, k );
        auto const clust_info = mkmeans.cluster_info( mass_trees );
        std::cout << "Finished after " << iterations << " iterations\n";

        // Write output.
        write_assignment_file( options, mkmeans.assignments(), clust_info, k );
        write_pkmeans_cluster_trees( options, mkmeans.centroids(), k );

        // Print some cluster info, and collect it for the overview file.
        auto const ci = print_cluster_info( options, mkmeans.assignments(), clust_info, k );
        overview.push_back( ci );
    }

    // Write the overview file for elbow plots etc.
    write_overview_file( options, overview );
}
