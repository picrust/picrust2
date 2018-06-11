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

#include "commands/analyze/kmeans_imbalance.hpp"

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

void setup_ikmeans( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto opt = std::make_shared<IkmeansOptions>();
    auto sub = app.add_subcommand(
        "imbalance-kmeans",
        "Run Imbalance k-means clustering on a set of samples."
    );

    // Setup common kmeans options.
    setup_kmeans( opt.get(), sub, "ikmeans_" );

    // Set the run function as callback to be called when this subcommand is issued.
    // Hand over the options by copy, so that their shared ptr stays alive in the lambda.
    sub->set_callback( [ opt ]() {
        run_ikmeans( *opt );
    });
}

// =================================================================================================
//      Helper Functions
// =================================================================================================

void write_ikmeans_cluster_trees(
    IkmeansOptions const& options,
    JplaceInputOptions::PlacementProfile const& profile,
    std::vector<size_t> const& assignments,
    size_t k
) {
    using namespace genesis;
    using namespace genesis::utils;

    // Assertions
    if( assignments.size() != profile.edge_masses.rows() ) {
        throw std::runtime_error(
            "Internal Error: Differing number of assignments (" + std::to_string( assignments.size() ) +
            ") and edge masses (" + std::to_string( profile.edge_masses.rows() ) + ")."
        );
    }

    // Get color map and norm.
    auto color_map  = options.color_map.color_map();
    auto color_norm = options.color_norm.get_sequential_norm();

    // Out base file name
    auto const base_fn = options.file_output.out_dir() + cluster_tree_basepath( options, k );

    // As we used the imbalances for the actual clustering, there is no mass tree that
    // we can use here yet. So, we need to add up masses for each cluster.
    // First, prepare storage for these.
    auto centroid_masses = std::vector<std::vector<double>>( k );
    for( size_t ik = 0; ik < k; ++ik ) {
        centroid_masses[ik] = std::vector<double>( profile.tree.edge_count(), 0.0 );
    }

    // Now, add up masses.
    for( size_t sample_idx = 0; sample_idx < assignments.size(); ++sample_idx ) {

        // Which cluster does the sample belong to?
        auto const assn = assignments[ sample_idx ];
        assert( assn < k );
        assert( sample_idx < profile.edge_masses.rows() );
        assert( profile.edge_masses.cols() == centroid_masses[assn].size() );

        // Add its masses to the centroid.
        for( size_t i = 0; i < profile.tree.edge_count(); ++i ) {
            centroid_masses[assn][i] += profile.edge_masses( sample_idx, i );
        }
    }

    // Now, each centroid contains the masses of all samples assigned to it.
    // Write them to tree files.
    for( size_t ci = 0; ci < k; ++ci ) {

        // Prepare colors
        auto const& masses = centroid_masses[ci];
        color_norm->autoscale_max( masses );

        // Now, make a color vector and write to files.
        auto const colors = color_map( *color_norm, masses );
        auto const cntr_fn = base_fn + std::to_string( ci );
        options.tree_output.write_tree_to_files(
            profile.tree,
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

void run_ikmeans( IkmeansOptions const& options )
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

    // Check for existing files.
    check_kmeans_output_files( options );

    // Read input data into imbalances matrix. Filter columns.
    auto profile = options.jplace_input.placement_profile( true );
    auto const columns = epca_filter_constant_columns( profile.edge_imbalances, 0.001 );

    // Move data into vectors, as this is what the kmeans needs.
    // It is a bit stupid, but otherweise the const column filtering would be weird...
    auto edge_imb_vec = std::vector<std::vector<double>>();
    edge_imb_vec.resize( profile.edge_imbalances.rows() );
    for( size_t i = 0; i < profile.edge_imbalances.rows(); ++i ) {
        edge_imb_vec[i].resize( profile.edge_imbalances.cols() );
        for( size_t j = 0; j < profile.edge_imbalances.cols(); ++j ) {
            edge_imb_vec[i][j] = profile.edge_imbalances(i, j);
        }
    }

    // Set up kmeans.
    auto ikmeans = EuclideanKmeans( profile.edge_imbalances.cols() );
    ikmeans.report_iteration = [&]( size_t iteration ){
        if( global_options.verbosity() >= 2 ) {
            std::cout << " - Iteration " << iteration << "\n";
        }
    };

    // Run kmeans for every specified k.
    auto const ks = get_k_values( options );
    std::vector<KmeansClusterOverview> overview;
    for( auto const& k : ks ) {

        // Run it.
        std::cout << "\nRunning Imbalance Kmeans with k=" << k << "\n";
        auto const iterations = ikmeans.run( edge_imb_vec, k );
        auto const clust_info = ikmeans.cluster_info( edge_imb_vec );
        std::cout << "Finished after " << iterations << " iterations\n";

        // Write output.
        write_assignment_file( options, ikmeans.assignments(), clust_info, k );
        write_ikmeans_cluster_trees( options, profile, ikmeans.assignments(), k );

        // Print some cluster info, and collect it for the overview file.
        auto const ci = print_cluster_info( options, ikmeans.assignments(), clust_info, k );
        overview.push_back( ci );
    }

    // Write the overview file for elbow plots etc.
    write_overview_file( options, overview );
}
