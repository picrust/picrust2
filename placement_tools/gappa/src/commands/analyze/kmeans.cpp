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

#include "commands/analyze/kmeans.hpp"

#include "options/global.hpp"

#include "CLI/CLI.hpp"

#include "genesis/utils/io/output_stream.hpp"
#include "genesis/utils/text/string.hpp"

#ifdef GENESIS_OPENMP
#   include <omp.h>
#endif

// =================================================================================================
//      Setup
// =================================================================================================

void setup_kmeans(
    KmeansOptions* opt,
    CLI::App* app,
    std::string const& file_prefix
) {
    // Sample input
    opt->jplace_input.add_jplace_input_opt_to_app( app );

    // Number of clusters to find.
    auto k_opt = app->add_option(
        "--k",
        opt->ks,
        "Number of clusters to find. Can be a comma-separated list of multiple values or "
        "ranges for k: 1-5,8,10,12",
        true
    );
    k_opt->group( "Settings" );
    k_opt->required();

    // Overview file.
    auto overview_file_opt = app->add_flag(
        "--write-overview-file",
        opt->overview_file,
        "If provided, a table file is written that summarizes the average distance and "
        "variance of the clusters for each k. Useful for elbow plots."
    );
    overview_file_opt->group( "Settings" );

    // Other jplace settings
    opt->jplace_input.add_point_mass_opt_to_app( app );
    opt->jplace_input.add_ignore_multiplicities_opt_to_app( app );

    // Color.
    opt->color_map.add_color_list_opt_to_app( app, "BuPuBk" );
    opt->color_norm.add_log_scaling_opt_to_app( app );

    // Output files.
    opt->tree_output.add_tree_output_opts_to_app( app );
    opt->file_output.add_output_dir_opt_to_app( app );
    opt->file_output.add_file_prefix_opt_to_app( app, "", file_prefix );
}

// =================================================================================================
//      Output Files
// =================================================================================================

std::string assignment_filepath(
    KmeansOptions const& options,
    size_t k
) {
    return options.file_output.file_prefix() + "k_" + std::to_string( k ) + "_assignments.csv";
}

std::string cluster_tree_basepath(
    KmeansOptions const& options,
    size_t k
) {
    return options.file_output.file_prefix() + "k_" + std::to_string( k ) + "_centroid_";
}

void check_kmeans_output_files(
    KmeansOptions const& options
) {
    auto const ks = get_k_values( options );

    std::vector<std::string> files_to_check;

    // Add assignemnt files.
    for( auto k : ks ) {
        files_to_check.push_back(
            assignment_filepath( options, k )
        );
    }

    // Add cluster tree files if specified.
    // For each k, we need to check all cluster numbers and all tree formats.
    for( auto k : ks ) {
        for( size_t ci = 0; ci < k; ++ci ) {
            for( auto const& e : options.tree_output.get_extensions() ) {
                files_to_check.push_back(
                    cluster_tree_basepath( options, k ) + std::to_string( ci ) + "\\." + e
                );
            }
        }
    }

    // Add overview file if needed.
    if( options.overview_file ) {
        files_to_check.push_back(
            options.file_output.file_prefix() + "overview.csv"
        );
    }

    options.file_output.check_nonexistent_output_files( files_to_check );
}

// =================================================================================================
//      Helpers
// =================================================================================================

std::vector<size_t> get_k_values( KmeansOptions const& options )
{
    using namespace genesis::utils;

    // Prepare an exception with nice text.
    auto excpt = CLI::ValidationError(
        "--k (" + options.ks +  ")",
        "Invalid list of values for k. Needs to be a comma-separated list of positive numbers or "
        "ranges, e.g., 5-10,12,15"
    );

    // Try to get the ks by splitting the list.
    std::vector<size_t> result;
    try {
        result = split_range_list( options.ks );
    } catch ( ... ) {
        throw excpt;
    }

    // Additional condition: need numbers, but no zero.
    if( result.size() == 0 || result[0] == 0 ) {
        throw excpt;
    }
    return result;
}

// =================================================================================================
//      Output
// =================================================================================================

void write_assignment_file(
    KmeansOptions const& options,
    std::vector<size_t> const& assignments,
    genesis::utils::KmeansClusteringInfo const& cluster_info,
    size_t k
) {
    auto const set_size = options.jplace_input.file_count();

    // Saftey
    if( assignments.size() != set_size || cluster_info.distances.size() != set_size ) {
        throw std::runtime_error(
            "Internal Error: Differing number of assignments (" + std::to_string( assignments.size() ) +
            ") and sample set size (" + std::to_string( set_size ) + ")."
        );
    }

    // Prepare assignments file.
    // TODO check with file overwrite settings
    auto const assm_fn = options.file_output.out_dir() + assignment_filepath( options, k );
    std::ofstream assm_os;
    genesis::utils::file_output_stream( assm_fn, assm_os );

    // Write header.
    assm_os << "Sample\tCluster\tDistance\n";

    // Write assignments
    for( size_t fi = 0; fi < set_size; ++fi ) {
        assm_os << options.jplace_input.base_file_name( fi );
        assm_os << "\t" << assignments[fi];
        assm_os << "\t" << cluster_info.distances[fi];
        assm_os << "\n";
    }
}

KmeansClusterOverview print_cluster_info(
    KmeansOptions const& options,
    std::vector<size_t> const& assignments,
    genesis::utils::KmeansClusteringInfo const& cluster_info,
    size_t k
) {
    // Currentky not needed
    (void) options;

    // Accumulate average distanace and variance of all clusters.
    double avg_dst = 0.0;
    double avg_var = 0.0;

    // Iterate clusters/centroids.
    for( size_t ik = 0; ik < k; ++ik ) {

        // Distance and variance of the current cluster.
        double cavg_dst = 0.0;
        size_t cavg_cnt = 0;

        // Do the accumulation for data points that belong to the current cluster.
        for( size_t ai = 0; ai < assignments.size(); ++ai ) {
            if( assignments[ai] != ik ) {
                continue;
            }
            avg_dst  += cluster_info.distances[ai];
            avg_var  += cluster_info.distances[ai] * cluster_info.distances[ai];
            cavg_dst += cluster_info.distances[ai];
            ++cavg_cnt;
        }
        cavg_dst /= static_cast<double>( cavg_cnt );

        std::cout << "Cluster " << ik << ": " << cluster_info.counts[ik] << " samples,";
        std::cout << " with a variance of " << cluster_info.variances[ik];
        std::cout << " and average distance " << cavg_dst << "\n";
    }
    avg_dst /= static_cast<double>( assignments.size() );
    avg_var /= static_cast<double>( assignments.size() );

    // Different calculation for result checking.
    // Both should yield the exact values as the ones we already have.
    // double avg_dst2 = 0.0;
    // double avg_var2 = 0.0;
    // for( size_t i = 0; i < cluster_info.distances.size(); ++i ) {
    //     avg_dst2 += cluster_info.distances[i];
    //     avg_var2  += cluster_info.distances[i] * cluster_info.distances[i];
    // }
    // avg_dst2 /= static_cast<double>( cluster_info.distances.size() );
    // avg_var2  /= static_cast<double>( cluster_info.distances.size() );

    std::cout << "Total average distance: " << avg_dst << "\n";
    std::cout << "Total average variance: " << avg_var << "\n";

    return { k, avg_dst, avg_var };
}

void write_overview_file(
    KmeansOptions const& options,
    std::vector<KmeansClusterOverview> const& overview
) {
    // Only write a tree file if user specified a path.
    if( ! options.overview_file ) {
        return;
    }

    // Prepare assignments file.
    // TODO check with file overwrite settings
    auto const ov_fn
        = options.file_output.out_dir() + options.file_output.file_prefix()
        + "overview.csv"
    ;
    std::ofstream ov_os;
    genesis::utils::file_output_stream( ov_fn, ov_os );

    // Write header.
    ov_os << "k\tDistance\tVariance\n";

    // Write overview
    for( auto const& l : overview ) {
        ov_os << l.k << "\t" << l.avg_distance << "\t" << l.avg_variance << "\n";
    }
}
