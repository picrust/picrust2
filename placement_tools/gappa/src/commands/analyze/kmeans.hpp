#ifndef GAPPA_COMMANDS_ANALYZE_KMEANS_H_
#define GAPPA_COMMANDS_ANALYZE_KMEANS_H_

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

#include "CLI/CLI.hpp"

#include "options/color/color_map.hpp"
#include "options/color/color_norm.hpp"
#include "options/file_output.hpp"
#include "options/jplace_input.hpp"
#include "options/tree_output.hpp"

#include "genesis/utils/math/kmeans.hpp"

#include <string>
#include <vector>

// =================================================================================================
//      Options
// =================================================================================================

class KmeansOptions
{
public:

    std::string ks;
    bool        overview_file;

    JplaceInputOptions jplace_input;
    ColorMapOptions    color_map;
    ColorNormOptions   color_norm;
    FileOutputOptions  file_output;
    TreeOutputOptions  tree_output;
};

struct KmeansClusterOverview
{
    size_t k;
    double avg_distance;
    double avg_variance;
};

// =================================================================================================
//      Functions
// =================================================================================================

void setup_kmeans(
    KmeansOptions* opt,
    CLI::App* app,
    std::string const& file_prefix
);

std::string assignment_filepath(
    KmeansOptions const& options,
    size_t k
);

std::string cluster_tree_basepath(
    KmeansOptions const& options,
    size_t k
);

void check_kmeans_output_files(
    KmeansOptions const& options
);

std::vector<size_t> get_k_values(
    KmeansOptions const& options
);

void write_assignment_file(
    KmeansOptions const& options,
    std::vector<size_t> const& assignments,
    genesis::utils::KmeansClusteringInfo const& cluster_info,
    size_t k
);

KmeansClusterOverview print_cluster_info(
    KmeansOptions const& options,
    std::vector<size_t> const& assignments,
    genesis::utils::KmeansClusteringInfo const& cluster_info,
    size_t k
);

void write_overview_file(
    KmeansOptions const& options,
    std::vector<KmeansClusterOverview> const& overview
);

#endif // include guard
