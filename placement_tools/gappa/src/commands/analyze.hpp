#ifndef GAPPA_COMMANDS_ANALYZE_H_
#define GAPPA_COMMANDS_ANALYZE_H_

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

#include "commands/analyze/assign.hpp"
#include "commands/analyze/correlation.hpp"
#include "commands/analyze/dispersion.hpp"
#include "commands/analyze/graft.hpp"
#include "commands/analyze/kmeans_imbalance.hpp"
#include "commands/analyze/kmeans_phylogenetic.hpp"
#include "commands/analyze/krd.hpp"
#include "commands/analyze/nhd.hpp"
#include "commands/analyze/squash.hpp"
#include "commands/analyze/visualize_color.hpp"

#include <string>
#include <vector>

// =================================================================================================
//      Functions
// =================================================================================================

void setup_analyze( CLI::App& app )
{
    // Create the module subcommand objects.
    auto sub = app.add_subcommand(
        "analyze",
        "Commands for analyzing and visualizing placement data."
    );
    sub->require_subcommand( 1 );

    // Add module subcommands.
    // setup_krd( *sub );
    // setup_nhd( *sub );

    setup_assign( *sub );
    setup_correlation( *sub );
    setup_dispersion( *sub );
    setup_graft( *sub );
    setup_ikmeans( *sub );
    setup_pkmeans( *sub );
    setup_squash( *sub );
    setup_visualize_color( *sub );
}

#endif // include guard
