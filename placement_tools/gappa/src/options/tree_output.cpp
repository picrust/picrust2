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

#include "options/tree_output.hpp"

#include "genesis/tree/drawing/functions.hpp"
#include "genesis/utils/text/string.hpp"
#include "genesis/utils/tools/color/functions.hpp"
#include "genesis/utils/tools/color/helpers.hpp"
#include "genesis/utils/tools/tickmarks.hpp"

#include <iostream>
#include <stdexcept>

// =================================================================================================
//      Setup Functions
// =================================================================================================

void TreeOutputOptions::add_tree_output_opts_to_app( CLI::App* sub )
{
    sub->add_flag(
        "--write-newick-tree",
        write_newick_tree_,
        "If set, the tree is written to a Newick file."
    )->group( "Tree Output" );
    sub->add_flag(
        "--write-nexus-tree",
        write_nexus_tree_,
        "If set, the tree is written to a Nexus file."
    )->group( "Tree Output" );
    sub->add_flag(
        "--write-phyloxml-tree",
        write_phyloxml_tree_,
        "If set, the tree is written to a Phyloxml file."
    )->group( "Tree Output" );
    sub->add_flag(
        "--write-svg-tree",
        write_svg_tree_,
        "If set, the tree is written to a Svg file."
    )->group( "Tree Output" );

    svg_tree_output.add_svg_tree_output_opts_to_app( sub );
}

// =================================================================================================
//      Run Functions
// =================================================================================================

void TreeOutputOptions::check_tree_formats() const
{
    if( ! write_newick_tree_ && ! write_nexus_tree_ && ! write_phyloxml_tree_ && ! write_svg_tree_ ) {
        std::cout << "Warning: You did not specify any tree output format. ";
        std::cout << "Thus, no tree files will be written. ";
        std::cout << "In order to specify the wanted formats, use the --write-...-tree options.\n";
    }
}

std::vector<std::string> TreeOutputOptions::get_extensions() const
{
    std::vector<std::string> res;

    if( write_newick_tree_ ) {
        res.push_back( "newick" );
    }
    if( write_nexus_tree_ ) {
        res.push_back( "nexus" );
    }
    if( write_phyloxml_tree_ ) {
        res.push_back( "phyloxml" );
    }
    if( write_svg_tree_ ) {
        res.push_back( "svg" );
    }

    return res;
}

void TreeOutputOptions::write_tree_to_files(
    genesis::tree::DefaultTree const&         tree,
    std::string const&                        file_path_prefix
) const {
    using namespace genesis::tree;

    if( write_newick_tree_ ) {
        write_tree_to_newick_file( tree, file_path_prefix + ".newick" );
    }

    if( write_nexus_tree_ ) {
        write_tree_to_nexus_file( tree, file_path_prefix + ".nexus" );
    }

    if( write_phyloxml_tree_ ) {
        write_tree_to_phyloxml_file( tree, file_path_prefix + ".phyloxml" );
    }

    if( write_svg_tree_ ) {
        write_tree_to_svg_file(
            tree,
            svg_tree_output.layout_parameters(),
            file_path_prefix + ".svg"
        );
    }
}

void TreeOutputOptions::write_tree_to_files(
    genesis::tree::DefaultTree const&         tree,
    std::vector<genesis::utils::Color> const& color_per_branch,
    std::string const&                        file_path_prefix
) const {
    using namespace genesis::tree;

    if( write_newick_tree_ ) {
        std::cout << "Warning: Option --write-newick-tree is set, but the output contains colors, ";
        std::cout << "which are not available in the Newick format. ";
        std::cout << "The Newick tree only contains the topology of the tree with names and branch lengths. ";
        std::cout << "Use another format to get a colored tree!\n";

        write_tree_to_newick_file( tree, file_path_prefix + ".newick" );
    }

    if( write_nexus_tree_ ) {
        write_color_tree_to_nexus_file( tree, color_per_branch, file_path_prefix + ".nexus" );
    }

    if( write_phyloxml_tree_ ) {
        write_color_tree_to_phyloxml_file( tree, color_per_branch, file_path_prefix + ".phyloxml" );
    }

    if( write_svg_tree_ ) {
        write_color_tree_to_svg_file(
            tree,
            svg_tree_output.layout_parameters(),
            color_per_branch,
            file_path_prefix + ".svg"
        );
    }
}

void TreeOutputOptions::write_tree_to_files(
    genesis::tree::DefaultTree const&         tree,
    std::vector<genesis::utils::Color> const& color_per_branch,
    genesis::utils::ColorMap const&           color_map,
    genesis::utils::ColorNormalization const& color_norm,
    std::string const&                        file_path_prefix
) const {
    using namespace genesis::tree;
    using namespace genesis::utils;

    // In case we output a non svg tree, we need to report colors and tickmarks,
    // as they are not available in the other formats.
    bool print_legend = false;

    if( write_newick_tree_ ) {
        std::cout << "Warning: Option --write-newick-tree is set, but the output contains colors, ";
        std::cout << "which are not available in the Newick format. ";
        std::cout << "The Newick tree only contains the topology of the tree with names and branch lengths. ";
        std::cout << "Use another format to get a colored tree!\n";

        write_tree_to_newick_file( tree, file_path_prefix + ".newick" );
    }

    if( write_nexus_tree_ ) {
        write_color_tree_to_nexus_file( tree, color_per_branch, file_path_prefix + ".nexus" );
        print_legend = true;
    }

    if( write_phyloxml_tree_ ) {
        write_color_tree_to_phyloxml_file( tree, color_per_branch, file_path_prefix + ".phyloxml" );
        print_legend = true;
    }

    if( write_svg_tree_ ) {
        write_color_tree_to_svg_file(
            tree,
            svg_tree_output.layout_parameters(),
            color_per_branch,
            color_map,
            color_norm,
            file_path_prefix + ".svg"
        );
    }

    if( print_legend ) {
        // TODO maybe make the num ticks changable. if so, also use it for the svg output!
        auto const tickmarks = color_tickmarks( color_norm, 5 );

        std::cout << "Output options --write-nexus-tree and --write-phyloxml-tree produce trees ";
        std::cout << "with colored branches; these formats are however not able to store the legend, ";
        std::cout << "that is, which color represents which value. ";
        std::cout << "Thus, use to following positions to create a legend. ";
        std::cout << "These positions range from 0.0 (lowest) to 1.0 (heighest), and are labelled ";
        std::cout << "with the values and colors represented by those positions.\n";

        for( auto const& tick : tickmarks ) {
            auto const rel_pos = tick.first;
            auto label = tick.second;

            if( rel_pos == 0.0 && color_map.clip_under() ) {
                label = "≤ " + label;
            }
            if( rel_pos == 1.0 && color_map.clip_over() ) {
                label = "≥ " + label;
            }

            auto const col_str = color_to_hex( color_map( rel_pos ));
            std::cout << "    At " << to_string_precise( rel_pos, 3 ) << ": Label '" << label << "', Color " << col_str << "\n";
        }

        std::cout << "Alternatively, use the option --write-svg-tree to create an Svg file ";
        std::cout << "from which the color legend can be copied.\n";
    }
}
