#ifndef GAPPA_OPTIONS_TREE_OUTPUT_SVG_H_
#define GAPPA_OPTIONS_TREE_OUTPUT_SVG_H_

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

#include "genesis/tree/drawing/functions.hpp"

#include <string>
#include <vector>

// =================================================================================================
//      Svg Tree Output Options
// =================================================================================================

/**
 * @brief
 */
class SvgTreeOutputOptions
{
public:

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    SvgTreeOutputOptions()  = default;
    ~SvgTreeOutputOptions() = default;

    SvgTreeOutputOptions( SvgTreeOutputOptions const& other ) = default;
    SvgTreeOutputOptions( SvgTreeOutputOptions&& )            = default;

    SvgTreeOutputOptions& operator= ( SvgTreeOutputOptions const& other ) = default;
    SvgTreeOutputOptions& operator= ( SvgTreeOutputOptions&& )            = default;

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    void add_svg_tree_output_opts_to_app( CLI::App* sub );

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

    genesis::tree::LayoutParameters layout_parameters() const;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    std::string shape_   = "circular";
    std::string type_    = "cladogram";
    double stroke_width_ = 5.0;
    bool ladderize_      = true;

};

#endif // include guard
