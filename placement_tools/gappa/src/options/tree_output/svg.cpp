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

#include "options/tree_output/svg.hpp"

#include <iostream>
#include <stdexcept>

// =================================================================================================
//      Setup Functions
// =================================================================================================

void SvgTreeOutputOptions::add_svg_tree_output_opts_to_app( CLI::App* sub )
{
    sub->add_set_ignore_case(
        "--svg-tree-shape",
        shape_,
        { "circular", "rectangular" },
        "Shape of the tree.",
        // "Shape of the tree, 'circular' or 'rectangular'.",
        true
    )->group( "Svg Tree Output" );

    sub->add_set_ignore_case(
        "--svg-tree-type",
        type_,
        { "cladogram", "phylogram" },
        "Type of the tree.",
        // "Type of the tree, 'cladogram' or 'phylogram'.",
        true
    )->group( "Svg Tree Output" );

    sub->add_option(
        "--svg-tree-stroke-width",
        stroke_width_,
        "Svg stroke width for the branches of the tree.",
        true
    )->group( "Svg Tree Output" );

    sub->add_flag(
        "--svg-tree-ladderize",
        ladderize_,
        "If set, the tree is ladderized."
    )->group( "Svg Tree Output" );
}

// =================================================================================================
//      Run Functions
// =================================================================================================

genesis::tree::LayoutParameters SvgTreeOutputOptions::layout_parameters() const
{
    using namespace genesis;
    using namespace genesis::tree;

    LayoutParameters res;

    if( shape_ == "circular" ) {
        res.shape = LayoutShape::kCircular;
    } else if( shape_ == "rectangular" ) {
        res.shape = LayoutShape::kRectangular;
    } else {
        throw CLI::ValidationError(
            "--svg-tree-shape", "Invalid shape '" + shape_ + "'."
        );
    }

    if( type_ == "cladogram" ) {
        res.type = LayoutType::kCladogram;
    } else if( type_ == "phylogram" ) {
        res.type = LayoutType::kPhylogram;
    } else {
        throw CLI::ValidationError(
            "--svg-tree-type", "Invalid type '" + type_ + "'."
        );
    }

    if( stroke_width_ <= 0.0 ) {
        throw CLI::ValidationError(
            "--svg-tree-stroke-width",
            "Svg stroke width has to be positive."
        );
    }

    res.stroke.width = stroke_width_;
    res.ladderize = ladderize_;
    return res;
}
