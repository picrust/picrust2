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

#include "options/matrix_output.hpp"

#include "genesis/utils/containers/matrix/operators.hpp"
#include "genesis/utils/core/fs.hpp"
#include "genesis/utils/text/string.hpp"

#include <iostream>
#include <stdexcept>

// =================================================================================================
//      Setup Functions
// =================================================================================================

void MatrixOutputOptions::add_matrix_output_opts_to_app( CLI::App* sub, std::string const& name )
{
    add_output_dir_opt_to_app( sub, name, ".", "Matrix Output" );
    add_file_prefix_opt_to_app( sub, name, name, "Matrix Output" );

    // sub->add_set_ignore_case(
    //     "--" + name + "-matrix-format",
    //     format_,
    //     { "matrix", "triangular", "list" },
    //     "Format of the output file.",
    //     true
    // )->group( "Matrix Output" );

    // sub->add_flag(
    //     "--omit-" + name + "-matrix-labels",
    //     omit_labels_,
    //     "If set, the matrix is written without column and row labels."
    // );
}

// =================================================================================================
//      Run Functions
// =================================================================================================

std::string MatrixOutputOptions::output_filename() const
{
    return file_prefix() + ".csv";
}

void MatrixOutputOptions::write_matrix( genesis::utils::Matrix<double> const& mat ) const
{
    using namespace genesis;
    using namespace genesis::utils;

    // TODO offer function for checking nonexistent out file beforehand
    // TODO offer function for dataframe, for other matrix types etc
    // TODO add double presicison
    // TODO add separator char
    // TODO implement other output formats (see below)

    auto const filename = out_dir() + output_filename();

    if( format_ == "matrix" ) {

        // TODO directly write to stream instead of intermediate string!
        utils::file_write( utils::to_string( mat ), filename );

    // } else if( format_ == "triangular" ) {
    //
    //
    // } else if( format_ == "list" ) {


    } else {
        throw CLI::ValidationError(
            "--svg-tree-shape", "Invalid format '" + format_ + "'."
        );
    }
}
