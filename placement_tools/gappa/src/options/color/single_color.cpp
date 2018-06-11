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

#include "options/color/single_color.hpp"

#include "options/global.hpp"

#include "genesis/utils/tools/color/functions.hpp"
#include "genesis/utils/text/string.hpp"
#include "genesis/utils/tools/color/names.hpp"

#include <stdexcept>

// =================================================================================================
//      Setup Functions
// =================================================================================================

CLI::Option* SingleColorOptions::add_single_color_opt_to_app(
    CLI::App* sub,
    std::string const& name,
    std::string const& default_color
) {
    // Correct setup check.
    if( color_option != nullptr ) {
        throw std::domain_error( "Cannot use the same SingleColorOptions object multiple times." );
    }

    // Set Default
    color_param_ = default_color;
    name_ = name;

    // Color List
    color_option = sub->add_option(
        "--" + name + "-color",
        color_param_,
        "Colors to use for " + name + ".",
        true
    );
    color_option->group( "Color" );

    return color_option;
}

// =================================================================================================
//      Run Functions
// =================================================================================================

genesis::utils::Color SingleColorOptions::color() const
{
    try {
        return genesis::utils::resolve_color_string( color_param_ );
    } catch( std::exception& ex ) {
        throw CLI::ValidationError(
            name_, "Invalid color '" + color_param_ + "': " +
            std::string( ex.what() )
        );
    }
}
