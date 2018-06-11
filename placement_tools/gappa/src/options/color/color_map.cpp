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

#include "options/color/color_map.hpp"

#include "options/global.hpp"

#include "genesis/utils/tools/color/diverging_lists.hpp"
#include "genesis/utils/tools/color/qualitative_lists.hpp"
#include "genesis/utils/tools/color/sequential_lists.hpp"
#include "genesis/utils/tools/color/functions.hpp"
#include "genesis/utils/text/string.hpp"
#include "genesis/utils/core/fs.hpp"
#include "genesis/utils/tools/color/names.hpp"

#include <stdexcept>

// =================================================================================================
//      Constructor
// =================================================================================================

ColorMapOptions::ColorMapOptions()
{
    using namespace genesis::utils;
    mask_color_param_ = color_to_hex( color_map_.mask_color() );
    over_color_param_ = color_to_hex( color_map_.over_color() );
    under_color_param_ = color_to_hex( color_map_.under_color() );
}

// =================================================================================================
//      Setup Functions
// =================================================================================================

CLI::Option* ColorMapOptions::add_color_list_opt_to_app(
    CLI::App* sub,
    std::string const& default_color_list
) {
    // Correct setup check.
    if( color_list_option != nullptr ) {
        throw std::domain_error( "Cannot use the same ColorMapOptions object multiple times." );
    }

    // Set Default
    palette_param_ = default_color_list;

    // Color List
    color_list_option = sub->add_option(
        "--color-list",
        palette_param_,
        "List of colors to use for the palette. Can either be the name of a color list, "
        "a file containing one color per line, or an actual list of colors.",
        true
    );
    color_list_option->group( "Color" );

    // Reverse
    reverse_color_list_option = sub->add_flag_function(
        "--reverse-color-list",
        [this]( size_t ){
            color_map_.reverse( true );
        },
        "If set, the --color-list is reversed."
    )->group( "Color" );

    return color_list_option;
}

CLI::Option* ColorMapOptions::add_under_color_opt_to_app( CLI::App* sub )
{
    // Clip Under
    clip_under_option = sub->add_flag_function(
        "--clip-under",
        [this]( size_t ){
            color_map_.clip_under( true );
        },
        "Clip (clamp) values less than min to be inside [ min, max ]. "
        "If set, --under-color is not used to indicate values out of range."
    )->group( "Color" );

    // Under Color
    under_color_option = sub->add_option(
        "--under-color",
        under_color_param_,
        "Color used to indicate values below min.",
        true
    )->group( "Color" );

    // Special: If we also use over color, we can offer a clip option shortcut.
    if( clip_over_option ) {
        clip_option = sub->add_flag_function(
            "--clip",
            [this]( size_t ){
                color_map_.clip( true );
            },
            "Clip (clamp) values to be inside [ min, max ]. "
            "This option is a shortcut to set --clip-under and --clip-over at once."
        )->group( "Color" );
    }

    return under_color_option;
}

CLI::Option* ColorMapOptions::add_over_color_opt_to_app( CLI::App* sub )
{
    // Clip Over
    clip_over_option = sub->add_flag_function(
        "--clip-over",
        [this]( size_t ){
            color_map_.clip_over( true );
        },
        "Clip (clamp) values greater than max to be inside [ min, max ]. "
        "If set, --over-color is not used to indicate values out of range."
    )->group( "Color" );

    // Over Color
    over_color_option = sub->add_option(
        "--over-color",
        over_color_param_,
        "Color used to indicate values above max.",
        true
    )->group( "Color" );

    // Special: If we also use under color, we can offer a clip option shortcut.
    if( clip_under_option ) {
        clip_option = sub->add_flag_function(
            "--clip",
            [this]( size_t ){
                color_map_.clip( true );
            },
            "Clip (clamp) values to be inside [ min, max ]. "
            "This option is a shortcut to set --clip-under and --clip-over at once."
        )->group( "Color" );
    }

    return over_color_option;
}

CLI::Option* ColorMapOptions::add_mask_color_opt_to_app( CLI::App* sub )
{
    // Mask Color
    mask_color_option = sub->add_option(
        "--mask-color",
        mask_color_param_,
        "Color used to indicate masked values.",
        true
    )->group( "Color" );

    return mask_color_option;
}

// =================================================================================================
//      Helper Functions
// =================================================================================================

genesis::utils::Color ColorMapOptions::resolve_color_string( std::string color_str, std::string const& param_name ) const
{
    try {
        return genesis::utils::resolve_color_string( color_str );
    } catch( std::exception& ex ) {
        throw CLI::ValidationError(
            param_name, "Invalid color '" + color_str + "': " +
            std::string( ex.what() )
        );
    }
}

std::vector<genesis::utils::Color> ColorMapOptions::resolve_color_list(
    std::vector<std::string> const& list,
    std::string const& param_name
) const {
    std::vector<genesis::utils::Color> result;

    for( auto const& color_str : list ) {
        result.push_back( resolve_color_string( color_str, param_name ));
    }

    return result;
}

// =================================================================================================
//      Run Functions
// =================================================================================================

genesis::utils::ColorMap const& ColorMapOptions::color_map() const
{
    using namespace genesis::utils;

    // If map was already filled in previous call, just return it.
    if( ! color_map_.empty() ) {
        return color_map_;
    }

    // Resolve special colors.
    color_map_.under_color( resolve_color_string( under_color_param_, "--under-color" ) );
    color_map_.over_color( resolve_color_string( over_color_param_, "--over-color" ) );
    color_map_.mask_color( resolve_color_string( mask_color_param_, "--mask-color" ) );

    // Now resolve the actual color list.
    // Try color list names first.
    if( contains_ci( diverging_color_list_names(), palette_param_ )) {
        color_map_.palette( diverging_color_list( palette_param_ ));
        return color_map_;
    }
    if( contains_ci( qualitative_color_list_names(), palette_param_ )) {
        color_map_.palette( qualitative_color_list( palette_param_ ));
        return color_map_;
    }
    if( contains_ci( sequential_color_list_names(), palette_param_ )) {
        color_map_.palette( sequential_color_list( palette_param_ ));
        return color_map_;
    }

    // Now check if it is a file with colors.
    if( is_file( palette_param_ ) ) {
        auto const list = split( file_read( palette_param_ ), "\n\r", true );
        color_map_.palette( resolve_color_list( list, "--color-list" ));
        return color_map_;
    }

    // Finally, treat it as a comma separated list of colors.
    auto const list = split( palette_param_, ",", true );
    color_map_.palette( resolve_color_list( list, "--color-list" ));
    return color_map_;
}
