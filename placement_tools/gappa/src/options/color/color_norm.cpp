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

#include "options/color/color_norm.hpp"

#include "options/global.hpp"

#include "genesis/utils/core/std.hpp"

#include <stdexcept>

// =================================================================================================
//      Setup Functions
// =================================================================================================

CLI::Option* ColorNormOptions::add_log_scaling_opt_to_app( CLI::App* sub )
{
    // Correct setup check.
    if( log_scaling_option != nullptr ) {
        throw std::domain_error( "Cannot use the same ColorNormOptions object multiple times." );
    }

    sub->add_flag(
        "--log-scaling",
        log_scaling_,
        "If set, the sequential color list is logarithmically scaled instead of linearily."
    )->group( "Color" );

    return log_scaling_option;
}

CLI::Option* ColorNormOptions::add_min_value_opt_to_app( CLI::App* sub )
{
    // Correct setup check.
    if( min_value_option != nullptr ) {
        throw std::domain_error( "Cannot use the same ColorNormOptions object multiple times." );
    }

    // Min
    min_value_option = sub->add_option(
        "--min-value",
        min_value_,
        "Minimum value that is represented by the color scale. "
        "If not set, the minimum value in the data is used."
    )->group( "Color" );

    return min_value_option;
}

CLI::Option* ColorNormOptions::add_mid_value_opt_to_app( CLI::App* sub )
{
    // Correct setup check.
    if( mid_value_option != nullptr ) {
        throw std::domain_error( "Cannot use the same ColorNormOptions object multiple times." );
    }

    // Min
    mid_value_option = sub->add_option(
        "--mid-value",
        mid_value_,
        "Mid value that is represented by the diverging color scale. "
        "If not set, the mid value in the data is used."
    )->group( "Color" );

    return mid_value_option;
}

CLI::Option* ColorNormOptions::add_max_value_opt_to_app( CLI::App* sub )
{
    // Correct setup check.
    if( max_value_option != nullptr ) {
        throw std::domain_error( "Cannot use the same ColorNormOptions object multiple times." );
    }

    // Max
    max_value_option = sub->add_option(
        "--max-value",
        max_value_,
        "Maximum value that is represented by the color scale. "
        "If not set, the maximum value in the data is used."
    )->group( "Color" );

    return max_value_option;
}

CLI::Option* ColorNormOptions::add_mask_value_opt_to_app( CLI::App* sub )
{
    // Correct setup check.
    if( mask_value_option != nullptr ) {
        throw std::domain_error( "Cannot use the same ColorNormOptions object multiple times." );
    }

    // Mask
    mask_value_option = sub->add_option(
        "--mask-value",
        mask_value_,
        "Mask value that identifies invalid values. "
        "Value in the data that compare equal to the mask value are colored using --mask-color. "
        "This is meant as a simple means of filtering and visualizing invalid values. "
        "If not set, no masking value is applied."
    )->group( "Color" );

    return mask_value_option;
}

// =================================================================================================
//      Run Functions
// =================================================================================================

std::unique_ptr<genesis::utils::ColorNormalizationLinear> ColorNormOptions::get_sequential_norm() const
{
    using namespace genesis::utils;
    std::unique_ptr<ColorNormalizationLinear> res;

    if( log_scaling_ ) {
        res = make_unique<ColorNormalizationLogarithmic>();
        apply_options( static_cast<ColorNormalizationLogarithmic&>( *res ) );
    } else {
        res = make_unique<ColorNormalizationLinear>();
        apply_options( static_cast<ColorNormalizationLinear&>( *res ));
    }

    return res;
}

genesis::utils::ColorNormalizationDiverging ColorNormOptions::get_diverging_norm() const
{
    using namespace genesis::utils;
    auto res = ColorNormalizationDiverging();
    apply_options( res );
    return res;
}

void ColorNormOptions::apply_options( genesis::utils::ColorNormalizationLinear& norm ) const
{
    // CLI objects evaluate to true if the option was passed by the user.
    // So here, we first test wheter the option object was actually created,
    // that is, whether the command uses it, and then, whether the user also specified a value.
    // Only if both are true, we use the value to overwrite the norm value.

    if( min_value_option && *min_value_option ) {
        norm.min_value( min_value_ );
    }
    if( max_value_option && *max_value_option ) {
        norm.max_value( max_value_ );
    }
    if( mask_value_option && *mask_value_option ) {
        norm.mask_value( mask_value_ );
    }
}

void ColorNormOptions::apply_options( genesis::utils::ColorNormalizationLogarithmic& norm ) const
{
    apply_options( static_cast<genesis::utils::ColorNormalizationLinear&>( norm ));
}

void ColorNormOptions::apply_options( genesis::utils::ColorNormalizationDiverging& norm ) const
{
    // First apply base class.
    apply_options( static_cast<genesis::utils::ColorNormalizationLinear&>( norm ));

    // Then special divergent options.
    if( mid_value_option && *mid_value_option ) {
        norm.mid_value( mid_value_ );
    }
}
