#ifndef GAPPA_OPTIONS_COLOR_COLOR_NORM_H_
#define GAPPA_OPTIONS_COLOR_COLOR_NORM_H_

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

#include "genesis/utils/tools/color.hpp"
#include "genesis/utils/tools/color/normalization.hpp"
#include "genesis/utils/tools/color/norm_diverging.hpp"
#include "genesis/utils/tools/color/norm_linear.hpp"
#include "genesis/utils/tools/color/norm_logarithmic.hpp"

#include <limits>
#include <memory>
#include <string>
#include <vector>

// =================================================================================================
//      Color Norm Options
// =================================================================================================

/**
 * @brief Helper class to add command line parameter to use a color normalization.
 */
class ColorNormOptions
{
public:

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    ColorNormOptions()  = default;
    ~ColorNormOptions() = default;

    ColorNormOptions( ColorNormOptions const& other ) = default;
    ColorNormOptions( ColorNormOptions&& )            = default;

    ColorNormOptions& operator= ( ColorNormOptions const& other ) = default;
    ColorNormOptions& operator= ( ColorNormOptions&& )            = default;

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    CLI::Option* add_log_scaling_opt_to_app( CLI::App* sub );
    CLI::Option* add_min_value_opt_to_app( CLI::App* sub );
    CLI::Option* add_mid_value_opt_to_app( CLI::App* sub );
    CLI::Option* add_max_value_opt_to_app( CLI::App* sub );
    CLI::Option* add_mask_value_opt_to_app( CLI::App* sub );

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

public:

    bool log_scaling() const
    {
        return log_scaling_;
    }

    std::unique_ptr<genesis::utils::ColorNormalizationLinear> get_sequential_norm() const;
    genesis::utils::ColorNormalizationDiverging get_diverging_norm() const;

    void apply_options( genesis::utils::ColorNormalizationLinear& norm ) const;
    void apply_options( genesis::utils::ColorNormalizationLogarithmic& norm ) const;
    void apply_options( genesis::utils::ColorNormalizationDiverging& norm ) const;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    bool log_scaling_ = false;
    double min_value_ = 0.0;
    double mid_value_ = 0.5;
    double max_value_ = 1.0;
    double mask_value_ = std::numeric_limits<double>::quiet_NaN();

public:

    CLI::Option* log_scaling_option = nullptr;
    CLI::Option* min_value_option = nullptr;
    CLI::Option* mid_value_option = nullptr;
    CLI::Option* max_value_option = nullptr;
    CLI::Option* mask_value_option = nullptr;

};

#endif // include guard
