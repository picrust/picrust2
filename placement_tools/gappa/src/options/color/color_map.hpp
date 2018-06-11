#ifndef GAPPA_OPTIONS_COLOR_COLOR_MAP_H_
#define GAPPA_OPTIONS_COLOR_COLOR_MAP_H_

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
#include "genesis/utils/tools/color/map.hpp"

#include <limits>
#include <memory>
#include <string>
#include <vector>

// =================================================================================================
//      Color Map Options
// =================================================================================================

/**
 * @brief Helper class to add command line parameter to use a color map,
 * that is, to select colors for output.
 */
class ColorMapOptions
{
public:

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    ColorMapOptions();
    ~ColorMapOptions() = default;

    ColorMapOptions( ColorMapOptions const& other ) = default;
    ColorMapOptions( ColorMapOptions&& )            = default;

    ColorMapOptions& operator= ( ColorMapOptions const& other ) = default;
    ColorMapOptions& operator= ( ColorMapOptions&& )            = default;

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    /**
     * @brief Add an option to input a list of colors to an app, e.g., in order to customize
     * the color scheme of gradients.
     *
     * The function takes either a name of a list, a file containing colors, or a comma separated
     * list of colors. Colors can be specified as names (web and xkcd), or as hex, using a leading #.
     */
    CLI::Option* add_color_list_opt_to_app(
        CLI::App* sub,
        std::string const& default_color_list
    );

    CLI::Option* add_under_color_opt_to_app( CLI::App* sub );
    CLI::Option* add_over_color_opt_to_app( CLI::App* sub );
    CLI::Option* add_mask_color_opt_to_app( CLI::App* sub );

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

public:

    /**
     * @brief Get the color map with all settings applied that were provided by the user.
     */
    genesis::utils::ColorMap const& color_map() const;

    // -------------------------------------------------------------------------
    //     Helper Functions
    // -------------------------------------------------------------------------

private:

    /**
     * @brief Helper function that wraps the genesis function of the same name,
     * but offers a nicer error feedback.
     */
    genesis::utils::Color resolve_color_string(
        std::string color_str,
        std::string const& param_name
    ) const;

    /**
     * @brief Same as resolve_color_string(), but for a whole list of colors.
     */
    std::vector<genesis::utils::Color> resolve_color_list(
        std::vector<std::string> const& list,
        std::string const& param_name
    ) const;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    // Helper members that store the user input for (list of) colors.
    // We need this, because this cannot bind directly to the properties of the color objects,
    // and because we take colors as strings of differnet format and need to convert first.
    std::string palette_param_;
    std::string under_color_param_;
    std::string over_color_param_;
    std::string mask_color_param_;

    mutable genesis::utils::ColorMap color_map_;

public:

    CLI::Option* color_list_option = nullptr;
    CLI::Option* reverse_color_list_option = nullptr;

    CLI::Option* under_color_option = nullptr;
    CLI::Option* clip_under_option = nullptr;

    CLI::Option* over_color_option = nullptr;
    CLI::Option* clip_over_option = nullptr;

    CLI::Option* clip_option = nullptr;
    CLI::Option* mask_color_option = nullptr;

};

#endif // include guard
