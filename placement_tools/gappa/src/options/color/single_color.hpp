#ifndef GAPPA_OPTIONS_COLOR_SINGLE_H_
#define GAPPA_OPTIONS_COLOR_SINGLE_H_

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

#include <string>
#include <vector>

// =================================================================================================
//      Color Options
// =================================================================================================

/**
 * @brief Helper class to add command line parameter to specify a single color.
 */
class SingleColorOptions
{
public:

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    SingleColorOptions() = default;
    ~SingleColorOptions() = default;

    SingleColorOptions( SingleColorOptions const& other ) = default;
    SingleColorOptions( SingleColorOptions&& )            = default;

    SingleColorOptions& operator= ( SingleColorOptions const& other ) = default;
    SingleColorOptions& operator= ( SingleColorOptions&& )            = default;

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    /**
     * @brief Add an option `--name-color` to the app.
     */
    CLI::Option* add_single_color_opt_to_app(
        CLI::App* sub,
        std::string const& name,
        std::string const& default_color
    );

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

public:

    /**
     * @brief Get the color that was provided by the user.
     */
    genesis::utils::Color color() const;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    std::string name_;
    std::string color_param_;

public:

    CLI::Option* color_option = nullptr;

};

#endif // include guard
