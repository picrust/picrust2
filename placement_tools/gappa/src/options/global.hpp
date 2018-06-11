#ifndef GAPPA_OPTIONS_GLOBAL_H_
#define GAPPA_OPTIONS_GLOBAL_H_

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

#include "tools/version.hpp"

#include <string>
#include <vector>

// =================================================================================================
//      Global Options
// =================================================================================================

class GlobalOptions
{
public:

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    void add_to_app( CLI::App& app );

    void set_command_line_args( int const argc, char const* const* argv );

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

    void run_global( CLI::App const& app );

    void init();
    void print( CLI::App const& app ) const;

    std::string command_line() const;

    size_t verbosity() const;
    size_t threads() const;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    size_t verbosity_ = 1;
    size_t verbosity_cnt_ = 0;
    bool   verbosity_quiet_ = false;
    size_t threads_ = 0;

    std::vector<std::string> command_line_;

};

// =================================================================================================
//      Global Instance
// =================================================================================================

extern GlobalOptions global_options;

#endif // include guard
