#ifndef GAPPA_OPTIONS_METADATA_INPUT_H_
#define GAPPA_OPTIONS_METADATA_INPUT_H_

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

#include "genesis/utils/containers/dataframe.hpp"

#include <string>
#include <vector>

// =================================================================================================
//      Metadata Input Options
// =================================================================================================

/**
 * @brief Helper class to add command line parameter to input a metadata table.
 */
class MetadataInputOptions
{
public:

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    MetadataInputOptions()  = default;
    ~MetadataInputOptions() = default;

    MetadataInputOptions( MetadataInputOptions const& other ) = default;
    MetadataInputOptions( MetadataInputOptions&& )            = default;

    MetadataInputOptions& operator= ( MetadataInputOptions const& other ) = default;
    MetadataInputOptions& operator= ( MetadataInputOptions&& )            = default;

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    /**
     * @brief Add the option to an App.
     */
    CLI::Option* add_metadata_input_opt_to_app(
        CLI::App* sub
    );

    /**
     * @brief Return the CLI11 option for the metadata file.
     */
    CLI::Option* file_option()
    {
        return file_option_;
    }

    /**
     * @brief Return the CLI11 option for the metadata file.
     */
    CLI::Option const* file_option() const
    {
        return file_option_;
    }

    /**
     * @brief Return the CLI11 option for the metadata fields.
     */
    CLI::Option* fields_option()
    {
        return fields_option_;
    }

    /**
     * @brief Return the CLI11 option for the metadata fields.
     */
    CLI::Option const* fields_option() const
    {
        return fields_option_;
    }

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

public:

    /**
     * @brief Get the metadata file path.
     */
    std::string const& metadata_file() const
    {
        return metadata_file_;
    }

    /**
     * @brief Get metadata table.
     */
    genesis::utils::Dataframe<double> read_metadata() const;

    /**
     * @brief Return whether the row names of a dataframe are the same as the given list of names,
     * order independently.
     */
    static bool check_row_names(
        genesis::utils::Dataframe<double> const& df,
        std::vector<std::string> const& row_names
    );

    /**
     * @brief Sort the rows of a dataframe by a given order.
     */
    static genesis::utils::Dataframe<double> sort_rows(
        genesis::utils::Dataframe<double> const& df,
        std::vector<std::string> const&          row_name_order
    );

    /**
    * @brief Print some user output related to the option.
    */
    void print() const;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    std::string metadata_file_;
    std::string metadata_fields_;

    CLI::Option* file_option_ = nullptr;
    CLI::Option* fields_option_ = nullptr;

};

#endif // include guard
