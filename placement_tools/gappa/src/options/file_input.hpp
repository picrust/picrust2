#ifndef GAPPA_OPTIONS_FILE_INPUT_H_
#define GAPPA_OPTIONS_FILE_INPUT_H_

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

#include <string>
#include <vector>

// =================================================================================================
//      File Input Options
// =================================================================================================

/**
 * @brief Helper class to add command line parameter to input mutliple files in a uniform way.
 *
 * This class is meant to be stored as a member of the actual options collection,
 * or to be derived from for special file type input options classes.
 * It offers a function to add mutliple file input to a command line interface,
 * and has several convenience functions to use the files that were passes by a user.
 */
class FileInputOptions
{
public:

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    FileInputOptions()  = default;
    virtual ~FileInputOptions() = default;

    FileInputOptions( FileInputOptions const& other ) = default;
    FileInputOptions( FileInputOptions&& )            = default;

    FileInputOptions& operator= ( FileInputOptions const& other ) = default;
    FileInputOptions& operator= ( FileInputOptions&& )            = default;

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    /**
     * @brief Add the options to an App.
     *
     * Takes a file type used for help messages, and an extension for valid files.
     * The extenion can be a regex, e.g., `(fas|fasta)`.
     */
    CLI::Option* add_multi_file_input_opt_to_app(
        CLI::App* sub,
        std::string const& type,
        std::string const& extension,
        bool               required = true,
        std::string const& group = "Input"
    );

    /**
     * @brief Return the CLI11 option that this object belongs to.
     */
    CLI::Option* option()
    {
        return option_;
    }

    /**
     * @brief Return the CLI11 option that this object belongs to.
     */
    CLI::Option const* option() const
    {
        return option_;
    }

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

public:

    /**
     * @brief Get the number of files that were provided by the user.
     *
     * Returns `file_paths().size()`.
     */
    size_t file_count() const;

    /**
     * @brief Get the resolved full file paths of all files provided by the user.
     *
     * This function uses the list of paths given by the user.
     * For each of them, it checks whether it is a file or a directory.
     * Files are immediately added to the result list, while directories are scanned
     * for any files with the extension in them, which are then added to the result list.
     *
     * This allows the users to hand over either their own files (no matter their extension), or
     * whole directory paths for convenience, which then however only use files ending in the
     * extension.
     */
    std::vector<std::string> const& file_paths() const;

    /**
     * @brief Get a specific file from the list.
     *
     * Returns `file_paths().at( index )`.
     */
    std::string const& file_path( size_t index ) const;

    /**
     * @brief Return the list of paths as provided by the user, that is, without processing.
     */
    std::vector<std::string> const& raw_file_paths() const;

    /**
     * @brief Get the file names of the provided files, i.e., without directory and ending.
     *
     * This function calls genesis::utils::file_basename() and genesis::utils::file_filename() for
     * all paths. The result is for example useful for user output.
     */
    std::vector<std::string> base_file_names() const;

    /**
     * @brief Get the file name of the file at the index, i.e., without directory and ending.
     *
     * This function is the same as base_file_names(), just for one file.
     */
    std::string base_file_name( size_t index ) const;

    /**
     * @brief Print some user output related to these options.
     */
    virtual void print() const;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    std::vector<std::string> raw_paths_;
    mutable std::vector<std::string> resolved_paths_;

    std::string file_type_;
    std::string file_ext_;

    CLI::Option* option_ = nullptr;

};

#endif // include guard
