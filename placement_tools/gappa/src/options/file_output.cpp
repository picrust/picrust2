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

#include "options/file_output.hpp"

#include "genesis/utils/core/fs.hpp"

#include <algorithm>
#include <stdexcept>

// =================================================================================================
//      Setup Functions
// =================================================================================================

CLI::Option* FileOutputOptions::add_output_dir_opt_to_app( CLI::App* sub )
{
    return add_output_dir_opt_to_app( sub, "" );
}

CLI::Option* FileOutputOptions::add_output_dir_opt_to_app(
    CLI::App* sub,
    std::string const& name,
    std::string const& initial_value,
    std::string const& group
) {
    // Correct setup check.
    if( out_dir_option != nullptr ) {
        throw std::domain_error( "Cannot use the same FileOutputOptions object multiple times." );
    }

    // Setup.
    auto const optname = "--" + name + ( name.empty() ? "" : "-" ) + "out-dir";
    name_ = name;
    out_dir_ = initial_value;

    // Add option
    out_dir_option = sub->add_option(
        optname,
        out_dir_,
        "Directory to write " + name + ( name.empty() ? "" : " " ) + "files to",
        true
    );
    // out_dir_option->check( CLI::ExistingDirectory );
    out_dir_option->group( group );

    // TODO add function to overwrite files, which sets the genesis option for this. add this to global!

    return out_dir_option;
}

CLI::Option* FileOutputOptions::add_file_prefix_opt_to_app(
    CLI::App* sub,
    std::string const& name,
    std::string const& initial_value,
    std::string const& group
) {
    // Correct setup check.
    if( prefix_option != nullptr ) {
        throw std::domain_error( "Cannot use the same FileOutputOptions object multiple times." );
    }

    // Setup.
    auto const optname = "--" + name + ( name.empty() ? "" : "-" ) + "file-prefix";
    prefix_ = initial_value;

    // Add option
    prefix_option = sub->add_option(
        optname,
        prefix_,
        "File prefix for " + ( name.empty() ? "output" : name ) + " files",
        true
    );
    prefix_option->check([]( std::string const& prefix ){
        if( ! genesis::utils::is_valid_filname( prefix ) ) {
            return std::string(
                "File prefix contains invalid characters (<>:\"\\/|?*) or surrounding whitespace."
            );
        }
        return std::string();
    });
    prefix_option->group( group );

    return prefix_option;
}

// =================================================================================================
//      Run Functions
// =================================================================================================

std::string FileOutputOptions::out_dir() const
{
    // Create dir if needed. This might create the dir also in cases were something failes later,
    // so we end up with an empty dir. This is however common in many other programs as well,
    // so let's not bother with this.
    genesis::utils::dir_create( out_dir_, true );

    return genesis::utils::dir_normalize_path( out_dir_ );
}

std::string FileOutputOptions::file_prefix() const
{
    return prefix_;
}

void FileOutputOptions::check_nonexistent_output_files(
    std::vector<std::string> const& filenames
) const {
    using namespace genesis::utils;

    // Shortcut: if the dir is not created yet, there cannot be any existing files in it.
    // We do this check here, so that we can be sure later in this function that the dir
    // is there, so that listing it contents etc actually works.
    if( ! genesis::utils::dir_exists( out_dir_ ) ) {
        return;
    }

    // Get basic strings
    auto const optname = "--" + name_ + ( name_.empty() ? "" : "-" ) + "out-dir";

    // Check if any of the files exists. Old version without regex.
    // std::string const dir = dir_normalize_path( out_dir_ );
    // for( auto const& file : filenames ) {
    //     if( file_exists( dir + file ) ) {
    //         throw CLI::ValidationError(
    //             "--out-dir (" + out_dir_ +  ")", "Output file already exists: " + file
    //         );
    //     }
    // }

    // TODO if file overwrite option is added, this check should become a warning!

    // TODO using regexes gives weird user output if any of the files already exists.
    // for example, dispersion_imbalances_sd\.svg
    // better use normal wildcarts, and use some string replacement here before calling
    // dir list contents.

    // Check if any of the files exists.
    for( auto const& file : filenames ) {
        auto const dir_cont = dir_list_contents( out_dir_, true, file );
        if( ! dir_cont.empty() ) {
            throw CLI::ValidationError(
                optname + " (" + out_dir_ +  ")", "Output path already exists: " + file
            );
        }
    }

    // Check if any file name is duplicated.
    // If so, we will have a problem after the first file has been written.
    auto cpy = filenames;
    std::sort( cpy.begin(), cpy.end() );
    auto const adj = std::adjacent_find( cpy.begin(), cpy.end() ) ;
    if( adj != cpy.end() ) {
        throw CLI::ValidationError(
            optname, "Output file name used multiple times: " + ( *adj )
        );
    }

    // TODO there is a change that multiple named output dirs are set to the same real dir,
    // and that then files with the same names are written. well, this would be comples to
    // check, so not now...
}
