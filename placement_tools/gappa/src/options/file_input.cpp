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

#include "options/file_input.hpp"

#include "options/global.hpp"

#include "genesis/utils/core/fs.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>

// Relative path resolution for printing.
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <linux/limits.h>

#ifdef GENESIS_OPENMP
#   include <omp.h>
#endif

// =================================================================================================
//      Setup Functions
// =================================================================================================

CLI::Option* FileInputOptions::add_multi_file_input_opt_to_app(
    CLI::App* sub,
    std::string const& type,
    std::string const& extension,
    bool               required,
    std::string const& group
){
    // Correct setup check.
    if( option_ != nullptr ) {
        throw std::domain_error( "Cannot use the same FileInputOptions object multiple times." );
    }

    // Store file type info.
    file_type_ = type;
    file_ext_  = extension;

    // Input files.
    option_ = sub->add_option(
        "--" + type + "-path",
        raw_paths_,
        "List of " + type + " files or directories to process. " +
        "For directories, only files with the extension ." + extension + " are processed."
    );
    if( required ) {
        option_->required();
    }

    // Check if it is a path.
    option_->check([]( std::string const& path ){
        if( ! genesis::utils::path_exists( path ) ) {
            return std::string( "Path is neither a file nor a directory: " + path );
        }
        return std::string();
    });
    option_->group( group );

    return option_;
}

// =================================================================================================
//      Run Functions
// =================================================================================================

size_t FileInputOptions::file_count() const
{
    return file_paths().size();
}

std::vector<std::string> const& FileInputOptions::file_paths() const
{
    #pragma omp critical(GAPPA_FILE_INPUT_PATHS)
    {
        if( resolved_paths_.empty() ) {
            using namespace genesis::utils;
            for( auto const& path : raw_paths_ ) {
                if( is_file( path ) ) {

                    resolved_paths_.push_back( path );

                } else if( is_dir( path ) ) {

                    // Get all files in dir.
                    auto list = dir_list_files( path, true, ".*\\." + file_ext_ + "$" );
                    for( auto const& jplace : list ) {
                        resolved_paths_.push_back( jplace );
                    }

                } else {
                    // throw std::runtime_error( "Not a valid file or directory: " + path );
                    throw CLI::ValidationError(
                        "--" + file_type_ + "-path", "Not a valid file or directory: " + path
                    );
                }
            }

            // If required, we actually need files!
            if( resolved_paths_.empty() && option_->get_required() ) {
                throw CLI::ValidationError(
                    "--" + file_type_ + "-path", "No files found."
                );
            }

            // We sort them to get reproducible order.
            std::sort( resolved_paths_.begin(), resolved_paths_.end() );
        }
    }

    return resolved_paths_;
}

std::string const& FileInputOptions::file_path( size_t index ) const
{
    auto const& files = file_paths();
    if( index >= files.size() ) {
        throw std::runtime_error( "Invalid file index." );
    }
    return files[ index ];
}

std::vector<std::string> const& FileInputOptions::raw_file_paths() const
{
    return raw_paths_;
}

std::vector<std::string> FileInputOptions::base_file_names() const
{
    using namespace genesis::utils;

    auto paths = file_paths();
    for( auto& path : paths ) {
        path = file_filename( file_basename( path ));
    }
    return paths;
}

std::string FileInputOptions::base_file_name( size_t index ) const
{
    using namespace genesis::utils;
    return file_filename( file_basename( file_path( index )));
}

void FileInputOptions::print() const
{
    std::string type = file_type_;
    if( ! type.empty() ) {
        type = " " + type;
    }

    // Print list of files, depending on verbosity.
    auto const& files = file_paths();
    if( global_options.verbosity() == 0 ) {
        return;
    } else if( global_options.verbosity() == 1 ) {
        std::cout << "Found " << files.size() << type << " files.\n";
    } else if( global_options.verbosity() == 2 ) {
        std::cout << "Found " << files.size() << type << " files: ";
        for( auto const& file : files ) {
            if( &file != &files[0] ) {
                std::cout << ", ";
            }
            std::cout << genesis::utils::file_basename( file );
        }
        std::cout << "\n";
    } else {
        std::cout << "Found " << files.size() << type << " files:\n";

        char resolved_path[PATH_MAX];
        for( auto const& file : files ) {
            auto ptr = realpath( file.c_str(), resolved_path );
            if( errno == 0 ) {
                std::cout << "  - " << ptr << "\n";
            } else {
                std::cout << "  - " << file << "\n";
                // use std::strerror(errno) to get error message
                errno = 0;
            }
        }
    }
}
