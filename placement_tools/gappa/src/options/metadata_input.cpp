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

#include "options/metadata_input.hpp"

#include "options/global.hpp"

#include "genesis/utils/containers/dataframe.hpp"
#include "genesis/utils/containers/dataframe/reader.hpp"
#include "genesis/utils/core/algorithm.hpp"
#include "genesis/utils/core/fs.hpp"
#include "genesis/utils/formats/csv/reader.hpp"
#include "genesis/utils/text/string.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <unordered_set>

// =================================================================================================
//      Setup Functions
// =================================================================================================

CLI::Option* MetadataInputOptions::add_metadata_input_opt_to_app(
    CLI::App* sub
){
    // Correct setup check.
    if( file_option_ != nullptr || fields_option_ != nullptr ) {
        throw std::domain_error( "Cannot use the same MetadataInputOptions object multiple times." );
    }

    // Metadata file
    file_option_= sub->add_option(
        "--metadata-file",
        metadata_file_,
        "Csv file with the metadata columns to use."
    );
    file_option_->check( CLI::ExistingFile );
    file_option_->group( "Input" );

    // Metadata fileds
    fields_option_ = sub->add_option(
        "--metadata-fields",
        metadata_fields_,
        "Metadata fields to use, separated by commata. If not provided, all are used."
    )->group( "Input" );

    return file_option_;
}

// =================================================================================================
//      Run Functions
// =================================================================================================

genesis::utils::Dataframe<double> MetadataInputOptions::read_metadata() const
{
    using namespace genesis::utils;

    // Prepare a reader that can convert anything to double.
    // We filter out later.
    DataframeReader<double> reader;
    reader.parse_value_functor( []( std::string const& cell ){
        double v;
        try{
            v = std::stod( cell );
        } catch( ... ) {
            v = std::numeric_limits<double>::quiet_NaN();
        }
        return v;
    });

    // Do the reading.
    auto df = reader.from_file( metadata_file_ );

    // If the user specified to use only certain columns, remove the others.
    if( ! metadata_fields_.empty() ) {

        // We use a csv reader to make sure that we properly parse the fields option.
        auto const fields_arr = CsvReader().from_string( metadata_fields_ );
        if( fields_arr.size() != 1 ) {
            throw CLI::ValidationError(
                "--metadata-fields (" + metadata_fields_ +  ")",
                "Invalid metadata fields list. Needs to be a comma-separated list of fields."
            );
        }

        // Turn it into a convenient set.
        std::unordered_set<std::string> fields;
        for( auto const& e : fields_arr[0] ) {
            if( ! e.empty() ) {
                fields.emplace( e );
            }
        }

        // Remove the columns that are not in the fields from the dataframe,
        // and remove the ones that are from the fields set.
        for( auto const& cn : df.col_names() ) {
            if( fields.count( cn ) > 0 ) {
                fields.erase( cn );
            } else {
                df.remove_col( cn );
            }
        }

        // Now, warn if there are remaining fields, that is, ones that were not found in the dataframe.
        if( fields.size() > 0 ) {
            std::cout << "Warning: You specified to use only certain metadata fields to be used ";
            std::cout << "via the --metadata-fields option. However, the following specified fields ";
            std::cout << "were not found in the metadata file: " << join( fields, ", " ) << "\n";
        }
    }

    // Now check for any "empty" columns that just contains zeros or invalid values.
    // Those can result from metadata columns that are not numbers,
    // and cannot be used for our methods. So, remove them.
    std::vector<std::string> rem_names;
    size_t i = 0;
    while( i < df.cols() ) {
        auto const& col = df[i];
        auto const bad_vals = std::all_of( col.begin(), col.end(), []( double v ){
            return v == 0.0 || ! std::isfinite( v );
        });
        if( bad_vals ) {
            rem_names.push_back( col.name() );
            df.remove_col( i );
        } else {
            ++i;
        }
    }

    // Some user warning if we removed columns.
    if( rem_names.size() > 0 ) {
        std::cout << "Warning: The following columns of the metadata file contained non-numerical ";
        std::cout << "data or only invalid values, which cannot be used here, and are hence ignored: ";
        std::cout << join( rem_names, ", " ) << "\n";
    }

    // User output
    if( global_options.verbosity() == 1 ) {
        std::cout << "Using " << df.col_names().size() << " metadata field";
        std::cout << ( df.col_names().size() == 1 ? "" : "s" ) << ".\n";
    }
    if( global_options.verbosity() == 2 ) {
        std::cout << "Using metadata fields: " << join( df.col_names(), ", " ) << "\n";
    }
    if( global_options.verbosity() >= 3 ) {
        std::cout << "Using metadata fields: \n";
        for( size_t i = 0; i < df.cols(); ++i ) {
            size_t const cnt = std::count_if( df[i].begin(), df[i].end(), []( double v ){
                return std::isfinite( v );
            });
            if( cnt != df.rows() ) {
                std::cout << " - " << df[i].name() << " (" << cnt << " of " << df.rows() << " valid values)\n";
            } else {
                std::cout << " - " << df[i].name() << "\n";
            }
        }
    }

    // Last check after removing all unused columns and printing the remaining onces.
    auto col_names = df.col_names();
    std::sort( col_names.begin(), col_names.end() );
    if( std::adjacent_find( col_names.begin(), col_names.end() ) != col_names.end() ) {
        throw std::runtime_error(
            "The metadata file contains duplicate fields (column names)."
        );
    }

    return df;
}

bool MetadataInputOptions::check_row_names(
    genesis::utils::Dataframe<double> const& df,
    std::vector<std::string> const&          row_names
) {

    // Helper function to sort a vector.
    auto sort_vec = []( std::vector<std::string> vec ){
        // std::for_each( vec.begin(), vec.end(), []( std::string& s ){ s = to_lower(s); });
        std::sort( vec.begin(), vec.end() );
        return vec;
    };

    // Check if the filenames match the metadata data rows. Compare vecs order-independently.
    auto const df_sci = sort_vec( df.row_names() );
    auto const bn_sci = sort_vec( row_names );

    return df_sci == bn_sci;
}

genesis::utils::Dataframe<double> MetadataInputOptions::sort_rows(
    genesis::utils::Dataframe<double> const& df,
    std::vector<std::string> const&          row_name_order
) {
    // We here simply make a sorted copy of the dataframe, because sorting inline is nasty.
    // This is not really a nice solution, but works for now.

    // Make a dataframe with the correct columns.
    genesis::utils::Dataframe<double> res;
    for( auto const& col : df ) {
        res.add_col( col.name() );
    }

    // Add the rows in the correct order, and fill in the values.
    for( auto const& row_name : row_name_order ) {
        res.add_row( row_name );
        auto const ridx = res.row_index( row_name );
        assert( ridx == res.rows() - 1 );

        auto const old_ridx = df.row_index( row_name );

        for( size_t cidx = 0; cidx < res.cols(); ++cidx ) {
            res( ridx, cidx ) = df( old_ridx, cidx );
        }
    }

    return res;
}

void MetadataInputOptions::print() const
{
    // User output
    if( global_options.verbosity() >= 2 ) {
        std::cout << "Metadata file: " << metadata_file_ << "\n";
    }
}
