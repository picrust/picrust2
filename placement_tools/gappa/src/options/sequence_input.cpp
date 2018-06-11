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

#include "options/sequence_input.hpp"

#include "genesis/sequence/sequence_set.hpp"
#include "genesis/utils/core/fs.hpp"

#include <iostream>
#include <stdexcept>

// =================================================================================================
//      Constructor and Constants
// =================================================================================================

// Fasta extensions: https://en.wikipedia.org/wiki/FASTA_format#File_extension
const std::string SequenceInputOptions::fasta_extensions_ = "fasta|fas|fsa|fna|ffn|faa|frn";
const std::string SequenceInputOptions::phylip_extensions_ = "phylip|phy";

SequenceInputOptions::SequenceInputOptions()
{
    fasta_reader_.site_casing( genesis::sequence::FastaReader::SiteCasing::kUnchanged );
    phylip_reader_.site_casing( genesis::sequence::PhylipReader::SiteCasing::kUnchanged );
    phylip_reader_.mode( genesis::sequence::PhylipReader::Mode::kAutomatic );
}

// =================================================================================================
//      Setup Functions
// =================================================================================================

CLI::Option* SequenceInputOptions::add_sequence_input_opt_to_app( CLI::App* sub, bool required )
{
    return FileInputOptions::add_multi_file_input_opt_to_app(
        sub, "sequence", "(" + fasta_extensions_ + "|" + phylip_extensions_ + ")", required
    )->group( "Input" );
}

CLI::Option* SequenceInputOptions::add_fasta_input_opt_to_app( CLI::App* sub, bool required )
{
    return FileInputOptions::add_multi_file_input_opt_to_app(
        sub, "fasta", "(" + fasta_extensions_ + ")", required
    )->group( "Input" );
}

// =================================================================================================
//      Run Functions
// =================================================================================================

genesis::sequence::SequenceSet SequenceInputOptions::sequence_set( size_t index ) const
{
    using namespace genesis::sequence;

    // Prepare.
    SequenceSet result;
    auto const& file_name = file_path( index );
    auto const ext = genesis::utils::file_extension( file_name );

    // Store the exception message thrown on the first attempt.
    // Then, if all fails, throw an exeption with the message again.
    // For example, try phylip, fail, store message, try fasta, fail too, throw message again.
    // Thus, in the end, we throw the error message that fits with the file extension,
    // in the hope that this is most useful to the user.
    std::string error_message;

    if( ext == "phylip" || ext == "phy" ) {

        // Try phylip if extension says so. Return if successfull.
        try{
            phylip_reader_.from_file( file_name, result );
            return result;
        } catch ( std::exception& ex ) {
            error_message = ex.what();
        }

        // Otherwise try fasta, again returning on success.
        try{
            result.clear();
            fasta_reader_.from_file( file_name, result );
            return result;
        } catch ( ... ) {}

    } else {

        // For all other extensions, try fasta first. Return if successfull.
        try{
            fasta_reader_.from_file( file_name, result );
            return result;
        } catch ( std::exception& ex ) {
            error_message = ex.what();
        }

        // If this does not work, try phylip.
        try{
            result.clear();
            phylip_reader_.from_file( file_name, result );
            return result;
        } catch ( ... ) {}
    }

    // If we are here, none of the above worked.
    throw std::runtime_error(
        "Input file " + file_name + " cannot be read as either fasta or phylip. "
        "Error message: " + error_message
    );
}

genesis::sequence::SequenceSet SequenceInputOptions::sequence_set_all() const
{
    using namespace genesis::sequence;
    SequenceSet result;
    for( size_t i = 0; i < file_count(); ++i ) {
        auto tmp = sequence_set( i );
        for( auto& seq : tmp ) {
            result.add( std::move( seq ));
        }
    }
    return result;
}
