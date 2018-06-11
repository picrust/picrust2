#ifndef GAPPA_OPTIONS_SEQUENCE_INPUT_H_
#define GAPPA_OPTIONS_SEQUENCE_INPUT_H_

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

#include "options/file_input.hpp"

#include "genesis/sequence/formats/fasta_reader.hpp"
#include "genesis/sequence/formats/phylip_reader.hpp"

#include <string>
#include <vector>

// =================================================================================================
//      Sequence Input Options
// =================================================================================================

/**
 * @brief Helper class for adding input options for sequence data files to a command.
 *
 * The class offers options to read fasta and/or phylip files. Depending on the file type, certain
 * reading options are (not) available. For example, sequential reading of the sequences in a file
 * is currently only supported for fasta files, but not for phylip files.
 */
class SequenceInputOptions
    : public FileInputOptions
{
public:

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    SequenceInputOptions();
    virtual ~SequenceInputOptions() = default;

    SequenceInputOptions( SequenceInputOptions const& other ) = default;
    SequenceInputOptions( SequenceInputOptions&& )            = default;

    SequenceInputOptions& operator= ( SequenceInputOptions const& other ) = default;
    SequenceInputOptions& operator= ( SequenceInputOptions&& )            = default;

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    /**
     * @brief Add options that allow reading any sequence file (fasta or phylip).
     */
    CLI::Option* add_sequence_input_opt_to_app( CLI::App* sub, bool required = true );

    /**
     * @brief Add options that allow to read fasta files only.
     */
    CLI::Option* add_fasta_input_opt_to_app( CLI::App* sub, bool required = true );

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

    /**
     * @brief Read in the sequence file at @p index in the list of input files and return it.
     *
     * The function tries both fasta and phylip format, using the extension as a hint, and if this
     * fails, tries the respective other format.
     *
     * See FileInputOptions::file_count() for the number of input files (valid range for the index)
     * and FileInputOptions::file_paths() for their list.
     */
    genesis::sequence::SequenceSet sequence_set( size_t index ) const;

    /**
     * @brief Read all sequences of the provided files into one set.
     */
    genesis::sequence::SequenceSet sequence_set_all() const;

    /**
     * @brief Get the FastaReader used when processing fasta files.
     *
     * Use this to change the behaviour of the reader if needed for a command.
     */
    genesis::sequence::FastaReader& fasta_reader()
    {
        return fasta_reader_;
    }

    genesis::sequence::FastaReader const& fasta_reader() const
    {
        return fasta_reader_;
    }

    /**
     * @brief Get the PhylipReader used when processing fasta files.
     *
     * Use this to change the behaviour of the reader if needed for a command.
     */
    genesis::sequence::PhylipReader& phylip_reader()
    {
        return phylip_reader_;
    }

    genesis::sequence::PhylipReader const& phylip_reader() const
    {
        return phylip_reader_;
    }

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    genesis::sequence::FastaReader  fasta_reader_;
    genesis::sequence::PhylipReader phylip_reader_;

    static const std::string fasta_extensions_;
    static const std::string phylip_extensions_;

};

#endif // include guard
