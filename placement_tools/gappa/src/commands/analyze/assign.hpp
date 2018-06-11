#ifndef GAPPA_COMMANDS_ANALYZE_ASSIGN_H_
#define GAPPA_COMMANDS_ANALYZE_ASSIGN_H_

/*
    gappa - Genesis Applications for Phylogenetic Placement Analysis
    Copyright (C) 2017-2018 Pierre Barbera, Lucas Czech and HITS gGmbH

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

#include "options/jplace_input.hpp"
#include "options/file_output.hpp"

#include "genesis/taxonomy/taxon_data.hpp"

#include <string>
#include <vector>

// =================================================================================================
//      Options
// =================================================================================================

class AssignOptions
{
public:

    std::string         taxon_file;
    std::string         sub_taxopath;
    JplaceInputOptions  jplace_input;

    double  dist_ratio = -1.0;

    FileOutputOptions   output_dir;
};

// =================================================================================================
//      Functions
// =================================================================================================

void setup_assign( CLI::App& app );
void run_assign( AssignOptions const& options );

// =================================================================================================
// Weird Classes
// =================================================================================================

class AssignTaxonData : public genesis::taxonomy::BaseTaxonData
{
    // -------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------

public:

    virtual ~AssignTaxonData() = default;

    // Move ctor and assignment.
    AssignTaxonData( AssignTaxonData&& )             = delete;
    AssignTaxonData& operator= ( AssignTaxonData&& ) = delete;

protected:

    AssignTaxonData() = default;

    // Copy ctor and assignment.
    AssignTaxonData( AssignTaxonData const& )             = default;
    AssignTaxonData& operator= ( AssignTaxonData const& ) = default;

public:

    static std::unique_ptr< AssignTaxonData > create()
    {
        return std::unique_ptr< AssignTaxonData >( new AssignTaxonData() );
    }

    virtual std::unique_ptr< BaseTaxonData > clone() const override
    {
        return std::unique_ptr< AssignTaxonData >( new AssignTaxonData( *this ) );
    }

    // -----------------------------------------------------
    //     Data Members
    // -----------------------------------------------------

    double aLWR = 0.0;
    double LWR  = 0.0;
};

#endif // include guard