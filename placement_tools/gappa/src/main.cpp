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

#include "commands/analyze.hpp"
// #include "commands/edit.hpp"
#include "commands/prepare.hpp"

#include "options/global.hpp"
#include "tools/help.hpp"
#include "tools/version.hpp"

// =================================================================================================
//      Main Program
// =================================================================================================

int main( int argc, char** argv )
{
    // Activate logging.
    // utils::Logging::log_to_stdout();
    // utils::Logging::details.time = true;
    //
    // utils::Options::get().number_of_threads( 4 );
    // LOG_BOLD << utils::Options::get().info();

    CLI::App app{ gappa_header() };
    app.require_subcommand( 1 );
    app.fallthrough( true );
    app.set_name( "gappa" );

    // Add app-wide global options.
    global_options.set_command_line_args( argc, argv );
    global_options.add_to_app( app );

    // Set up all subcommands.
    setup_prepare( app );
    // setup_edit( app );
    setup_analyze( app );
    setup_wiki( app );

    try {
        app.parse( argc, argv );
    } catch ( CLI::ParseError const& e ) {
        return app.exit( e );
    }

    return 0;
}
