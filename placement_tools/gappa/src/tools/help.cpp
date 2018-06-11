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

#include "tools/help.hpp"

#include "genesis/utils/core/fs.hpp"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// =================================================================================================
//      Wiki Pages
// =================================================================================================

// -------------------------------------------------------------------------
//     Setup
// -------------------------------------------------------------------------

void setup_wiki( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto options = std::make_shared<WikiOptions>();
    auto sub = app.add_subcommand(
        "wiki",
        "Create wiki help pages for gappa."
    )->group( "" );

    // Need to capture the main app, as the wiki needs this to run.
    options->app = &app;
    while( options->app->get_parent() ) {
        options->app = options->app->get_parent();
    }

    // Markdown dir option.
    auto md_dir_opt = sub->add_option(
        "--md-dir",
        options->md_dir,
        "Directory with the Markdown files documenting the gappa commands."
    );
    md_dir_opt->group( "Settings" );
    md_dir_opt->check( CLI::ExistingDirectory );
    // md_dir_opt->required();

    // Out dir option.
    auto out_dir_opt = sub->add_option(
        "--out-dir",
        options->out_dir,
        "Directory to write Wiki files to. Should be a git clone of the wiki repository."
    );
    out_dir_opt->group( "Settings" );
    out_dir_opt->check( CLI::ExistingDirectory );
    // out_dir_opt->required();

    // Set the run function as callback to be called when this subcommand is issued.
    // Hand over the options by copy, so that their shared ptr stays alive in the lambda.
    sub->set_callback( [ options ]() {
        run_wiki( *options );
    });
}

// -------------------------------------------------------------------------
//     App Subcommand Helpers
// -------------------------------------------------------------------------

/**
 * @brief Get the immediate subcommands of an App, sorted by name.
 */
std::vector<CLI::App*> get_sorted_subcommands( CLI::App const* app )
{
    auto subcomms = app->get_subcommands( false );
    std::sort(
        subcomms.begin(), subcomms.end(),
        []( CLI::App* lhs, CLI::App* rhs ){
            return lhs->get_name() < rhs->get_name();
        }
    );
    return subcomms;
}

/**
 * @brief Get all subcommands, recursively, of an App, sorted by name.
 */
std::vector<CLI::App const*> get_all_subcommands( CLI::App const* app )
{
    std::vector<CLI::App const*> result;

    // Fill with subcommands.
    auto subcomms = get_sorted_subcommands( app );
    for( auto const& sc : subcomms ) {
        result.push_back( sc );

        // Recurse. The copying is wasteful, but it's not that much, and I am lazy today.
        auto subsubcomms = get_all_subcommands( sc );
        for( auto const& ssc : subsubcomms ) {
            result.push_back( ssc );
        }
    }

    return result;
}

/**
 * @brief Add the contents of a file to a stream.
 */
void add_markdown_content( WikiOptions const& options, std::string const& md_file, std::ostream& os )
{
    using namespace genesis::utils;

    // Add markdown file content.
    std::string const fn = dir_normalize_path( options.md_dir ) + md_file + ".md";
    if( file_exists( fn ) ) {
        std::ifstream mds( fn );
        os << mds.rdbuf();
    } else {
        std::cout << " - No documentation markdown found: " << md_file << "\n";
    }
}

// -------------------------------------------------------------------------
//     Make Options Table
// -------------------------------------------------------------------------

void make_options_table( std::vector<CLI::Option*> options, std::ostream& os )
{
    // map from group name to table contents.
    // we use a vec to keep order.
    // std::map<std::string, std::string> opt_helps;
    std::vector<std::pair<std::string, std::string>> opt_helps;

    // Add lines for each group.
    for( auto const& opt : options ) {

        // Do not add help option.
        if( opt->get_name() == "-h,--help" ) {
            continue;
        }

        // Write to temporary stream.
        std::stringstream tmp_os;

        // Simple version that makes a markdown table.
        // tmp_os << "| <nobr>`" << opt->get_name() << "`</nobr> ";
        // tmp_os << "|";
        // if( opt->get_required() ) {
        //     tmp_os << " **Required.**";
        // }
        // if( ! opt->help_aftername().empty() ) {
        //     // print stuff without leading space.
        //     tmp_os << " `" << opt->help_aftername().substr( 1 ) << "`<br>";
        // }
        //
        // auto descr = opt->get_description();
        // tmp_os << " " << descr << " |\n";
        // // tmp_os << " " << opt->get_description() << " |\n";
        // // tmp_os << "| " << opt->get_description() << " |\n";

        // Add content to the group help.
        tmp_os << "<tr><td><code>" << opt->get_name() << "</code></td>";
        tmp_os << "<td>";
        if( opt->get_required() ) {
            tmp_os << "<strong>Required.</strong>";
        }
        if( ! opt->help_aftername().empty() ) {
            // print stuff without leading space.
            auto han = opt->help_aftername().substr( 1 );
            auto const rp = han.find( " (REQUIRED)" );
            if( rp != std::string::npos ) {
                han.erase( rp,  11 );
            }

            tmp_os << " <code>" << han << "</code><br>";
        }
        auto descr = opt->get_description();
        if( descr.substr( 0, 10 ) == "Required. " ) {
            descr = descr.substr( 10 );
        }
        tmp_os << " " << descr << "</td></tr>\n";
        // tmp_os << " " << opt->get_description() << " |\n";
        // tmp_os << "| " << opt->get_description() << " |\n";

        // Add content to the group help.
        // first check if the group was already used, and if not add it.
        auto get_group_content = [&]( std::string const& name ) -> std::string& {
            for( auto& elem : opt_helps ) {
                if( elem.first == name ) {
                    return elem.second;
                }
            }
            opt_helps.push_back({ name, "" });
            return opt_helps.back().second;
        };
        get_group_content( opt->get_group() ) += tmp_os.str();
        // opt_helps[ opt->get_group() ] += tmp_os.str();
    }

    // Simple markdown verison to print the groups and their tables.
    // for( auto const& gr : opt_helps ) {
    //     os << "**" << gr.first << ":**\n\n";
    //     os << "| Option  | Description |\n";
    //     os << "| ------- | ----------- |\n";
    //     // os << "| Option  | Type | Description |\n";
    //     // os << "| ------- | ---- | ----------- |\n";
    //     os << gr.second << "\n";
    // }

    // Print the groups and their tables
    os << "<table>\n";
    bool done_first_group = false;
    for( auto const& gr : opt_helps ) {
        if( done_first_group ) {
            os << "<tr height=30px></tr>\n";
        }
        os << "<thead><tr><th colspan=\"2\" align=\"left\">" << gr.first << "</th></tr></thead>\n";
        os << "<tbody>\n";
        os << gr.second;
        os << "</tbody>\n";
        done_first_group = true;
    }
    os << "</table>\n\n";
}

// -------------------------------------------------------------------------
//     Make Subcommands Table
// -------------------------------------------------------------------------

void make_subcommands_table( std::vector<CLI::App*> subcomms, std::ostream& os )
{
    os << "| Subcommand  | Description |\n";
    os << "| ----------- | ----------- |\n";

    for( auto const& subcomm : subcomms ) {
        os << "| [" << subcomm->get_name() << "](../wiki/Subcommand:-" << subcomm->get_name() << ") ";
        os << "| " << subcomm->get_description() << " |\n";
    }
    os << "\n";
}

// -------------------------------------------------------------------------
//     Make Wiki Page
// -------------------------------------------------------------------------

void make_wiki_command_page( WikiOptions const& options, CLI::App const& command )
{
    using namespace genesis::utils;

    // User output.
    std::cout << "Subcommand: " << command.get_name() << "\n";

    // Get stuff of this command.
    auto const subcomms = command.get_subcommands( false );
    auto const opts = command.get_options();

    // Open out file stream.
    std::string const out_file
        = dir_normalize_path( options.out_dir )
        + "Subcommand:-" + command.get_name() + ".md"
    ;
    if( ! file_exists( out_file )) {
        std::cout << " - No existing wiki file!\n";
    }
    std::ofstream os( out_file );

    // Get the usage line.
    std::string usage = command.get_name();
    auto parent = const_cast< CLI::App& >( command ).get_parent();
    while( parent ) {
        usage  = parent->get_name() + " " + usage;
        parent = parent->get_parent();
    }

    // We do not count the help option, so we need to manually check if there are any others.
    bool has_options = false;
    for( auto const& opt : opts ) {
        if( opt->get_name() != "-h,--help" ) {
            has_options = true;
            break;
        }
    }

    // Write command header.
    os << command.get_description() << "\n\n";
    os << "Usage: `" << usage;
    if( has_options ) {
        os << " [options]";
    }
    if( ! subcomms.empty() ) {
        if( command.get_require_subcommand_min() > 0 ) {
            os << " subcommand";
        } else {
            os << " [subcommand]";
        }
    }
    os << "`\n\n";

    // Print the options of thus command.
    if( has_options ) {
        os << "## Options\n\n";
        make_options_table( opts, os );
    }

    // Print the subcommands of this command.
    if( ! subcomms.empty() ) {
        os << "## Subcommands\n\n";
        make_subcommands_table( subcomms, os );
    }

    // Add markdown file content and finish.
    add_markdown_content( options, command.get_name(), os );
    os.close();
}

// -------------------------------------------------------------------------
//     Make Wiki Home Page
// -------------------------------------------------------------------------

void make_wiki_home_page( WikiOptions const& options )
{
    using namespace genesis::utils;

    // Make Home page.
    std::cout << "Home\n";

    // Open stream
    std::string const out_file = dir_normalize_path( options.out_dir ) + "Home.md";
    if( ! file_exists( out_file )) {
        std::cout << " - No existing wiki file!\n";
    }
    std::ofstream os( out_file );

    // Add home header.
    add_markdown_content( options, "Home_header", os );
    os << "\n";

    // Add submodule lists.
    auto subcomms = get_sorted_subcommands( options.app );
    for( auto const& sc : subcomms ) {
        if( sc->get_group() == "" ) {
            continue;
        }

        os << "### Module `" << sc->get_name() << "`\n\n";
        os << sc->get_description() << "\n\n";
        make_subcommands_table( get_sorted_subcommands( sc ), os );
    }

    // Add home footer.
    add_markdown_content( options, "Home_footer", os );
    os.close();
}

// -------------------------------------------------------------------------
//     Run Function
// -------------------------------------------------------------------------

void run_wiki( WikiOptions const& options )
{
    // Get all subcommands of the main app and make wiki pages for all subcommands.
    // auto subcomms = get_all_subcommands( options.app );
    // for( auto const& subcomm : subcomms ) {
    //     make_wiki_command_page( options, *subcomm );
    // }

    // Home page.
    make_wiki_home_page( options );

    // Now, make pages for the commands of the modules.
    auto subcomms = get_sorted_subcommands( options.app );
    for( auto const& sc : subcomms ) {
        for( auto const& ssc : get_sorted_subcommands( sc )) {
            make_wiki_command_page( options, *ssc );
        }
    }
}
