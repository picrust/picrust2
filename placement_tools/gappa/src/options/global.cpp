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

#include "options/global.hpp"

#include "tools/version.hpp"
#include "tools/help.hpp"

#include "genesis/utils/core/options.hpp"

#include <thread>

// =================================================================================================
//      Setup Functions
// =================================================================================================

void GlobalOptions::add_to_app( CLI::App& app )
{
    // Verbosity
    auto opt_verb_l = app.add_option(
        "--verbosity",
        verbosity_,
        "Verbosity level [0-3]",
        true
    );
    auto opt_verb_s = app.add_flag(
        "-v",
        verbosity_cnt_,
        "Verbosity; add multiple times for more (-vvv)"
    );
    auto opt_verb_q = app.add_flag(
        "--quiet",
        verbosity_quiet_,
        "Set verbosity to 0, that is, only report errors and warnings"
    );
    opt_verb_l->excludes( opt_verb_s );
    opt_verb_l->excludes( opt_verb_q );
    opt_verb_s->excludes( opt_verb_l );
    opt_verb_s->excludes( opt_verb_q );
    opt_verb_q->excludes( opt_verb_l );
    opt_verb_q->excludes( opt_verb_s );

    // Threads
    app.add_option(
        "--threads",
        threads_,
        "Number of threads to use for calculations"
    );

    // TODO add random seed option
    // TODO add global file overwrite option.
    // TODO in order to run callbacks for certain options, use ther full functional form!
    // for example, the allow overwrite option (yet to do), or threads or the like can use this.
    // then, init is no longer needed

    // Run the app wide callback
    app.set_callback([ this, &app ](){
        run_global( app );
    });

    // Footer
    app.set_footer( gappa_title() );
}

void GlobalOptions::set_command_line_args( int const argc, char const* const* argv )
{
    // Store all arguments in the array.
    command_line_.clear();
    for (int i = 0; i < argc; i++) {
        command_line_.push_back(argv[i]);
    }
}

// =================================================================================================
//      Run Functions
// =================================================================================================

void GlobalOptions::run_global( CLI::App const& app )
{
    init();
    print( app );
}

void GlobalOptions::init()
{
    // If user did not provide number, use hardware value.
    if( threads_ == 0 ) {
        threads_ = std::thread::hardware_concurrency();
    }

    // If hardware value is not available, just use 1 thread.
    // This is executed if the call to hardware_concurrency fails.
    if( threads_ == 0 ) {
        threads_ = 1;
    }

    // Set number of threads for genesis.
    genesis::utils::Options::get().number_of_threads( threads_ );
}

void GlobalOptions::print( CLI::App const& app ) const
{
    if( verbosity() == 0 ) {
        return;
    }

    // Print our nice header.
    std::cout << gappa_header() << "\n";
    if( verbosity() == 1 ) {
        return;
    }

    // TODO sub sub commands are not printed here:
    // More verbose output.
    std::cout << "Invocation:        " << command_line() << "\n";
    for( auto const& sub : app.get_subcommands() ) {
        std::cout << "Subcommand:        " << sub->get_name() << "\n";
    }
    std::cout << "Threads:           " << threads() << "\n";
    std::cout << "\n";
}

std::string GlobalOptions::command_line() const
{
    std::string ret = "";
    for (size_t i = 0; i < command_line_.size(); ++i) {
        ret += ( i==0 ? "" : " " ) + command_line_[i];
    }
    return ret;
}

size_t GlobalOptions::verbosity() const
{
    return ( verbosity_quiet_ == true ) ? 0 : (( verbosity_cnt_ > 0 ) ? verbosity_cnt_ + 1 : verbosity_ );
}

size_t GlobalOptions::threads() const
{
    return threads_;
}

// =================================================================================================
//      Global Instance
// =================================================================================================

GlobalOptions global_options;
