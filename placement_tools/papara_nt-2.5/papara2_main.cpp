/*
 * Copyright (C) 2009-2012 Simon A. Berger
 * 
 * This file is part of papara.
 * 
 *  papara is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  papara is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with papara.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "blast_partassign.h"

#include "ivymike/getopt.h"
#include "ivymike/time.h"

#include "papara.h"

using namespace papara;

// namespace papara {
// log_stream lout;
// }

void print_banner( std::ostream &os ) {
    os << "______    ______    ______        ___\n";
    os << "| ___ \\   | ___ \\   | ___ \\   |\\ | |   \n";
    os << "| |_/ /_ _| |_/ /_ _| |_/ /__ | \\| | \n";
    os << "|  __/ _` |  __/ _` |    // _` |\n";
    os << "| | | (_| | | | (_| | |\\ \\ (_| |\n";
    os << "\\_|  \\__,_\\_|  \\__,_\\_| \\_\\__,_|\n";
	os << "  papara_core version " << papara::get_version_string() << "\n";

}

void pad( std::ostream &os, size_t w ) {
    for( size_t i = 0; i < w; ++i ) {
        os << " ";
    }
}

void print_help( std::ostream &os, const std::vector<std::string> &options, const std::vector<std::string> &text ) {
    assert( options.size() == text.size() );

    size_t max_opt_size = 0;

    for( size_t i = 0; i < options.size(); ++i ) {
        max_opt_size = std::max( max_opt_size, options[i].size() );
    }


    size_t pad_right = 3;
    size_t col_width = 2 + max_opt_size + pad_right;
    for( size_t i = 0; i < options.size(); ++i ) {
        pad( os, 2 );
        os << options[i];
        pad( os, col_width - 2 - options[i].size() );
        //os << text[i];
        const std::string &t = text[i];


        for( size_t j = 0; j < t.size(); ++j ) {
            if( t[j] == '@') {
                os << "\n";
                pad( os, col_width );
            } else {
                os << t[j];
            }
        }

        os << "\n";
    }
    
    
}


void print_help( std::ostream &os ) {
    std::vector<std::string> options;
    std::vector<std::string> text;

    options.push_back( "-t <tree file>" );
    text.push_back( "File name of the reference tree (newick format)");

    options.push_back( "-s <ref alignment>" );
    text.push_back( "File name of the reference alignment (phylip format)");

    options.push_back( "-q <query seqs.>" );
    text.push_back( "File name of the query sequences (fasta format)");

    options.push_back( "-a" );
    text.push_back( "Sequences are protein data");

    options.push_back( "-j <num threads>" );
    text.push_back( "Specify number of threads (default:1)");

    options.push_back( "-n <run name>" );
    text.push_back( "Specify filename suffix of the output@files (default: \"default\")");

    options.push_back( "-r" );
    text.push_back( "Turn of writing RA-side gaps in the output file.");

    options.push_back( "-p" );
    text.push_back( "User defined scoring scheme: <open>:<extend>:<match>:<match cg>@The default scores correspond to '-p -3:-1:2:-3'" );
    
    print_help( os, options, text );

}




template<typename pvec_t, typename seq_tag>
void run_papara( const std::string &qs_name, const std::string &alignment_name, const std::string &tree_name, size_t num_threads, const std::string &run_name, bool ref_gaps, const papara_score_parameters &sp, bool write_fasta, partassign::part_assignment *part_assign, const std::pair<size_t,size_t> &fixed_qs_bounds ) {

    ivy_mike::perf_timer t1;

    queries<seq_tag> qs(qs_name.c_str());

    
    
    
    t1.add_int();
    references<pvec_t,seq_tag> refs( tree_name.c_str(), alignment_name.c_str(), &qs );

    
    
    t1.add_int();

    qs.preprocess();
   
    t1.add_int();

    refs.remove_full_gaps();
    refs.build_ref_vecs();

    if( part_assign != 0 ) {
        if( ref_gaps ) {
            std::cout << "REMARK: using per-gene alignment deactivates reference-side gaps!\n";
            ref_gaps = false;
        }
        
        
        //qs.init_partition_assignments( *part_assign );
        std::vector<std::pair<size_t,size_t> > qs_bounds = partassign::resolve_qs_bounds( refs, qs, *part_assign );
        
        
        qs.set_per_qs_bounds( qs_bounds );
        
        
    } else if( fixed_qs_bounds.first != size_t(-1) ) {
	ref_gaps = false;
	std::cout << "fixed bounds " << fixed_qs_bounds.first << " " << fixed_qs_bounds.second << "\n";
        std::vector<std::pair<size_t,size_t> > qs_bounds( qs.size(), fixed_qs_bounds );
        qs.set_per_qs_bounds( qs_bounds );
    }
    
    
    
    
    t1.add_int();

//     t1.print();

    const size_t num_candidates = 0;

    scoring_results res( qs.size(), scoring_results::candidates(num_candidates) );


    lout << "scoring scheme: " << sp.gap_open << " " << sp.gap_extend << " " << sp.match << " " << sp.match_cgap << "\n";

    driver<pvec_t,seq_tag>::calc_scores(num_threads, refs, qs, &res, sp );

    std::string score_file(filename(run_name, "alignment"));
    std::string quality_file(filename(run_name, "quality"));
    std::string cands_file(filename(run_name, "cands"));


    size_t pad = 1 + std::max(qs.max_name_length(), refs.max_name_length());

//     std::ofstream os( score_file.c_str() );
//     assert( os.good() );

    std::ofstream os_qual( quality_file.c_str() );
    assert( os_qual.good() );


    std::ofstream os_cands;
    if( num_candidates > 1 ) {
        os_cands.open( cands_file.c_str() );
        assert( os_cands.good());
    }



    std::unique_ptr<papara::output_alignment> oa;
    if( write_fasta ) {
        oa.reset( new papara::output_alignment_fasta( score_file.c_str() ));
    } else {
        oa.reset( new papara::output_alignment_phylip( score_file.c_str() ));
    }
    
    //refs.write_seqs(os, pad);
    //     driver<pvec_t,seq_tag>::align_best_scores( os, os_qual, os_cands, qs, refs, res, pad, ref_gaps, sp );
    driver<pvec_t,seq_tag>::align_best_scores_oa( oa.get(), qs, refs, res, pad, ref_gaps, sp );
    
}


void print_commandline( std::ostream &os, char **argv, int argc ) {

	std::stringstream ss;
	std::copy( argv, argv + argc, std::ostream_iterator<char *>( os, " " ) );
	os << std::endl;
}


int main( int argc, char *argv[] ) {

    
    namespace igo = ivy_mike::getopt;

    ivy_mike::getopt::parser igp;

    std::string opt_tree_name;
    std::string opt_alignment_name;
    std::string opt_qs_name;
    std::string opt_user_parameters;
    
    std::string opt_blast_hits;
    std::string opt_partitions;
    std::string opt_partition_name;
    
    bool opt_use_cgap;
    int opt_num_threads;
    std::string opt_run_name;
    bool opt_write_testbench;
    bool opt_force_overwrite;
    bool opt_aa;
    bool opt_no_ref_gaps;
    bool opt_print_help;
    bool opt_write_fasta;
    
    igp.add_opt( 't', igo::value<std::string>(opt_tree_name) );
    igp.add_opt( 's', igo::value<std::string>(opt_alignment_name) );
    igp.add_opt( 'q', igo::value<std::string>(opt_qs_name) );
    igp.add_opt( 'c', igo::value<bool>(opt_use_cgap, true).set_default(false) );
    igp.add_opt( 'a', igo::value<bool>(opt_aa, true).set_default(false) );
    igp.add_opt( 'j', igo::value<int>(opt_num_threads).set_default(1) );
    igp.add_opt( 'n', igo::value<std::string>(opt_run_name).set_default("default") );
    igp.add_opt( 'b', igo::value<bool>(opt_write_testbench, true).set_default(false) );
    igp.add_opt( 'f', igo::value<bool>(opt_force_overwrite, true).set_default(false) );
    igp.add_opt( 'r', igo::value<bool>(opt_no_ref_gaps, true).set_default(false) );
    igp.add_opt( 'p', igo::value<std::string>(opt_user_parameters).set_default("") );
    igp.add_opt( 'h', igo::value<bool>(opt_print_help, true).set_default(false) );
    igp.add_opt( 'g', igo::value<bool>(opt_write_fasta, true).set_default(false) );
    igp.add_opt( 'l', igo::value<std::string>(opt_blast_hits) );
    igp.add_opt( 'x', igo::value<std::string>(opt_partitions) );
    igp.add_opt( 'k', igo::value<std::string>(opt_partition_name) );
    
    igp.parse(argc,argv);

    if( opt_print_help ) {
        print_banner( std::cerr );
        std::cerr << "\nUser options:\n";
        print_help( std::cerr);
        return 0;
        
    }
         
    
    if( igp.opt_count('t') != 1 || igp.opt_count('s') != 1 ) {
        print_banner(std::cerr);

        std::cerr << "missing options -t and/or -s (-q is optional)\n";

        print_help( std::cerr );
        return 0;
    }
    
    // optional accelration by blast hits/partition file
    std::unique_ptr<partassign::part_assignment> part_assignment;
    
    std::pair<size_t,size_t> fixed_qs_bounds(-1,-1);
    
    if( igp.opt_count('l') == 1  ) {
        if( igp.opt_count('l') != igp.opt_count('x') ) {
            std::cerr << "options -l and -x have to be used together (or not at all)\n";
            print_help( std::cerr );
            return 0;
        }
        
        std::ifstream blast_is( opt_blast_hits.c_str() );
        if( !blast_is.good() ) {
            std::cerr << "can not open blast hits file\n";
            return 0;
        }
        
        std::ifstream part_is( opt_partitions.c_str() );
        if( !part_is.good() ) {
            std::cerr << "can not open partition file\n";
            return 0;
        }
        
        
        part_assignment.reset( new partassign::part_assignment( blast_is, part_is ));
    } else if(igp.opt_count('k') == 1 ) {
        std::ifstream part_is( opt_partitions.c_str() );
        if( !part_is.good() ) {
            std::cerr << "can not open partition file\n";
            return 0;
        }
        
        fixed_qs_bounds = partassign::partition_bounds( part_is, opt_partition_name );
        std::cout << fixed_qs_bounds.first << " " << fixed_qs_bounds.second << "\n";
        if( fixed_qs_bounds.first == size_t(-1) || fixed_qs_bounds.second == size_t(-1) ) {
            std::cerr << "partition not found: " + opt_partition_name + "\n";
            return 0;
        }
    } else if( igp.opt_count('x') == 1 ) {
        std::cerr << "options -x needs either -l or -k\n";
        print_help( std::cerr );
        return 0;
    }
        
    
    
    ivy_mike::timer t;

//    const char *qs_name = 0;
//    if( !opt_qs_name.empty() ) {
//        qs_name = opt_qs_name.c_str();
//    }
//
    std::string log_filename = filename( opt_run_name, "log" );

    if( opt_run_name != "default" && !opt_force_overwrite && file_exists(log_filename.c_str()) ) {
        std::cout << "log file already exists for run '" << opt_run_name << "'\n";
        return 0;
    }
    
    
    papara::add_log_tee papara_log_cout( std::cout );
    
    std::ofstream logs( log_filename.c_str());
    if( !logs ) {
        std::cout << "could not open logfile for writing: " << log_filename << std::endl;
        return 0;
    }
    
    
    papara::add_log_tee papara_log_file( logs );

    papara::lout << "papara called as:\n";
    print_commandline( papara::lout, argv, argc ); 

    const bool ref_gaps = !opt_no_ref_gaps;


    papara_score_parameters sp = papara_score_parameters::default_scores();
    if( !opt_user_parameters.empty() ) {
        std::cout << "meeeeep\n";
        
        sp = papara_score_parameters::parse_scores(opt_user_parameters.c_str());
    }
    
    
    if( opt_use_cgap ) {

        if( opt_aa ) {
            run_papara<pvec_cgap, tag_aa>( opt_qs_name, opt_alignment_name, opt_tree_name, opt_num_threads, opt_run_name, ref_gaps, sp, opt_write_fasta, part_assignment.get(), fixed_qs_bounds );
        } else {
            run_papara<pvec_cgap, tag_dna>( opt_qs_name, opt_alignment_name, opt_tree_name, opt_num_threads, opt_run_name, ref_gaps, sp, opt_write_fasta, part_assignment.get(), fixed_qs_bounds );
        }
    } else {
        if( opt_aa ) {
            run_papara<pvec_pgap, tag_aa>( opt_qs_name, opt_alignment_name, opt_tree_name, opt_num_threads, opt_run_name, ref_gaps, sp, opt_write_fasta, part_assignment.get(), fixed_qs_bounds );
        } else {
            run_papara<pvec_pgap, tag_dna>( opt_qs_name, opt_alignment_name, opt_tree_name, opt_num_threads, opt_run_name, ref_gaps, sp, opt_write_fasta, part_assignment.get(), fixed_qs_bounds );
        }
    }

    std::cout << t.elapsed() << std::endl;
    lout << "SUCCESS " << t.elapsed() << std::endl;

    
    
    
    return 0;
//     getchar();
}

