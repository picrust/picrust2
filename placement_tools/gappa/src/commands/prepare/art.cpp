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

#include "commands/prepare/art.hpp"

#include "options/global.hpp"

#include "CLI/CLI.hpp"

#include "genesis/sequence/counts.hpp"
#include "genesis/sequence/formats/fasta_input_iterator.hpp"
#include "genesis/sequence/formats/fasta_writer.hpp"
#include "genesis/sequence/functions/consensus.hpp"
#include "genesis/sequence/functions/entropy.hpp"
#include "genesis/sequence/functions/labels.hpp"
#include "genesis/sequence/sequence_set.hpp"
#include "genesis/sequence/sequence.hpp"

#include "genesis/taxonomy/formats/taxonomy_reader.hpp"
#include "genesis/taxonomy/formats/taxopath_generator.hpp"
#include "genesis/taxonomy/formats/taxopath_parser.hpp"
#include "genesis/taxonomy/functions/entropy_data.hpp"
#include "genesis/taxonomy/functions/entropy.hpp"
#include "genesis/taxonomy/functions/taxonomy.hpp"
#include "genesis/taxonomy/functions/taxopath.hpp"
#include "genesis/taxonomy/iterator/preorder.hpp"
#include "genesis/taxonomy/taxon.hpp"
#include "genesis/taxonomy/taxonomy.hpp"
#include "genesis/taxonomy/taxopath.hpp"

#include "genesis/utils/core/std.hpp"
#include "genesis/utils/io/output_stream.hpp"
#include "genesis/utils/text/string.hpp"

#ifdef GENESIS_OPENMP
#   include <omp.h>
#endif

#include <cctype>
#include <map>
#include <fstream>

// =================================================================================================
//      Typedefs
// =================================================================================================

// =================================================================================================
//      Setup
// =================================================================================================

void setup_art( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto opt = std::make_shared<ArtOptions>();
    auto sub = app.add_subcommand(
        "art",
        "Generate consensus sequences from a sequence database according to the ART method."
    );

    // -----------------------------------------------------------
    //     Input Data
    // -----------------------------------------------------------

    // Taxonomy file
    auto tax_file_opt = sub->add_option(
        "--taxonomy-file",
        opt->taxonomy_file,
        "File that lists the taxa of the database."
    );
    tax_file_opt->required();
    tax_file_opt->check( CLI::ExistingFile );
    tax_file_opt->group( "Input" );

    // Sequence file
    auto seq_file_opt = sub->add_option(
        "--sequence-file",
        opt->sequence_file,
        "Fasta file containing the sequences of the database."
    );
    seq_file_opt->required();
    seq_file_opt->check( CLI::ExistingFile );
    seq_file_opt->group( "Input" );

    // opt->sequence_input.add_fasta_input_opt_to_app( sub );

    // -----------------------------------------------------------
    //     Entropy pruning option
    // -----------------------------------------------------------

    // Target size
    auto target_size_opt = sub->add_option(
        "--target-size",
        opt->target_taxonomy_size,
        "Target size of how many taxa to select for building consensus sequences."
    );
    target_size_opt->required();
    target_size_opt->group( "Taxonomy Expansion" );

    // Sub Taxopath
    sub->add_option(
        "--sub-taxonomy",
        opt->sub_taxopath,
        "If a taxopath from the taxonomy is provided, only the respective sub-taxonomy is used."
    )->group( "Taxonomy Expansion" );

    // Min subclade size
    // auto min_subclade_size_opt =
    sub->add_option(
        "--min-subclade-size",
        opt->min_subclade_size,
        "Minimal size of sub-clades. Everything below is expanded."
    )->group( "Taxonomy Expansion" );

    // Max subclade size
    // auto max_subclade_size_opt =
    sub->add_option(
        "--max-subclade-size",
        opt->max_subclade_size,
        "Maximal size of a non-expanded sub-clades. Everything bigger is first expanded."
    )->group( "Taxonomy Expansion" );

    // Min tax level
    // auto min_tax_level_opt =
    sub->add_option(
        "--min-tax-level",
        opt->min_tax_level,
        "Minimal taxonomic level. Taxa below this level are always expanded."
    )->group( "Taxonomy Expansion" );

    // Allow approximation
    // auto allow_approx_opt =
    sub->add_flag(
        "--allow-approximation",
        opt->allow_approximation,
        "Allow to expand taxa that help getting closer to the --target-size, even if they are not "
        "the ones with the highest entropy."
    )->group( "Taxonomy Expansion" );

    // Allow approximation
    // auto allow_approx_opt =
    sub->add_flag(
        "--no-taxa-selection",
        opt->no_taxa_selection,
        "If set, no taxa selection using entropy is performed. Instead, all taxa on all levels/ranks are "
        "used and consensus sequences for all of them are calculated. This is useful for testing and "
        "to try out new ideas."
    )->group( "Taxonomy Expansion" );

    // -----------------------------------------------------------
    //     Consensus options
    // -----------------------------------------------------------

    // Consensus Method
    auto cons_meth_opt = sub->add_set_ignore_case(
        "--consensus-method",
        opt->consensus_method,
        { "majorities", "cavener", "threshold" },
        "Consensus method to use for combining sequences.",
        true
    )->group( "Consensus Method" );

    // Consensus Treshold
    auto cons_thresh_opt = sub->add_option(
        "--consensus-threshold",
        opt->consensus_threshold,
        "Threshold value to use with --consensus-method threshold. Has to be in [ 0.0, 1.0 ].",
        true
    );
    cons_thresh_opt->needs( cons_meth_opt );
    cons_thresh_opt->check( CLI::Range( 0.0, 1.0 ));
    cons_thresh_opt->group( "Consensus Method" );

    // -----------------------------------------------------------
    //     Output Options
    // -----------------------------------------------------------

    opt->output.add_output_dir_opt_to_app( sub );

    // Write info files
    sub->add_flag(
        "--write-info-files",
        opt->write_info_files,
        "If set, two additional info files are written, containing the new pruned taxonomy, "
        "as well as the entropy of all clades of the original taxonomy."
    )->group( "Output" );

    // -----------------------------------------------------------
    //     Callback
    // -----------------------------------------------------------

    // Set the run function as callback to be called when this subcommand is issued.
    // Hand over the options by copy, so that their shared ptr stays alive in the lambda.
    sub->set_callback( [ opt ]() {
        run_art( *opt );
    });
}

// =================================================================================================
//      Read Taxonomy
// =================================================================================================

genesis::taxonomy::Taxonomy read_taxonomy( ArtOptions const& options )
{
    using namespace genesis::sequence;
    using namespace genesis::taxonomy;

    // User output.
    if( global_options.verbosity() >= 1 ) {
        std::cout << "Reading taxonomy and preparing entropy calculations.\n";
    }

    // Get alignment length.
    auto it = FastaInputIterator();
    it.from_file( options.sequence_file );
    auto const seq_len = it->size();

    // Read from file.
    auto tax = Taxonomy();
    TaxonomyReader().from_file( options.taxonomy_file, tax );
    sort_by_name( tax );

    // Init a pointer to the whole taxonomy.
    Taxonomy* subtax = &tax;

    // If the user only wants a sub taxon, overwrite the pointer.
    if( ! options.sub_taxopath.empty() ) {

        // Find the sub taxon and assign it to the taxonomy as our new main taxon.
        auto const taxopath = TaxopathParser().from_string( options.sub_taxopath );
        auto subtaxon = find_taxon_by_taxopath( tax, taxopath );
        subtax = subtaxon;

        if( subtax == nullptr ) {
            throw std::runtime_error(
                "Taxon " + options.sub_taxopath + " not found in the taxonomy."
            );
        }

        // We need to set data for the selected main sub clade taxon,
        // as this is not done in the preorder iterator that comes next.
        subtaxon->reset_data( EntropyTaxonData::create() );
        subtaxon->data<EntropyTaxonData>().counts = SiteCounts( "ACGT", seq_len );
    }

    // Create a Sequence Count objeect for each taxon below the pointer.
    // This might allocate quite a lot of memory!
    auto add_sequence_counts_to_taxonomy = [&]( Taxon& taxon ){
        taxon.reset_data( EntropyTaxonData::create() );
        taxon.data<EntropyTaxonData>().counts = SiteCounts( "ACGT", seq_len );
    };
    assert( subtax );
    preorder_for_each( *subtax, add_sequence_counts_to_taxonomy );

    // User output.
    if( global_options.verbosity() >= 1 ) {
        std::cout << "Taxonomy contains a total of " << total_taxa_count( tax ) << " taxa, ";
        std::cout << "with " << taxa_count_lowest_levels( tax ) << " taxa at the lowest level.\n";

        if( ! options.sub_taxopath.empty() ) {
            std::cout << "The selected Subtaxonomy contains a total of " ;
            std::cout << total_taxa_count( *subtax ) << " taxa, ";
            std::cout << "with " << taxa_count_lowest_levels( *subtax );
            std::cout << " taxa at the lowest level.\n";
        }
    }

    return tax;
}

// =================================================================================================
//      Fill Site Counts
// =================================================================================================

void fill_site_counts( ArtOptions const& options, genesis::taxonomy::Taxonomy& tax )
{
    using namespace genesis::sequence;
    using namespace genesis::taxonomy;

    // User output.
    if( global_options.verbosity() >= 1 ) {
        std::cout << "Reading sequences.\n";
    }

    // User output prep. We count how often each char occurs in the sequences,
    // how many sequences weree processed in total, how many of those were not found at all in the
    // taxonomy, and how many were not part of the specified subtaxonomy (if specified at all).
    std::map<char, size_t> char_counts;
    size_t total_seqs_count = 0;
    size_t no_tax_seqs_count = 0;
    size_t not_subtax_seqs_count = 0;

    // Prepare helpers.
    auto taxopath_parser = TaxopathParser();
    auto fasta_reader = FastaReader();
    fasta_reader.site_casing( FastaReader::SiteCasing::kToUpper );

    // Iterate sequences.
    auto it = FastaInputIterator( fasta_reader );
    for( it.from_file( options.sequence_file ); it; ++it ) {

        // User output.
        // If we have verbose output, count characters.
        if( global_options.verbosity() >= 2 ) {
            for( auto const& s : *it ) {
                ++char_counts[s];
            }

            // Also, output a progress. Could be done nicer in the future.
            if( total_seqs_count % 100000 == 0 ) {
                std::cout << "At sequence " << total_seqs_count << "\n";
            }
        }
        ++total_seqs_count;

        // Get taxo path of the sequence.
        // We offer two versions: the sequence label is just the path, or it starts at the first space.
        std::string taxopath_str;
        auto const delim = it->label().find_first_of( " \t" );
        if( delim == std::string::npos ) {
            taxopath_str = it->label();
        } else {
            taxopath_str = it->label().substr( delim + 1 );
        }

        // Parse the taxo path and find it in the taxonomy.
        // If the first attempt fails, remove the last element (assumed to be species level),
        // and try again. If we fail again, we cannot use this sequence.
        auto taxopath = taxopath_parser( taxopath_str );
        auto taxp = find_taxon_by_taxopath( tax, taxopath );
        if( taxp == nullptr ) {
            taxopath.pop_back();
            taxp = find_taxon_by_taxopath( tax, taxopath );
        }
        if( taxp == nullptr ) {
            ++no_tax_seqs_count;
            if( global_options.verbosity() >= 3 ) {
                std::cout << "Sequence " << it->label() << " not found in the taxonomy!\n";
            }
            continue;
        }

        // Now that we have found the taxon of that sequence, check whether it is part of the
        // specified subtaxonomy. If no subtaxonomy was specified, all are valid.
        // We do this by testing whether the taxon has data, because read_taxonomy() only sets
        // data entries for the subtaxonomy.
        if( ! taxp->has_data() ) {
            ++not_subtax_seqs_count;
            if( global_options.verbosity() >= 3 ) {
                std::cout << "Sequence " << it->label() << " not part of the subtaxonomy.\n";
            }
            continue;
        }

        // Accummulate counts for all taxonomic ranks.
        // We go up in the taxonomy and add counts to all super-clades as well,
        // until we reach the super taxon that is not part of the selected sub-clade
        // (if a sub clade was specified. otherweise, it just goes all the way up).
        auto cur_tax = taxp;
        do {
            cur_tax->data<EntropyTaxonData>().counts.add_sequence( *it );
            cur_tax = cur_tax->parent();
        } while( cur_tax != nullptr && cur_tax->has_data() );
    }

    // User output.
    if( global_options.verbosity() >= 1 || no_tax_seqs_count > 0 ) {
        std::cout << "Processed " << total_seqs_count << " sequences.\n";
        if( no_tax_seqs_count > 0 ) {
            std::cout << "Thereof, " << no_tax_seqs_count << " sequences were not found in the taxonomy.\n";
        }
        if( not_subtax_seqs_count > 0 ) {
            if( no_tax_seqs_count == 0 ) {
                std::cout << "Thereof, ";
            } else {
                std::cout << "Furthermore, ";
            }

            std::cout << not_subtax_seqs_count << " sequences were skipped because they are ";
            std::cout << "not part of the specified subtaxonomy.\n";
        }
    }
    if( global_options.verbosity() >= 2 ) {
        std::cout << "Character counts in the sequences:\n";
        size_t sum = 0;
        for( auto const& count : char_counts ) {
            std::cout << "    " << count.first << ": " << count.second << "\n";
            sum += count.second;
        }

        auto amb = sum - char_counts['-'];
        amb -= char_counts['A'] + char_counts['C'];
        amb -= char_counts['G'] + char_counts['T'] + char_counts['U'];
        auto const amb_per = 100.0 * static_cast<double>( amb ) / static_cast<double>( sum );
        if( amb > 0 ) {
            std::cout << "There were " << amb << " (" << amb_per << "%) ambiguous sites, ";
            std::cout << "which were counted as gaps.\n";
        }
    }
    if( char_counts['U'] > char_counts['T'] ) {
        // The above condition might add U or T to the char counts map.
        // We don't use it any more afterwards, so that's okay.
        std::cout << "Warning: There are more 'U' sites in the sequences than 'T' sites. ";
        std::cout << "Are you sure that the sites are properly converted to 'T'?\n";
    }
}

// =================================================================================================
//      Calculate Entropy
// =================================================================================================

void calculate_entropy( ArtOptions const& options, genesis::taxonomy::Taxonomy& tax )
{
    using namespace genesis::sequence;
    using namespace genesis::taxonomy;

    if( options.no_taxa_selection ) {
        std::cout << "Skipping entropy calculation due to --no-taxa-selection being set.\n";
        return;
    }

    // User output.
    if( global_options.verbosity() >= 1 ) {
        std::cout << "Calculating entropy.\n";
    }

    // Set the default options.
    // auto opt = SiteEntropyOptions::kWeighted;
    auto opt = SiteEntropyOptions::kIncludeGaps;

    // Calculate! Skip those that do not have data, that is, which are not part of the subtaxonomy.
    // This is a simple way of testing for the subtaxonomy, instead of finding it again here.
    auto calc_entropies = [&]( Taxon& t ) {
        if( ! t.has_data() ) {
            return;
        }

        auto const& counts = t.data<EntropyTaxonData>().counts;
        t.data<EntropyTaxonData>().entropy = averaged_entropy( counts, false, opt );
    };
    preorder_for_each( tax, calc_entropies );
}

// =================================================================================================
//      Select Taxa
// =================================================================================================

void select_taxa( ArtOptions const& options, genesis::taxonomy::Taxonomy& tax )
{
    using namespace genesis::sequence;
    using namespace genesis::taxonomy;

    // User output.
    if( global_options.verbosity() >= 1 ) {
        std::cout << "Selecting taxa based on entropy.\n";
    }

    // Init a pointer to the whole taxonomy.
    Taxonomy* subtax = &tax;

    // If the user only wants a sub taxon, overwrite the pointer.
    if( ! options.sub_taxopath.empty() ) {
        auto const taxopath = TaxopathParser().from_string( options.sub_taxopath );
        subtax = find_taxon_by_taxopath( tax, taxopath );

        if( subtax == nullptr ) {
            throw std::runtime_error(
                "Taxon " + options.sub_taxopath + " not found in the taxonomy."
            );
        }
    }

    if( options.no_taxa_selection ) {

        // If we do not run taxa selection by entropy, make all taxa border, that is, select all.
        auto set_status_to_border = [&]( Taxon& t ) {
            if( ! t.has_data() ) {
                return;
            }

            t.data<EntropyTaxonData>().status = EntropyTaxonData::PruneStatus::kBorder;
        };
        preorder_for_each( *subtax, set_status_to_border );

    } else {

        // Get algo settings.
        PruneByEntropySettings prune_settings;
        prune_settings.min_subtaxonomy_size = options.min_subclade_size;
        prune_settings.max_subtaxonomy_size = options.max_subclade_size;
        prune_settings.min_border_level     = options.min_tax_level;
        prune_settings.allow_approximation  = options.allow_approximation;

        // Run Forrest, run!
        prune_by_entropy( *subtax, options.target_taxonomy_size, prune_settings );
        if( ! validate_pruned_taxonomy( *subtax ) ) {
            throw std::runtime_error( "Something went wrong, the selected taxa are inconsistent.\n" );
        }
    }

    // User output.
    if( global_options.verbosity() >= 1 ) {
        auto const border_cnt = count_taxa_with_prune_status(
            *subtax, EntropyTaxonData::PruneStatus::kBorder
        );
        std::cout << "Selected " << border_cnt << " taxa for which to build consensus sequences.\n";
    }
}

// =================================================================================================
//      Generate Consensus Sequences
// =================================================================================================

void generate_consensus_sequences( ArtOptions const& options, genesis::taxonomy::Taxonomy const& tax )
{
    using namespace genesis::sequence;
    using namespace genesis::taxonomy;

    // User output.
    if( global_options.verbosity() >= 1 ) {
        std::cout << "Generating consensus sequences.\n";
    }

    // Helper function.
    auto write_fasta_sequence = [] ( std::ofstream& out, std::string const& name, std::string const& sites ) {
        out << ">" << name << "\n";
        for (size_t i = 0; i < sites.length(); i += 80) {
            out << sites.substr(i, 80) << "\n";
        }
    };

    // Prepare output.
    // TODO check with file overwrite settings
    std::ofstream cons_seq_os;
    auto const fn = options.output.out_dir() + options.consensus_sequence_file;
    genesis::utils::file_output_stream( fn, cons_seq_os );
    auto const tax_gen = TaxopathGenerator();

    // Calculate consensus sequences and write them.
    auto make_consensus_sequences = [&]( Taxon const& t ) {

        // Skip taxa that are not in the subtaxonomy.
        if( ! t.has_data() ) {
            return;
        }

        // Only interested in the border taxa.
        if( t.data<EntropyTaxonData>().status != EntropyTaxonData::PruneStatus::kBorder ) {
            return;
        }

        // Prep.
        auto const name = sanitize_label( tax_gen(t) );
        auto const& counts = t.data<EntropyTaxonData>().counts;

        // Consensus sequence.
        std::string sites;
        if( options.consensus_method == "majorities" ) {
            sites = consensus_sequence_with_majorities( counts );
        } else if( options.consensus_method == "cavener" ) {
            sites = consensus_sequence_cavener( counts );
        } else if( options.consensus_method == "threshold" ) {
            sites = consensus_sequence_with_threshold( counts, options.consensus_threshold );
        } else {
            throw CLI::ConversionError( "Unknown consensus method: " + options.consensus_method );
        }
        write_fasta_sequence( cons_seq_os, name, sites );

    };
    preorder_for_each( tax, make_consensus_sequences );
    cons_seq_os.close();
}

// =================================================================================================
//      Write Taxonomy Info
// =================================================================================================

void write_info_files( ArtOptions const& options, genesis::taxonomy::Taxonomy const& tax )
{
    using namespace genesis::taxonomy;

    if( ! options.write_info_files ) {
        return;
    }

    // TODO check with file overwrite settings

    // Prepare entropy output.
    std::ofstream entropy_os;
    auto const entropy_fn = options.output.out_dir() + options.entropy_info_file;
    genesis::utils::file_output_stream( entropy_fn, entropy_os );
    entropy_os << "Taxon\tStatus\tChild_Taxa\tLowest_Level_Taxa\tTotal_Taxa\tSequences\tEntropy\n";

    // Prepare taxonomy output.
    std::ofstream taxonomy_os;
    auto const taxonomy_fn = options.output.out_dir() + options.taxonomy_info_file;
    genesis::utils::file_output_stream( taxonomy_fn, taxonomy_os );
    taxonomy_os << "Taxon\tChild_Taxa\tLowest_Level_Taxa\tTotal_Taxa\n";

    // Write to files.
    auto const gen = TaxopathGenerator();
    auto print_taxon_info = [&]( Taxon const& t )
    {
        // Skip taxa that are not in the subtaxonomy.
        if( ! t.has_data() ) {
            return;
        }

        // Calculate values.
        auto const name = gen( t );
        auto const total_chldrn = total_taxa_count( t );
        auto const lowest_chldrn = taxa_count_lowest_levels( t );
        auto const added_seqs = t.data<EntropyTaxonData>().counts.added_sequences_count();
        auto const entr = t.data<EntropyTaxonData>().entropy;

        // Status: was the taxon selected or not.
        std::string status = "-";
        if( t.data<EntropyTaxonData>().status == EntropyTaxonData::PruneStatus::kOutside ) {
            status = "Outside";
        }
        if( t.data<EntropyTaxonData>().status == EntropyTaxonData::PruneStatus::kBorder ) {
            status = "Selected";
        }
        if( t.data<EntropyTaxonData>().status == EntropyTaxonData::PruneStatus::kInside ) {
            status = "Inside";
        }

        // For all taxa, write out entropy info.
        entropy_os << name;
        entropy_os << "\t" << status;
        entropy_os << "\t" << t.size();
        entropy_os << "\t" << lowest_chldrn;
        entropy_os << "\t" << total_chldrn;
        entropy_os << "\t" << added_seqs;
        entropy_os << "\t" << entr;
        entropy_os << "\n";

        // Write all inner and border taxa to the taxnomy file.
        if( t.data<EntropyTaxonData>().status != EntropyTaxonData::PruneStatus::kOutside ) {
            taxonomy_os << name;
            taxonomy_os << "\t" << t.size();
            taxonomy_os << "\t" << lowest_chldrn;
            taxonomy_os << "\t" << total_chldrn;
            taxonomy_os << "\n";
        }
    };
    preorder_for_each( tax, print_taxon_info );

    // Wrap up.
    entropy_os.close();
    taxonomy_os.close();
}

// =================================================================================================
//      Run
// =================================================================================================

void run_art( ArtOptions const& options )
{
    using namespace genesis::utils;

    // User output.
    if( global_options.verbosity() >= 2 ) {
        std::cout << "Taxonomy file: " << options.taxonomy_file << "\n";
        std::cout << "Sequence file: " << options.sequence_file << "\n";

        std::cout << "Target taxonomy size: " << options.target_taxonomy_size << "\n";
        std::cout << "Consensus method: " << options.consensus_method;
        if( options.consensus_method == "threshold" ) {
            std::cout << " (" << options.consensus_threshold << ")";
        }
        std::cout << "\n";

        if( options.min_subclade_size > 0 ) {
            std::cout << "Min subclade size: " << options.min_subclade_size << "\n";
        }
        if( options.max_subclade_size > 0 ) {
            std::cout << "Max subclade size: " << options.max_subclade_size << "\n";
        }
        if( options.min_tax_level > 0 ) {
            std::cout << "Min taxonomic level: " << options.min_tax_level << "\n";
        }
        if( options.allow_approximation ) {
            std::cout << "Allowing approximation.\n";
        }
        if( ! options.sub_taxopath.empty() ) {
            std::cout << "Subtaxonomy: " << options.sub_taxopath << "\n";
        }

        std::cout << "\n";
    }

    // Check output files
    options.output.check_nonexistent_output_files({ options.consensus_sequence_file });
    if( options.write_info_files ) {
        options.output.check_nonexistent_output_files({
            options.entropy_info_file,
            options.taxonomy_info_file
        });
    }

    // Run the whole thing!
    auto taxonomy = read_taxonomy( options );
    fill_site_counts( options, taxonomy );
    calculate_entropy( options, taxonomy );
    select_taxa( options, taxonomy );
    generate_consensus_sequences( options, taxonomy );
    write_info_files( options, taxonomy );

    // User output.
    if( global_options.verbosity() >= 1 ) {
        std::cout << "Finished.\n";
    }
}
