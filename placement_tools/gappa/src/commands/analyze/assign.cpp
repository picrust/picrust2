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

#include "commands/analyze/assign.hpp"

#include "options/global.hpp"

#include "CLI/CLI.hpp"

#include "genesis/taxonomy/formats/taxopath_parser.hpp"
#include "genesis/taxonomy/formats/taxopath_generator.hpp"
#include "genesis/taxonomy/functions/taxopath.hpp"
#include "genesis/taxonomy/iterator/preorder.hpp"
#include "genesis/taxonomy/taxopath.hpp"
#include "genesis/taxonomy/taxonomy.hpp"

#include "genesis/utils/formats/csv/reader.hpp"
#include "genesis/utils/io/input_stream.hpp"
#include "genesis/utils/io/output_stream.hpp"

#include "genesis/tree/default/functions.hpp"
#include "genesis/tree/default/tree.hpp"
#include "genesis/tree/default/newick_writer.hpp"
#include "genesis/tree/iterator/postorder.hpp"
#include "genesis/tree/function/functions.hpp"
#include "genesis/tree/tree/edge.hpp"

#include "genesis/placement/sample.hpp"
#include "genesis/placement/pquery.hpp"
#include "genesis/placement/pquery/placement.hpp"


#include <fstream>

#ifdef GENESIS_OPENMP
#   include <omp.h>
#endif

constexpr char UNDETERMINED[] = "N/A";

// =================================================================================================
//      Setup
// =================================================================================================

void setup_assign( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto opt = std::make_shared<AssignOptions>();
    auto sub = app.add_subcommand(
        "assign",
        "Taxonomically assign placed query sequences and output tabulated summarization."
    );

    // Input
    opt->jplace_input.add_jplace_input_opt_to_app( sub )->required();

    sub->add_option(
        "--taxon-file",
        opt->taxon_file,
        "File containing a tab-separated list of taxon to taxonomic string assignments."
    )->check(CLI::ExistingFile)->required()->group("Input");

    // Add specific options.
    sub->add_option(
        "--sub-taxopath",
        opt->sub_taxopath,
        "Taxopath (example: Eukaryota;Animalia;Chordata) by which the high level summary should be filtered. "
        "Doesn't affect intermediate results, and an unfiltered verison will be printed as well."
    )->group("Settings");

    sub->add_option(
        "--distribution-ratio",
        opt->dist_ratio,
        "Ratio by which LWR is split between annotations if an edge has two possible annotations. "
        "Specifies the amount going to the proximal annotation. If not set program will determine "
        "the ratio automatically from the 'distal length' specified per placement."
    )->check(CLI::Range(0.0,1.0))->group("Settings");

    // Output
    opt->output_dir.add_output_dir_opt_to_app( sub );

    // Set the run function as callback to be called when this subcommand is issued.
    // Hand over the options by copy, so that their shared ptr stays alive in the lambda.
    sub->set_callback( [ opt ]() {
        run_assign( *opt );
    });
}

// =================================================================================================
//      Run
// =================================================================================================

using namespace genesis;
using namespace genesis::placement;
using namespace genesis::tree;
using namespace genesis::taxonomy;
using namespace genesis::utils;

static Taxopath intersect( Taxopath const& lhs, Taxopath const& rhs )
{
    Taxopath result;

    for (size_t i = 0; ( i < std::min( lhs.size(), rhs.size() ) ) and ( lhs[i] == rhs[i] ); ++i) {
        result.push_back( lhs[i] );
    }

    if ( result.empty() ) {
        result.push_back( UNDETERMINED );
    }

    return result;
}

// go through the tree in postorder fashion and label inner nodes according to the most common taxonomic rank of the children
static void postorder_label( PlacementTree const& tree, std::vector<Taxopath>& node_labels )
{
    for ( auto it : postorder(tree) ) {
        if ( it.node().is_inner() ) {
            auto const child_1_idx = it.node().link().next().outer().node().index();
            auto const child_2_idx = it.node().link().next().next().outer().node().index();

            assert( not node_labels[ child_1_idx ].empty() );

            node_labels[ it.node().index() ] = intersect( node_labels[ child_1_idx ], node_labels[ child_2_idx ] );
        }
    }
}

static void print_labelled( PlacementTree const& tree,
                            std::vector<Taxopath> const& node_labels,
                            std::string const& file_name )
{
    DefaultTreeNewickWriter writer;
    writer.node_to_element_plugins.push_back(
        [&]( TreeNode const& node, NewickBrokerElement& element ){
            element.comments.emplace_back(
                TaxopathGenerator().to_string( node_labels[node.index()] )
            );
        }
    );
    writer.to_file( tree, file_name );
}

static std::vector<Taxopath> assign_leaf_taxopaths( PlacementTree const& tree,
                                                    std::string const& taxon_file)
{
    TaxopathParser tpp;
    CsvReader csv_reader;
    csv_reader.separator_chars( "\t" );
    std::vector<Taxopath> node_labels(tree.node_count());

    utils::InputStream it( utils::make_unique< utils::FileInputSource >( taxon_file ));
    while (it) {
        auto fields = csv_reader.parse_line( it );

        if ( fields.size() != 2 ) {
            throw std::runtime_error{"A line in the taxon file didn't have two tab separated columns."};
        }

        auto name = fields[0];
        std::string tax_string = fields[1];

        auto node_ptr = find_node( tree, name );

        if ( node_ptr == nullptr ) {
            throw std::runtime_error{"Could not find node with name: " + name};
        }

        node_labels[ node_ptr->index() ] = tpp.from_string( tax_string );
    }

    // check if any leafs weren't assigned a Taxopath
    for ( auto const& node_it : tree.nodes() ) {
        if ( node_it->is_leaf() and node_labels[ node_it->index() ].empty() ) {
            auto name = node_it->data< DefaultNodeData >().name;
            throw std::runtime_error{"The leaf in the tree labelled '" + name
                + "' wasn't assigned a taxonomic path. Did you forget to include it in the taxon file?"};
        }
    }

    return node_labels;
}

static void add_lwr_to_taxonomy(const double lwr,
                                Taxopath const& path,
                                Taxonomy& taxonomy )
{
    auto cur_tax = &add_from_taxopath( taxonomy, path );

    bool first = true;

    do {
        if( not cur_tax->has_data() ) {
            cur_tax->reset_data( AssignTaxonData::create() );
        }

        if ( first ) {
            // add normal elw to this taxon
            cur_tax->data<AssignTaxonData>().LWR += lwr;
            first = false;
        }

        // add accumulated elw up the taxpath
        cur_tax->data<AssignTaxonData>().aLWR += lwr;
        cur_tax = cur_tax->parent();
    } while( cur_tax != nullptr );
}

void print_taxonomy_with_lwr( Taxonomy const& t, std::ostream& stream )
{
    // get total LWR as sum of all top level aLWR
    double sum = 0.0;
    for( auto const& ct : t ) {
        sum += ct.data<AssignTaxonData>().aLWR;
    }

    preorder_for_each( t, [&]( Taxon const& tax ){

        stream  << std::setprecision(4)
                << tax.data<AssignTaxonData>().LWR
                << "\t"
                << tax.data<AssignTaxonData>().LWR / sum
                << "\t"
                << tax.data<AssignTaxonData>().aLWR
                << "\t"
                << tax.data<AssignTaxonData>().aLWR / sum
                << "\t"
                << TaxopathGenerator().to_string( tax )
                << "\n";
    });
}

void print_taxonomy_table( Taxonomy const& t, std::ostream& stream )
{
    stream << "LWR\tfract\taLWR\tafract\ttaxopath\n";
    print_taxonomy_with_lwr( t, stream );
}

Taxonomy& get_subtaxonomy( Taxonomy tax, AssignOptions const& options )
{
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

    assert( subtax != nullptr );

    return *subtax;
}

static void assign( Sample const& sample,
                    std::vector<Taxopath> const& node_labels,
                    AssignOptions const& options,
                    std::string per_pquery_result_file = "" )
{
    bool    const   auto_ratio = ( options.dist_ratio < 0.0 );
    double  const   dist_ratio = options.dist_ratio;
    assert( auto_ratio or ( dist_ratio > 0.0 and dist_ratio <= 1.0 ) );

    auto const& tree = sample.tree();
    Taxonomy diversity;

    std::ofstream per_pquery_out_stream;
    bool const intermediate_results = (not per_pquery_result_file.empty());

    if ( intermediate_results ) {
        genesis::utils::file_output_stream( per_pquery_result_file, per_pquery_out_stream );
    }

    for ( auto const& pq : sample.pqueries() ) {
        Taxonomy per_pq_assignments;
        for ( auto const& p : pq.placements() ) {
            auto const lwr = p.like_weight_ratio;
            // get its adjacent nodes
            auto const& edge = tree.edge_at( p.edge().index() );
            auto const& proximal_node   = edge.primary_node();
            auto const& distal_node     = edge.secondary_node();

            // get the taxopaths
            auto const& proximal_tax    = node_labels[ proximal_node.index() ];
            auto const& distal_tax      = node_labels[ distal_node.index() ];

            double ratio = dist_ratio;
            // determine the ratio
            if ( auto_ratio ) {
                auto const position         = p.proximal_length;
                auto const branch_length    = edge.data<DefaultEdgeData>().branch_length;
                // in percent, how far toward the distal are we?
                auto const toward_distal    = (1 / branch_length) * position;
                // the ratio is effectively "how much lwr mass should go toward the PROXIMAL", so we need to flip it
                ratio = 1.0 - toward_distal;
            }

            // calculate lwr portions
            auto proximal_portion   = lwr * ratio;
            auto distal_portion     = lwr * (1.0 - ratio);


            // add LW to taxopaths of the nodes according to strategy
            // first to the local one
            if ( intermediate_results ) {
                add_lwr_to_taxonomy( proximal_portion,  proximal_tax,   per_pq_assignments );
                add_lwr_to_taxonomy( distal_portion,    distal_tax,     per_pq_assignments );
            }

            // then to the global one
            add_lwr_to_taxonomy( proximal_portion, proximal_tax, diversity );
            add_lwr_to_taxonomy( distal_portion, distal_tax, diversity );
        }

        if ( intermediate_results ) {
            for ( auto const& name : pq.names() ) {
                per_pquery_out_stream << name.name;
            }
            per_pquery_out_stream << std::endl;
            print_taxonomy_with_lwr( per_pq_assignments, per_pquery_out_stream );
        }
    }
    auto out_dir = options.output_dir.out_dir();

    // return diversity profile
    std::ofstream profile;
    genesis::utils::file_output_stream( out_dir + "profile.csv", profile );
    print_taxonomy_table( diversity, profile );

    // constrain to subtaxonomy if specified
    const auto subtaxonomy = get_subtaxonomy( diversity, options );
    // and print to file
    std::ofstream profile_filtered;
    genesis::utils::file_output_stream( out_dir + "profile_filtered.csv", profile_filtered );
    print_taxonomy_table( subtaxonomy, profile_filtered );
}

void run_assign( AssignOptions const& options )
{
    auto out_dir = options.output_dir.out_dir();

    options.jplace_input.print();
    auto sample = options.jplace_input.merged_samples();
    auto tree = sample.tree();

    if ( not is_bifurcating(tree) ) {
        throw std::runtime_error{"Supplied tree is not bifurcating."};
    }

    if( global_options.verbosity() >= 2 ) {
        std::cout << "Getting taxonomic information from: " << options.taxon_file << "\n";
    }

    // vector to hold the per node taxopaths
    // fill the per node taxon assignments
    auto node_labels = assign_leaf_taxopaths(tree, options.taxon_file);

    // assign taxpaths to inner nodes
    postorder_label( tree, node_labels );

    // print taxonomically labelled tree as intermediate result
    print_labelled( tree, node_labels, out_dir + "labelled_tree" );

    // per rank LWR score eval
    assign( sample, node_labels, options, out_dir + "per_pquery_assign" );

    if( global_options.verbosity() >= 1 ) {
        std::cout << "Finished.\n";
    }
}
