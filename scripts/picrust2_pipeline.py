#!/usr/bin/env python

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.0.0-b.6"

import argparse
from os import path
import sys
import time
from tempfile import TemporaryDirectory
from picrust2.place_seqs import place_seqs_pipeline
from picrust2.wrap_hsp import castor_hsp_workflow
from picrust2.default import (default_fasta, default_tree, default_tables,
                              default_map, default_regroup_map,
                              default_pathway_map)
from picrust2.run_minpath import run_minpath_pipeline
from picrust2.metagenome_pipeline import run_metagenome_pipeline
from picrust2.util import (make_output_dir, check_files_exist,
                           add_descrip_col)

HSP_METHODS = ['mp', 'emp_prob', 'pic', 'scp', 'subtree_average']

parser = argparse.ArgumentParser(

    description="Wrapper for full 16S PICRUSt2 pipeline. Will output "
                "predictions for E.C. numbers, KEGG orthologs, and MetaCyc "
                "pathway abundances and coverages by default. Note that this "
                "is a convenience script, but you can also run each step "
                "individually. Descriptions of gene families and pathways "
                "will be added automatically to output files unless custom "
                "trait tables are input. Please note that currently THIS SCRIPT IS EXPERIMENTAL.",

    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-s', '--study_fasta', metavar='PATH', required=True,
                    type=str, help='FASTA of unaligned study sequences')

parser.add_argument('-i', '--input', metavar='PATH', required=True, type=str,
                    help='Input table of sequence abundances (BIOM or TSV '
                         'format)')

parser.add_argument('-o', '--output', metavar='PATH', required=True,
                    type=str, help='Output folder for final files')

parser.add_argument('--threads', type=int, default=1,
                    help='Number of threads to use (default: %(default)d).')

parser.add_argument('-r', '--ref_msa', metavar='PATH', type=str,
                    default=default_fasta,
                    help='FASTA of aligned reference sequences (default: %(default)s).')

parser.add_argument('-t', '--tree', metavar='PATH', type=str,
                    default=default_tree,
                    help='Input tree based on aligned reference sequences. '
                         '(default: %(default)s).')

parser.add_argument('--in_traits', type=str.upper, default='EC,KO',
                    help='Comma-delimited list (with no spaces) of which gene '
                         'families to predict from this set: COG, EC, KO, '
                         'PFAM, TIGRFAM. Note that E.C. numbers will always '
                         'be predicted unless --no_pathways is set '
                         '(default: %(default)s).')

parser.add_argument('--custom_trait_tables', type=str, metavar='PATH', default=None,
                    help='Optional path to custom trait tables with gene families '
                         'as columns and genomes as rows (overrides '
                         '--in_traits setting) to be used for hidden-state '
                         'prediction. Multiple tables can be specified by '
                         'delimiting filenames by commas. Importantly, the '
                         'first custom table specified will be used for '
                         'inferring pathway abundances. Typically this command '
                         'would be used with a custom marker gene table '
                         '(--marker_gene_table) as well')

parser.add_argument('--marker_gene_table', type=str, metavar='PATH', 
                    default=default_tables["16S"],
                    help='Path to marker gene copy number table (16S copy '
                         'numbers by default).')

parser.add_argument('--pathway_map', metavar='MAP', type=str,
                    default=default_pathway_map,
                    help='MinPath mapfile. The default mapfile maps MetaCyc '
                         'reactions to prokaryotic pathways '
                         '(default: %(default)s).')

parser.add_argument('--no_pathways', default=False, action='store_true',
                    help='Flag to indicate that pathways should NOT be '
                         'inferred (otherwise they will be inferred by '
                         'default). Predicted E.C. number abundances are used '
                         'to infer pathways when default reference files are '
                         'used.')

parser.add_argument('--regroup_map', metavar='ID_MAP',
                    default=default_regroup_map, type=str,
                    help='Mapfile of ids to regroup gene families to before '
                         'running MinPath. The default mapfile is for '
                         'regrouping E. C. numbers to MetaCyc reactions '
                         '(default: %(default)s).')

parser.add_argument('--no_regroup', default=False, action="store_true",
                    help='Do not regroup input gene families to reactions '
                         'as specified in the regrouping mapfile. This option '
                         'should only be used if you are using custom '
                         'reference and/or mapping files.')

parser.add_argument('--stratified', default=False, action='store_true',
                    help='Flag to indicate that stratified tables should be '
                         'generated at all steps (will increase run-time)')

parser.add_argument('--max_nsti', metavar='INT', type=int, default=2,
                    help='Sequences with NSTI values above this value will '
                         'be excluded (default: %(default)d).')

parser.add_argument('--min_reads', metavar='INT', type=int, default=1,
                    help='Minimum number of reads across all samples for '
                         'each input ASV. ASVs below this cut-off will be '
                         'counted as part of the \"RARE\" category in the '
                         'stratified output (default: %(default)d).')

parser.add_argument('--min_samples', metavar='INT', type=int, default=1,
                    help='Minimum number of samples that an ASV needs to be '
                         'identfied within. ASVs below this cut-off will be '
                         'counted as part of the \"RARE\" category in the '
                         'stratified output (default: %(default)d).')

parser.add_argument('-m', '--hsp_method', default='mp',
                    choices=HSP_METHODS,
                    help='HSP method to use.' +
                    '"mp": predict discrete traits using max parsimony. '
                    '"emp_prob": predict discrete traits based on empirical '
                    'state probabilities across tips. "subtree_average": '
                    'predict continuous traits using subtree averaging. '
                    '"pic": predict continuous traits with phylogentic '
                    'independent contrast. "scp": reconstruct continuous '
                    'traits using squared-change parsimony (default: '
                    '%(default)s).')

parser.add_argument('-n', '--calculate_NSTI', default=False,
                    action='store_true',
                    help='Calculate NSTI and add to output ' +
                         'file')

parser.add_argument('-c', '--confidence', default=False, action='store_true',
                    help='Output 95 percent confidence ' +
                         'intervals (only possible for mk_model, emp_prob, '
                         'and mp settings)')

parser.add_argument('--seed', default=100, type=int,
                    help='Seed to make output reproducible, which is '
                         'necessary for the mp and emp_prob methods '
                         '(default: %(default)d).')

parser.add_argument('--no_gap_fill', default=False, action="store_true",
                    help='Do not perform gap filling before predicting '
                         'pathway abundances (Gap filling is on otherwise by '
                         'default.')

parser.add_argument('--per_sequence_contrib', default=False, action="store_true",
                    help='Run MinPath on the gene families contributed by '
                    'each sequence (i.e. a predicted genome) individually. '
                    'This will only matter --per_sequence_contrib is set. '
                    'Note this will GREATLY increase the runtime, but will '
                    'output the predicted pathway abundance contributed by the '
                    'predicted gene families in each predicted genome alone '
                    '(i.e. not the contribution to the community-wide '
                    'abundance). Pathway coverage stratified by contributing '
                    'sequence will also be output when this option is set '
                    '(default: %(default)d).')

parser.add_argument('--print_cmds', default=False, action='store_true',
                    help='If specified, print out wrapped commands to screen')

parser.add_argument('-v', '--version', default=False, action='version',
                    version="%(prog)s " + __version__)


def main():

    args = parser.parse_args()

    # Get start time.
    start_time = time.time()

    # Check that input files exist.
    check_files_exist([args.study_fasta, args.input])

    # Make output folder.
    make_output_dir(args.output)

    out_tree = path.join(args.output, "out.tre")

    if args.custom_trait_tables is None:

        # Check that specified functional categories are allowed.
        FUNC_TRAIT_OPTIONS = ['COG', 'EC', 'KO', 'PFAM', 'TIGRFAM']
        funcs = args.in_traits.split(",")
        for func in funcs:
            if func not in FUNC_TRAIT_OPTIONS:
                sys.exit("Error - specified category " + func + " is not one of "
                         "the default categories.")

        # Add EC to this set if pathways are to be predicted.
        if "EC" not in funcs and not args.no_pathways:
            funcs.append("EC")

        rxn_func = "EC"

        func_tables = default_tables

    else:
        funcs = []
        func_tables = {}

        table_i = 0

        for custom in args.custom_trait_tables.split(","):

            func_id = path.splitext(path.basename(custom))[0]
            funcs.append(func_id)
            func_tables[func_id] = custom

            if table_i == 0:
                rxn_func = func_id
                table_i += 1

    # Append marker as well, since this also needs to be run.
    funcs.append("marker")
    func_tables["marker"] = args.marker_gene_table

    # Methods for discrete trait prediction with CI enabled.
    discrete_set = set(['emp_prob', 'mp'])

    if args.confidence and args.hsp_method in discrete_set:
        ci_setting = True
    else:
        ci_setting = False

    gap_fill_opt = not args.no_gap_fill

    with TemporaryDirectory() as temp_dir:

        print("Placing sequences onto reference tree.")

        place_seqs_pipeline(study_fasta=args.study_fasta,
                            ref_msa=args.ref_msa,
                            tree=args.tree,
                            out_tree=out_tree,
                            threads=args.threads,
                            papara_output=None,
                            out_dir=temp_dir,
                            chunk_size=5000,
                            print_cmds=args.print_cmds)

        print("Finished placing sequences on output tree: " + out_tree)

    # Get predictions for all specified functions and keep track of outfiles.
    predicted_funcs = {}

    for func in funcs:

        # Only output NSTI in 16S table.
        nsti_setting = False
        if func == "marker" and args.calculate_NSTI:
            nsti_setting = True

        print("Running hidden-state prediction for " + func)

        hsp_table, ci_table = castor_hsp_workflow(tree_path=out_tree,
                                                  trait_table_path=func_tables[func],
                                                  hsp_method=args.hsp_method,
                                                  calc_nsti=nsti_setting,
                                                  calc_ci=ci_setting,
                                                  check_input=False,
                                                  num_proc=args.threads,
                                                  ran_seed=args.seed)

        count_outfile = path.join(args.output, func + "_predicted.tsv")

        # Add "_nsti" to filename if output.
        if nsti_setting:
            count_outfile = path.join(args.output, func + "_nsti_predicted.tsv")

        # Keep track of output file name for next step of pipeline.
        predicted_funcs[func] = count_outfile

        print("Writing out predicted gene family abundances to " + count_outfile)

        hsp_table.to_csv(path_or_buf=count_outfile, index_label="sequence", sep="\t")

        # Output the CI file as well if option set.
        if ci_setting:
            ci_outfile = path.join(args.output, func + "_predicted_ci.tsv")
            print("Writing out predicted gene family CIs to " + ci_outfile)
            ci_table.to_csv(path_or_buf=ci_outfile, index_label="sequence",
                            sep="\t")

    marker_infile = predicted_funcs["marker"]

    # Loop over each function again and run metagenome pipeline.
    for func in funcs:

        if func == "marker":
            continue

        func_infile = predicted_funcs[func]

        func_output_dir = path.join(args.output, func + "_metagenome_out")

        print("Running metagenome pipeline for " + func)

        # Infer metagenome abundances per-sample.
        with TemporaryDirectory() as temp_dir:

            # Pass arguments to key function and get predicted functions
            # stratified and unstratified by genomes.
            strat_pred, unstrat_pred = run_metagenome_pipeline(input_biom=args.input,
                                                               function=func_infile,
                                                               marker=marker_infile,
                                                               out_dir=func_output_dir,
                                                               max_nsti=args.max_nsti,
                                                               min_reads=args.min_reads,
                                                               min_samples=args.min_samples,
                                                               strat_out=args.stratified,
                                                               proc=args.threads,
                                                               output_normfile=True)

            print("Writing metagenome output files for " + func + " to: " +
                  func_output_dir)

            # Generate output table filepaths and write out pandas dataframe.
            unstrat_outfile = path.join(func_output_dir, "pred_metagenome_unstrat.tsv")

            unstrat_pred.index.name = "function"
            unstrat_pred.reset_index(inplace=True)

            if args.custom_trait_tables is None:
                unstrat_pred = add_descrip_col(inputfile=unstrat_pred,
                                               mapfile=default_map[func],
                                               in_df=True)

            unstrat_pred.to_csv(path_or_buf=unstrat_outfile, sep="\t", index=False)

            # Write out stratified table only if that option was specified.
            if args.stratified:
                strat_outfile = path.join(func_output_dir, "pred_metagenome_strat.tsv")
                strat_pred.reset_index(inplace=True)

                if args.custom_trait_tables is None:
                    strat_pred = add_descrip_col(inputfile=strat_pred,
                                                 mapfile=default_map[func],
                                                 in_df=True)

                strat_pred.to_csv(path_or_buf=strat_outfile, sep="\t", index=False)

    # Infer pathway abundances and coverages unless --no_pathways set.
    if not args.no_pathways:

        if args.stratified:
            in_metagenome = path.join(args.output,
                                      rxn_func + "_metagenome_out",
                                      "pred_metagenome_strat.tsv")
        else:
            in_metagenome = path.join(args.output,
                                      rxn_func + "_metagenome_out",
                                      "pred_metagenome_unstrat.tsv")

        print("Inferring MetaCyc pathways from predicted functions in this "
              "file: " + in_metagenome)

        with TemporaryDirectory() as temp_dir:
            unstrat_abun, unstrat_cov, strat_abun, strat_cov = run_minpath_pipeline(
                                                            inputfile=in_metagenome,
                                                            mapfile=default_pathway_map,
                                                            regroup_mapfile=default_regroup_map,
                                                            proc=args.threads,
                                                            out_dir=temp_dir,
                                                            gap_fill=gap_fill_opt,
                                                            per_sequence_contrib=args.per_sequence_contrib,
                                                            print_cmds=args.print_cmds)

        pathways_out = path.join(args.output, "pathways_out")

        make_output_dir(pathways_out)

        print("Writing predicted pathway abundances and coverages to " +
              pathways_out)

        # Write output files.
        unstrat_abun_outfile = path.join(pathways_out, "path_abun_unstrat.tsv")
        unstrat_abun.reset_index(inplace=True)

        if args.custom_trait_tables is None:
            unstrat_abun = add_descrip_col(inputfile=unstrat_abun,
                                           mapfile=default_map["METACYC"],
                                           in_df=True)

        unstrat_abun.to_csv(path_or_buf=unstrat_abun_outfile,  sep="\t",
                            index=False)

        unstrat_cov_outfile = path.join(pathways_out, "path_cov_unstrat.tsv")
        unstrat_cov.reset_index(inplace=True)

        if args.custom_trait_tables is None:
            unstrat_cov = add_descrip_col(inputfile=unstrat_cov,
                                        mapfile=default_map["METACYC"],
                                        in_df=True)

        unstrat_cov.to_csv(path_or_buf=unstrat_cov_outfile,  sep="\t",
                           index=False)

        # Write stratified output only if something besides None was returned.
        if strat_abun is not None:
            strat_abun_outfile = path.join(pathways_out, "path_abun_strat.tsv")

            if args.custom_trait_tables is None:
                strat_abun = add_descrip_col(inputfile=strat_abun,
                                             mapfile=default_map["METACYC"],
                                             in_df=True)
            strat_abun.to_csv(path_or_buf=strat_abun_outfile,  sep="\t",
                              index=False)

        if strat_cov is not None:
            strat_cov_outfile = path.join(pathways_out, "path_cov_strat.tsv")

            if args.custom_trait_tables is None:
                strat_cov = add_descrip_col(inputfile=strat_cov,
                                            mapfile=default_map["METACYC"],
                                            in_df=True)

            strat_cov.to_csv(path_or_buf=strat_cov_outfile,  sep="\t",
                             index=False)

    # Print out elapsed time.
    elapsed_time = time.time() - start_time
    print("Completed PICRUSt2 pipeline in " + "%.2f" % elapsed_time + " seconds.")


if __name__ == "__main__":
    main()
