#!/usr/bin/env python

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.0.0-b.7"

import argparse
from os import path
import sys
import time
from picrust2.default import (default_fasta, default_tree, default_tables,
                              default_map, default_regroup_map,
                              default_pathway_map)
from picrust2.util import make_output_dir
from picrust2.pipeline import full_pipeline

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

parser.add_argument('--no_descrip', default=False, action='store_true',
                    help='Do not add function descriptions to output tables.')

parser.add_argument('--verbose', default=False, action='store_true',
                    help='If specified, print out wrapped commands to screen')

parser.add_argument('-v', '--version', default=False, action='version',
                    version="%(prog)s " + __version__)


def main():

    # Get start time.
    start_time = time.time()

    args = parser.parse_args()
    
    func_dfs, unstrat_abun, unstrat_cov, strat_abun, strat_cov = full_pipeline(study_fasta=args.study_fasta,
                                                                              input_table=args.input,
                                                                              output_folder=args.output,
                                                                              threads=args.threads,
                                                                              ref_msa=args.ref_msa,
                                                                              tree=args.tree,
                                                                              in_traits=args.in_traits,
                                                                              custom_trait_tables=args.custom_trait_tables,
                                                                              marker_gene_table=args.marker_gene_table,
                                                                              pathway_map=args.pathway_map,
                                                                              no_pathways=args.no_pathways,
                                                                              regroup_map=args.regroup_map,
                                                                              no_regroup=args.no_regroup,
                                                                              stratified=args.stratified,
                                                                              max_nsti=args.max_nsti,
                                                                              min_reads=args.min_reads,
                                                                              min_samples=args.min_samples,
                                                                              hsp_method=args.hsp_method,
                                                                              calculate_NSTI=args.calculate_NSTI,
                                                                              confidence=args.confidence,
                                                                              seed=args.seed,
                                                                              no_gap_fill=args.no_gap_fill,
                                                                              per_sequence_contrib=args.per_sequence_contrib,
                                                                              no_descrip=args.no_descrip,
                                                                              verbose=args.verbose)

    for func in func_dfs.keys():

        func_output_dir = path.join(args.output, func + "_metagenome_out")

        print("Writing metagenome output files for " + func + " to: " + func_output_dir)

        unstrat_outfile = path.join(func_output_dir, "pred_metagenome_unstrat.tsv")
        func_dfs[func][0].to_csv(path_or_buf=unstrat_outfile, sep="\t", index=False)    
    
        if func_dfs[func][1] is not None:
            strat_outfile = path.join(func_output_dir, "pred_metagenome_strat.tsv")
            func_dfs[func][1].to_csv(path_or_buf=strat_outfile, sep="\t", index=False)

    if not args.no_pathways:
        pathways_out = path.join(args.output, "pathways_out")

        print("Writing predicted pathway abundances and coverages to " + pathways_out)

        make_output_dir(pathways_out)

        unstrat_abun_outfile = path.join(pathways_out, "path_abun_unstrat.tsv")
        unstrat_cov_outfile = path.join(pathways_out, "path_cov_unstrat.tsv")
        strat_abun_outfile = path.join(pathways_out, "path_abun_strat.tsv")
        strat_cov_outfile = path.join(pathways_out, "path_cov_strat.tsv")

        unstrat_abun.to_csv(path_or_buf=unstrat_abun_outfile,  sep="\t",
                            index=False)

        unstrat_cov.to_csv(path_or_buf=unstrat_cov_outfile,  sep="\t",
                           index=False)

        if strat_abun is not None:
            strat_abun.to_csv(path_or_buf=strat_abun_outfile,  sep="\t",
                              index=False)

        if strat_cov is not None:
            strat_cov.to_csv(path_or_buf=strat_cov_outfile,  sep="\t",
                             index=False)

    # Print out elapsed time.
    elapsed_time = time.time() - start_time
    print("Completed PICRUSt2 pipeline in " + "%.2f" % elapsed_time + " seconds.")

if __name__ == "__main__":
    main()
