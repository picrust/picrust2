#!/usr/bin/env python

__license__ = "GPL"
__version__ = "2-alpha.8"

import argparse
from picrust2.wrap_hsp import castor_hsp_wrapper
from picrust2.util import make_output_dir_for_file, get_picrust_project_dir

HSP_METHODS = ['mp', 'emp_prob', 'pic', 'scp', 'subtree_average']

TRAIT_OPTIONS = ['16S', 'COG', 'EC', 'KO', 'PFAM', 'TIGRFAM']

# Inititalize default trait table files.
project_dir = get_picrust_project_dir()

default_tables = {"16S" : project_dir + \
                  "/precalculated/prokaryotic/16S_counts_mean_round_var.txt",
                  "COG" : project_dir + \
                  "/precalculated/prokaryotic/cog_counts_mean_round_var.txt",
                  "EC" : project_dir + \
                  "/precalculated/prokaryotic/ec_counts_mean_round_var.txt",
                  "KO" : project_dir + \
                  "/precalculated/prokaryotic/ko_counts_mean_round_var.txt",
                  "PFAM" : project_dir + \
                  "/precalculated/prokaryotic/pfam_counts_mean_round_var.txt",
                  "TIGRFAM": project_dir + \
                  "/precalculated/prokaryotic/tfam_counts_mean_round_var.txt"}

parser = argparse.ArgumentParser(

    description="This script performs hidden state " +
                "prediction on tips in the input tree " +
                "with unknown trait values. Typically " +
                "this script is used to predict the " +
                "abundance of gene families present in " +
                "each taxon, given a tree and a set of " +
                "known trait values. This script " +
                "outputs a table of trait predictions. " +
                "Note that this script assumes that the input " +
                "trait values will include \"0\" counts.",

    epilog='''

Usage example:
hsp.py -t study.tre -i 16S -o 16S_predicted_traits

''', formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-t', '--tree', metavar='PATH', required=True, type=str,
                    help='The full reference tree, in newick format')

parser.add_argument('-o', '--output_prefix', metavar='PATH', type=str,
                    required=True,
                    help='Prefix for output filenames (RDS, predicted count ' +
                         'table and optionally a table of CIs)')

parser.add_argument('-i', '--in_trait', type=str.upper, choices=TRAIT_OPTIONS,
                    help='Specifies which default trait table should be ' +
                          'used. Use the --observed_trait_table option ' +
                          'to input a non-default trait table.')

parser.add_argument('--observed_trait_table', metavar='PATH', type=str,
                    help='The input trait table describing directly ' +
                         'observed traits (e.g. sequenced genomes) in ' +
                         'tab-delimited format. Necessary if you want to ' +
                         'use a custom table.')

parser.add_argument('-m', '--hsp_method', default='mp',
                    choices=HSP_METHODS,
                    help='HSP method to use.' +
                    '"mp": predict discrete traits using max parsimony. ' +
                    '"emp_prob": predict discrete traits based on empirical ' +
                    'state probabilities across tips. "subtree_average": ' +
                    'predict continuous traits using subtree averaging. ' +
                    '"pic": predict continuous traits with phylogentic ' +
                    'independent contrast. "scp": reconstruct continuous ' +
                    'traits using squared-change parsimony')

parser.add_argument('-n', '--calculate_NSTI', default=False,
                    action='store_true',
                    help='Calculate NSTI and add to output ' +
                         'file')

parser.add_argument('-c', '--confidence', default=False, action='store_true',
                    help='Output 95 percent confidence ' +
                         'intervals (only possible for mk_model, emp_prob, ' +
                         'and mp settings)')

parser.add_argument('--check', default=False, action='store_true',
                    help='Check input trait table before HSP')

parser.add_argument('-p', '--processes', default=1, type=int,
                    help='Number of processes to run in parallel')

parser.add_argument('--seed', default=None, type=int,
                    help='Seed to make output reproducible ' +
                         '(necessary for mp and emp_prob methods)')

parser.add_argument('--debug', default=False, action='store_true',
                    help='Run in debugging mode')

parser.add_argument('-v', '--version', default=False, action='version',
                    version="%(prog)s " + __version__)


def main():

    args = parser.parse_args()

    # Determine which input trait table was specified. If neither a default 
    # or custom table was specified then throw an error.
    if args.in_trait:
        trait_table = default_tables[args.in_trait]
    elif args.observed_trait_table:
        trait_table = args.observed_trait_table
    else:
        raise RuntimeError(
            "A default input trait table needs to be specified with the " +
            "--in_trait option, or alternatively a custom table can be " +
            "specified with the --observed_trait_table option")

    # Methods for discrete trait prediction with CI enabled.
    discrete_set = set(['emp_prob', 'mp'])

    if args.confidence and args.hsp_method in discrete_set:
        ci_setting = True
    else:
        ci_setting = False

    rds_outfile = args.output_prefix + ".rds"
    count_outfile = args.output_prefix + ".tsv"
    ci_outfile = args.output_prefix + "_ci.tsv"

    hsp_table, ci_table = castor_hsp_wrapper(tree_path=args.tree,
                                             trait_table_path=trait_table,
                                             hsp_method=args.hsp_method,
                                             calc_nsti=args.calculate_NSTI,
                                             calc_ci=ci_setting,
                                             check_input=args.check,
                                             num_cores=args.processes,
                                             rds_outfile=rds_outfile,
                                             ran_seed=args.seed,
                                             HALT_EXEC=args.debug)

    # Output the table to file.
    make_output_dir_for_file(count_outfile)
    hsp_table.to_csv(path_or_buf=count_outfile, index_label="tips",
                     sep="\t")

    # Output the CI file as well if option set.
    if ci_setting:
        make_output_dir_for_file(ci_outfile)
        ci_table.to_csv(path_or_buf=ci_outfile, index_label="tips", sep="\t")


if __name__ == "__main__":
    main()
