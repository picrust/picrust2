#!/usr/bin/env python

__license__ = "GPL"
__version__ = "2-alpha.6"

import argparse
from picrust2.wrap_hsp import castor_hsp_wrapper
from picrust2.util import make_output_dir_for_file

HSP_METHODS = ['mp', 'emp_prob', 'pic', 'scp', 'subtree_average']

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
hsp.py -i precalculated/prokaryotic/ec_counts_mean_round_var.txt -t study.tre

''', formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-i', '--observed_trait_table', metavar='PATH',
                    required=True, type=str,
                    help='The input trait table describing directly ' +
                         'observed traits (e.g. sequenced genomes) in ' +
                         'tab-delimited format')

parser.add_argument('-t', '--tree', metavar='PATH', required=True, type=str,
                    help='The full reference tree, in newick format')


parser.add_argument('-o', '--output_trait_table', metavar='PATH', type=str,
                    default='predicted_traits.tsv',
                    help='The output filepath for trait predictions')

parser.add_argument('-r', '--rds_outfile', metavar='PATH', type=str,
                    default='trait_state_probs.rds',
                    help='The output filepath for the R object ' +
                         'containing predicted trait state ' +
                         'probabilities')

parser.add_argument('--ci_out', metavar='PATH', type=str,
                    default='predicted_traits_ci.tsv',
                    help='The output filepath for confidence intervals ' +
                         'trait predictions (if -c option set)')

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

    # Methods for discrete trait prediction with CI enabled.
    discrete_set = set(['emp_prob', 'mp'])

    if args.confidence and args.hsp_method in discrete_set:
        ci_setting = True
    else:
        ci_setting = False

    hsp_table, ci_table = castor_hsp_wrapper(tree_path=args.tree,
                                             trait_table_path=args.observed_trait_table,
                                             hsp_method=args.hsp_method,
                                             calc_nsti=args.calculate_NSTI,
                                             calc_ci=ci_setting,
                                             check_input=args.check,
                                             num_cores=args.processes,
                                             rds_outfile=args.rds_outfile,
                                             ran_seed=args.seed,
                                             HALT_EXEC=args.debug)

    # Output the table to file.
    make_output_dir_for_file(args.output_trait_table)
    hsp_table.to_csv(path_or_buf=args.output_trait_table, index_label="tips",
                     sep="\t")

    # Output the CI file as well if option set.
    if ci_setting:
        make_output_dir_for_file(args.ci_out)
        ci_table.to_csv(path_or_buf=args.ci_out, index_label="tips", sep="\t")


if __name__ == "__main__":
    main()
