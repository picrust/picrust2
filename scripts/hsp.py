#!/usr/bin/env python

import argparse
from picrust2.wrap_hsp import castor_hsp_workflow
from picrust2.util import make_output_dir_for_file, check_files_exist
from picrust2.default import default_tables

HSP_METHODS = ['mp', 'emp_prob', 'pic', 'scp', 'subtree_average']

TRAIT_OPTIONS = ['16S', 'COG', 'EC', 'KO', 'PFAM', 'TIGRFAM', 'PHENO']

parser = argparse.ArgumentParser(

    description="This script performs hidden state prediction on tips in "
                "the input tree with unknown trait values. Typically this "
                "script is used to predict the copy number of gene families "
                "present in the predicted genome for each amplicon sequence "
                "variant, given a tree and a set of known trait values. "
                "This script outputs a table of trait predictions.",
    epilog='''
Usage example:
hsp.py -n -t out.tre -i 16S -o 16S_predicted_and_nsti.tsv.gz --processes 1

''', formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-t', '--tree', metavar='PATH', required=True, type=str,
                    help='The full reference tree in newick format containing '
                         'both study sequences (i.e. ASVs or OTUs) and '
                         'reference sequences.')

parser.add_argument('-o', '--output', metavar='PATH', type=str, required=True,
                    help='Output table with predicted abundances per study '
                         'sequence in input tree. If the extension \".gz\" '
                         'is added the table will automatically be gzipped.')

parser.add_argument('-i', '--in_trait', type=str.upper, choices=TRAIT_OPTIONS,
                    help='Specifies which default trait table should be '
                          'used. Use the --observed_trait_table option '
                          'to input a non-default trait table.')

parser.add_argument('--observed_trait_table', metavar='PATH', type=str,
                    help='The input trait table describing directly '
                         'observed traits (e.g. sequenced genomes) in '
                         'tab-delimited format. Necessary if you want to '
                         'use a custom table.')

parser.add_argument('-e', '--edge_exponent', default=0.5, type=float,
                    help='Setting for maximum parisomony hidden-state '
                          'prediction. Specifies weighting transition costs '
                          'by the inverse length of edge lengths. If 0, then '
                          'edge lengths do not influence predictions. Must be '
                          'a non-negative real-valued number (default: '
                          '%(default)f).')

parser.add_argument('--chunk_size', default=500, type=int,
                    help='Number of functions to run at a time on one '
                         'processor. Note that you should consider how many '
                         'processes you have specified before changing this '
                         'option. E.g. if you specify the chunk_size to be '
                         'the total number of functions, 1 processor will '
                         'be used even if you specified more so the job will '
                         'be substantially slower (default: %(default)d).')

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
                    help='Calculate NSTI and add to output file.')

parser.add_argument('--check', default=False, action='store_true',
                    help='Check input trait table before HSP.')

parser.add_argument('-p', '--processes', default=1, type=int,
                    help='Number of processes to run in parallel (default: '
                    '%(default)d).')

parser.add_argument('--seed', default=100, type=int,
                    help='Seed to make output reproducible, which is '
                         'necessary for the emp_prob method '
                         '(default: %(default)d).')

parser.add_argument('--verbose', default=False, action='store_true',
                    help='If specified, print out wrapped commands and other '
                         'details to screen.')

parser.add_argument('-v', '--version', default=False, action='version',
                    version="%(prog)s " + __version__)


def main():

    args = parser.parse_args()

    # Determine which input trait table was specified. If neither a default
    # or custom table was specified then throw an error.
    if args.in_trait and args.observed_trait_table:
        raise RuntimeError(
            "Only one of the arguments --in_trait and --observed_trait_table "
            "can be specified, but currently both are set.")
    elif args.in_trait:
        trait_table = default_tables[args.in_trait]
    elif args.observed_trait_table:
        trait_table = args.observed_trait_table
    else:
        raise RuntimeError(
            "A default input trait table needs to be specified with the "
            "--in_trait option, or alternatively a custom table can be "
            "specified with the --observed_trait_table option")

    # Check that input filenames exist.
    check_files_exist([args.tree, trait_table])

    # No longer support outputting CIs with this script.
    ci_setting = False

    hsp_table, ci_table = castor_hsp_workflow(tree_path=args.tree,
                                              trait_table_path=trait_table,
                                              hsp_method=args.hsp_method,
                                              edge_exponent=args.edge_exponent,
                                              chunk_size=args.chunk_size,
                                              calc_nsti=args.calculate_NSTI,
                                              calc_ci=ci_setting,
                                              check_input=args.check,
                                              num_proc=args.processes,
                                              ran_seed=args.seed,
                                              verbose=args.verbose)

    # Output the table to file.
    make_output_dir_for_file(args.output)
    hsp_table.to_csv(path_or_buf=args.output, index_label="sequence",
                     sep="\t", compression="infer")


if __name__ == "__main__":
    main()
