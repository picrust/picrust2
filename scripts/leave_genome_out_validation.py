#!/usr/bin/env python

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2-alpha.10"

import argparse
from picrust2.wrap_hsp import castor_hsp_loocv_wrapper

parser = argparse.ArgumentParser(

    description="Given a set of tip names in tree will run HSP when leaving " +
                "that tip out. Expected and predicted values will be output " +
                "along with the Spearman rho between the expected and " +
                "observed and the NSTI value for that tip. Currently " +
                "only HSP with maximum parsimony is enabled",

    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-i', '--observed_trait_table', metavar='PATH',
                    required=True, type=str,
                    help='The input trait table describing directly ' +
                         'observed traits (e.g. sequenced genomes) in ' +
                         'tab-delimited format')

parser.add_argument('-t', '--tree', metavar='PATH', required=True, type=str,
                    help='The full reference tree, in newick format')

parser.add_argument('-n', '--names', metavar="PATH", required=True, type=str,
                    help='File with tip names to leave out - one per line')

parser.add_argument('-o', '--output_trait_table', metavar='PATH', type=str,
                    default='predicted_traits.tsv',
                    help='The output filepath for trait predictions ' +
                         '(default: %(default)s).')

parser.add_argument('--exp_out', metavar="PATH", default='expected_traits.tsv',
                    help='The output filepath for expected trait values ' +
                         '(default: %(default)s).')

parser.add_argument('--pred_out', metavar="PATH",
                    default='expected_traits.tsv',
                    help='The output filepath for predicted trait values ' +
                         '(default: %(default)s).')

parser.add_argument('-m', '--metrics_out', metavar="PATH",
                    default='rho_nsti_out.tsv',
                    help='The output filepath for Rho and NSTI values per ' +
                         'left out genome (default: %(default)s).')

parser.add_argument('-p', '--processes', default=1, type=int,
                    help='Number of processes to run in parallel ' +
                         '(default: %(default)d).')

parser.add_argument('--debug', default=False, action='store_true',
                    help='Run in debugging mode')

parser.add_argument('-v', '--version', default=False, action='version',
                    version="%(prog)s " + __version__)


def main():

    args = parser.parse_args()

    castor_hsp_loocv_wrapper(tree_path=args.tree,
                             trait_table_path=args.observed_trait_table,
                             tips_path=args.names,
                             hsp_method="mp",
                             expected_out_path=args.exp_out,
                             predicted_out_path=args.pred_out,
                             metrics_out_path=args.metrics_out,
                             num_cores=args.processes,
                             HALT_EXEC=args.debug)


if __name__ == "__main__":
    main()
