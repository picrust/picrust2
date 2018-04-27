#!/usr/bin/env python

from __future__ import division

__license__ = "GPL"
__version__ = "2-alpha.7"

import argparse
from picrust2.run_minpath import run_minpath_pipeline
import pandas as pd

parser = argparse.ArgumentParser(

    description="Runs MinPath on a table of E.C. numbers",

    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-i', '--input', metavar='PATH', required=True, type=str,
                    help='Input BIOM table of E.C. number abundances')

parser.add_argument('-m', '--map', metavar='PATH', required=True, type=str,
                    help='MinPath mapfile')

parser.add_argument('-o', '--output', metavar='PATH', required=True, type=str,
                    help='Output file containing pathway abundances')

parser.add_argument('--keep_tmp', default=False, action="store_true",
                    help='If specified, keep temporary folder')

parser.add_argument('-t', '--threads', default=1, type=int,
                    help='Number of threads')

parser.add_argument('--tmp_dir', metavar='PATH', type=str,
                    help='Temporary directory for running MinPath')

parser.add_argument('--print_cmds', default=False, action="store_true",
                    help='If specified, print out wrapped commands to screen')


def main():

    args = parser.parse_args()

    metacyc_predictions = run_minpath_pipeline(inputfile=args.input,
                                               mapfile=args.map,
                                               keep_tmp=args.keep_tmp,
                                               threads=args.threads,
                                               tmp_dir=args.tmp_dir,
                                               print_cmds=args.print_cmds)

    metacyc_predictions.to_csv(path_or_buf=args.output,
                               sep="\t",
                               index_label="pathway")


if __name__ == "__main__":
    main()
