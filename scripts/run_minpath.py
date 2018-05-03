#!/usr/bin/env python

from __future__ import division

__license__ = "GPL"
__version__ = "2-alpha.8"

import argparse
from picrust2.run_minpath import run_minpath_pipeline
from tempfile import TemporaryDirectory
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

parser.add_argument('--out_dir', metavar='PATH', type=str, default=None,
                    help='Output folder for intermediate files (wont be ' +
                         'kept unless this option is set.')

parser.add_argument('-t', '--threads', default=1, type=int,
                    help='Number of threads')

parser.add_argument('--print_cmds', default=False, action="store_true",
                    help='If specified, print out wrapped commands to screen')


def main():

    args = parser.parse_args()

    # If intermediate output directory set then create and output there.
    # Otherwise make a temporary directory for the intermediate files.

    if args.out_dir:
        make_output_dir(args.out_dir)
        
        metacyc_predictions = run_minpath_pipeline(inputfile=args.input,
                                                   mapfile=args.map,
                                                   threads=args.threads,
                                                   out_dir=args.out_dir,
                                                   print_cmds=args.print_cmds)
    else:
        with TemporaryDirectory() as temp_dir:
                metacyc_predictions = run_minpath_pipeline(inputfile=args.input,
                                                           mapfile=args.map,
                                                           threads=args.threads,
                                                           out_dir=temp_dir,
                                                           print_cmds=args.print_cmds)

    metacyc_predictions.to_csv(path_or_buf=args.output,
                               sep="\t",
                               index_label="pathway")

if __name__ == "__main__":
    main()
