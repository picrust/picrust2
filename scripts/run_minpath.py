#!/usr/bin/env python

from __future__ import division

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2-alpha.10"

import argparse
from picrust2.run_minpath import run_minpath_pipeline
from tempfile import TemporaryDirectory
from picrust2.util import make_output_dir, check_files_exist

parser = argparse.ArgumentParser(

    description="Runs MinPath on a table of E.C. numbers",

    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-i', '--input', metavar='IN_TABLE', required=True,
                    type=str,
                    help='Input TSV table of E.C. number abundances')

parser.add_argument('-m', '--map', metavar='MAP', required=True, type=str,
                    help='MinPath mapfile')

parser.add_argument('-o', '--out_prefix', metavar='PREFIX', required=True,
                    type=str, help='Prefix for stratified and unstratified ' +
                                   'pathway abundances.')

parser.add_argument('--intermediate', metavar='DIR', type=str, default=None,
                    help='Output folder for intermediate files (wont be ' +
                         'kept unless this option is set.')

parser.add_argument('-p', '--proc', default=1, type=int,
                    help='Number of threads')

parser.add_argument('--print_cmds', default=False, action="store_true",
                    help='If specified, print out wrapped commands to screen')


def main():

    args = parser.parse_args()

    # Check that input files exist.
    check_files_exist([args.input, args.map])

    # If intermediate output directory set then create and output there.
    # Otherwise make a temporary directory for the intermediate files.
    if args.intermediate:
        make_output_dir(args.intermediate)

        unstrat_out, strat_out = run_minpath_pipeline(inputfile=args.input,
                                                      mapfile=args.map,
                                                      proc=args.proc,
                                                      out_dir=args.intermediate,
                                                      print_cmds=args.print_cmds)
    else:
        with TemporaryDirectory() as temp_dir:
                unstrat_out, strat_out = run_minpath_pipeline(inputfile=args.input,
                                                              mapfile=args.map,
                                                              proc=args.proc,
                                                              out_dir=temp_dir,
                                                              print_cmds=args.print_cmds)

    # Write output files.
    unstrat_outfile = args.out_prefix + "_unstrat_path.tsv"
    strat_outfile = args.out_prefix + "_strat_path.tsv"

    unstrat_out.to_csv(path_or_buf=unstrat_outfile,  sep="\t",
                       index_label="pathway")

    strat_out.to_csv(path_or_buf=strat_outfile,  sep="\t", index=False)


if __name__ == "__main__":
    main()
