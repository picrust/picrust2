#!/usr/bin/env python

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.0.0-b.4"

import argparse
from picrust2.run_minpath import run_minpath_pipeline
from tempfile import TemporaryDirectory
from picrust2.util import make_output_dir, check_files_exist
from os import path

parser = argparse.ArgumentParser(

    description="Wrapper for MinPath to infer which pathways are present "
                "given gene family abundances in a sample. Gene families can "
                "be regrouped to a different gene family or reaction type if "
                "needed. Typical usage would be to input E.C. numbers and "
                "regroup them to MetaCyc reations to infer MetaCyc pathways. "
                "Currently only unstructured pathways are supported. Pathways "
                "that are called as present are given an abundance based on "
                "the mean abundance of the top 50% most abundant reactions in "
                "the pathway. Stratified output will also be output if a "
                "stratified metagenome is input NOTE: STRATIFIED ABUNDANCES "
                "ARE BASED ON HOW MUCH THAT PREDICTED GENOME (E.G. SEQUENCE) "
                "CONTRIBUTES TO THE COMMUNITY-WIDE ABUNDANCE, NOT THE "
                "ABUNDANCE OF THE PATHWAY BASED ON THE PREDICTED GENES IN "
                "THAT GENOME ALONE. In other words, a predicted genome might "
                "be contributing a lot to the community-wide pathway "
                "abundance simply because one required gene for that pathway "
                "is at extremely high abundance in that genome even though no "
                "other required genes for that pathway are present.",

    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-i', '--input', metavar='IN_TABLE', required=True,
                    type=str,
                    help='Input TSV table of gene family abundances')

parser.add_argument('-m', '--map', metavar='MAP', required=True, type=str,
                    help='MinPath mapfile')

parser.add_argument('-o', '--out_dir', metavar='DIRECTORY', required=True,
                    type=str, help='Output folder for pathway abundance output')

parser.add_argument('--intermediate', metavar='DIR', type=str, default=None,
                    help='Output folder for intermediate files (wont be ' +
                         'kept unless this option is set.')

parser.add_argument('-p', '--proc', default=1, type=int,
                    help='Number of processes to run (default: %(default)d).')

parser.add_argument('--print_cmds', default=False, action="store_true",
                    help='If specified, print out wrapped commands to screen')

parser.add_argument('-v', '--version', default=False, action='version',
                    version="%(prog)s " + __version__)

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

    make_output_dir(args.out_dir)

    # Write output files.
    unstrat_outfile = path.join(args.out_dir, "path_abun_unstrat.tsv")
    unstrat_out.to_csv(path_or_buf=unstrat_outfile,  sep="\t",
                       index_label="pathway")

    # Write stratified output only if something besides None was returned.
    if strat_out is not None:
        strat_outfile = path.join(args.out_dir, "path_abun_strat.tsv")
        strat_out.to_csv(path_or_buf=strat_outfile,  sep="\t", index=False)


if __name__ == "__main__":
    main()
