#!/usr/bin/env python

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.0.0-b.1"

import argparse
from os import path
from picrust2.metagenome_pipeline import run_metagenome_pipeline
from picrust2.util import check_files_exist

parser = argparse.ArgumentParser(

    description="Per-sample metagenome functional profiles are generated " +
                "based on the predicted functions for each study sequence. " +
                "Note that typically these sequences correspond to OTUs or " +
                "ASVs. The specified sequence abundance table will be normalized " +
                "by the predicted number of marker gene copies. Two main " +
                "output files will be generated: one stratified and one " +
                "non-stratified by contributing taxa",

    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-i', '--input', metavar='PATH',
                    required=True, type=str,
                    help='Input table of sequence abundances (BIOM or TSV ' +
                         'format)')

parser.add_argument('-f', '--function', metavar='PATH',
                    required=True, type=str,
                    help='Table of predicted gene family copy numbers ' +
                         '(output of hsp.py)')

parser.add_argument('-m', '--marker', metavar='PATH',
                    required=True, type=str,
                    help='Table of predicted marker gene copy numbers ' +
                         '(output of hsp.py, typically for 16S)')

parser.add_argument('--max_nsti', metavar='INT', type=int, default=2,
                    help='Sequences with NSTI values above this value will ' +
                         'be excluded (default: %(default)d).')

parser.add_argument('-p', '--proc', metavar='INT', type=int, default=1,
                    help='Number of processes to run in parallel ' +
                         '(default: %(default)d).')

parser.add_argument('-o', '--out_dir', metavar='PATH', type=str,
                    default='metagenome_out',
                    help='Output directory for metagenome predictions. ' +
                         '(default: %(default)s).')

parser.add_argument('-v', '--version', default=False, action='version',
                    version="%(prog)s " + __version__)

def main():

    args = parser.parse_args()

    # Check that input files exist.
    check_files_exist([args.input, args.function, args.marker])

    # Pass arguments to key function and get predicted functions
    # stratified and unstratified by genomes.
    strat_pred, unstrat_pred = run_metagenome_pipeline(input_biom=args.input,
                                                       function=args.function,
                                                       marker=args.marker,
                                                       out_dir=args.out_dir,
                                                       max_nsti=args.max_nsti,
                                                       proc=args.proc,
                                                       output_normfile=True)

    # Generate output table filepaths and write out pandas dataframes.
    strat_outfile = path.join(args.out_dir, "pred_metagenome_strat.tsv")
    unstrat_outfile = path.join(args.out_dir, "pred_metagenome_unstrat.tsv")

    # Note that no index labels are written for stratified output.
    strat_pred.to_csv(path_or_buf=strat_outfile, sep="\t", index=False)
    unstrat_pred.to_csv(path_or_buf=unstrat_outfile, sep="\t")


if __name__ == "__main__":
    main()
