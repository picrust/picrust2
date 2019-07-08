#!/usr/bin/env python

__copyright__ = "Copyright 2018-2019, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.1.4-b"

import argparse
from os import path
from picrust2.metagenome_pipeline import run_metagenome_pipeline
from picrust2.util import check_files_exist

parser = argparse.ArgumentParser(

    description="Per-sample metagenome functional profiles are generated " +
                "based on the predicted functions for each study sequence. " +
                "Note that typically these sequences correspond to OTUs or " +
                "ASVs. The specified sequence abundance table will be " +
                "normalized by the predicted number of marker gene copies " +
                "before outputting the final files. Two main output files " +
                "will be generated: one stratified and one non-stratified " +
                "by contributing taxa",

    epilog='''Usage example:
metagenome_pipeline.py -i seqabun.biom -f predicted_EC.tsv.gz -m predicted_16S.tsv.gz --max_nsti 2.0 -o metagenome_out
''',
    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-i', '--input', metavar='PATH',
                    required=True, type=str,
                    help='Input table of sequence abundances (BIOM, TSV, or ' +
                         'mothur shared file format).')

parser.add_argument('-f', '--function', metavar='PATH',
                    required=True, type=str,
                    help='Table of predicted gene family copy numbers ' +
                         '(output of hsp.py).')

parser.add_argument('-m', '--marker', metavar='PATH',
                    required=True, type=str,
                    help='Table of predicted marker gene copy numbers ' +
                         '(output of hsp.py, typically for 16S).')

parser.add_argument('--max_nsti', metavar='FLOAT', type=float, default=2.0,
                    help='Sequences with NSTI values above this value will ' +
                         'be excluded (default: %(default)d).')

parser.add_argument('--min_reads', metavar='INT', type=int, default=1,
                    help='Minimum number of reads across all samples for ' +
                         'each input ASV. ASVs below this cut-off will be ' +
                         'counted as part of the \"RARE\" category in the ' +
                         'stratified output (default: %(default)d).')

parser.add_argument('--min_samples', metavar='INT', type=int, default=1,
                    help='Minimum number of samples that an ASV needs to be ' +
                         'identfied within. ASVs below this cut-off will be ' +
                         'counted as part of the \"RARE\" category in the ' +
                         'stratified output (default: %(default)d).')

parser.add_argument('--strat_out', default=False, action='store_true',
                    help='Output table stratified by sequences as well. By '
                         'default this will be in \"contributional\" format '
                         '(i.e. long-format) unless the \"--wide_table\" '
                         'option is set. The startified outfile is named '
                         '\"metagenome_contrib.tsv.gz\" when in long-format.')

parser.add_argument('--wide_table', default=False, action='store_true',
                    help='Output wide-format stratified table of metagenome '
                         'predictions when \"--strat_out\" is set. This is '
                         'the deprecated method of generating stratified '
                         'tables since it is extremely memory intensive. The '
                         'startified outfile is named '
                         '\"pred_metagenome_strat.tsv.gz\" when this option '
                         'is set.')

parser.add_argument('-o', '--out_dir', metavar='PATH', type=str,
                    default='metagenome_out',
                    help='Output directory for metagenome predictions. ' +
                         '(default: %(default)s).')

parser.add_argument('-v', '--version', default=False, action='version',
                    version="%(prog)s " + __version__)


def main():

    args = parser.parse_args()

    check_files_exist([args.input, args.function, args.marker])

    strat_pred, unstrat_pred = run_metagenome_pipeline(
                                            input_seqabun=args.input,
                                            function=args.function,
                                            marker=args.marker,
                                            out_dir=args.out_dir,
                                            max_nsti=args.max_nsti,
                                            min_reads=args.min_reads,
                                            min_samples=args.min_samples,
                                            strat_out=args.strat_out,
                                            wide_table=args.wide_table)

    unstrat_outfile = path.join(args.out_dir, "pred_metagenome_unstrat.tsv.gz")
    unstrat_pred.to_csv(path_or_buf=unstrat_outfile, sep="\t", index=True,
                        index_label="function", compression="gzip")

    if args.strat_out and not args.wide_table:
        strat_outfile = path.join(args.out_dir, "pred_metagenome_contrib.tsv.gz")
        strat_pred.to_csv(path_or_buf=strat_outfile, sep="\t", index=False,
                          compression="gzip")

    elif args.strat_out and args.wide_table:
        strat_outfile = path.join(args.out_dir, "pred_metagenome_strat.tsv.gz")
        strat_pred.to_csv(path_or_buf=strat_outfile, sep="\t", index=True,
                          compression="gzip")


if __name__ == "__main__":
    main()
