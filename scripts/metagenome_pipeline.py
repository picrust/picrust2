#!/usr/bin/env python

__license__ = "GPL"
__version__ = "2-alpha.7"

import argparse
from picrust2.util import system_call_check

parser = argparse.ArgumentParser(

    description="Wrapper for the three metagenome prediction pipeline scripts",

    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-i', '--input', metavar='PATH',
                    required=True, type=str,
                    help='Input BIOM table')

parser.add_argument('-c', '--input_copies', metavar='PATH',
                    required=True, type=str,
                    help='Table of predicted marker gene copy numbers')

parser.add_argument('-f', '--input_function', metavar='PATH',
                    required=True, type=str,
                    help='Table of predicted gene family copy numbers')

parser.add_argument('-o', '--out_prefix', metavar='PREFIX', type=str,
                    default='pipeline_out',
                    help='Prefix for output file names')

parser.add_argument('--tsv', action="store_true",
                    help='If specified, also output tables in TSV format')


def main():

    args = parser.parse_args()

    norm_out = args.out_prefix + ".norm.biom"
    meta_out = args.out_prefix + ".genefamilies.biom"

    norm_cmd = "normalize_by_copy_number.py -i " + args.input + " -c " +\
               args.input_copies + " -o " + norm_out

    meta_cmd = "predict_metagenomes.py -i " + norm_out + " -c " +\
               args.input_function + " -o " + meta_out

    print(norm_cmd)
    process = system_call_check(norm_cmd.split(" "))

    print(meta_cmd)
    process = system_call_check(meta_cmd.split(" "))

    if args.tsv:
        norm_out_tsv = args.out_prefix + ".norm.biom.tsv"
        meta_out_tsv = args.out_prefix + ".genefamilies.biom.tsv"

        norm_convert_cmd = "biom convert -i " + norm_out + " -o " +\
                           norm_out_tsv + " --to-tsv"

        meta_convert_cmd = "biom convert -i " + meta_out + " -o " +\
                           meta_out_tsv + " --to-tsv"

        print(norm_convert_cmd)
        process = system_call_check(norm_convert_cmd.split(" "))

        print(meta_convert_cmd)
        process = system_call_check(meta_convert_cmd.split(" "))


if __name__ == "__main__":
    main()
