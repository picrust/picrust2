#!/usr/bin/env python

__author__ = "Gavin Douglas"
__copyright__ = "Copyright 2018, The PICRUSt Project"
__credits__ = ["Gavin Douglas", "Morgan Langille"]
__license__ = "GPL"
__version__ = "2-alpha.3"

from cogent.util.option_parsing import parse_command_line_parameters, make_option
from picrust.util import system_call_check

script_info = {}
script_info['brief_description'] = "Will run full metagenome prediction " +\
                                   "pipeline."

script_info['script_description'] = "Wrapper for the three metagenome " +\
                                    "prediction pipeline scripts"

# Define command-line interface
script_info['output_description'] = "Output is a tab-delimited table of " +\
                                   "predicted character states"
script_info['required_options'] = [

  make_option('-i', '--input', type="existing_filepath",
              help='input biom table'),

  make_option('-c', '--input_copies', type="existing_filepath",
              help='Table of predicted marker gene copy numbers'),

  make_option('-f', '--input_function', type="existing_filepath",
              help='Table of predicted gene family copy numbers')

]

script_info['optional_options'] = [

  make_option('-o', '--out_prefix', type="new_filepath",
              default='pipeline_out',
              help='prefix for output file names ' +
                   '[default: %default]'),

  make_option('--tsv', default=False, action="store_true",
              help='if specified, also output tables in TSV format ' +
                   '[default: %default]')
]

script_info['version'] = __version__


def main():

    option_parser, opts, args = parse_command_line_parameters(**script_info)

    norm_out = opts.out_prefix + ".norm.biom"
    meta_out = opts.out_prefix + ".genefamilies.biom"

    norm_cmd = "normalize_by_copy_number.py -i " + opts.input + " -c " +\
               opts.input_copies + " -o " + norm_out

    meta_cmd = "predict_metagenomes.py -i " + norm_out + " -c " +\
               opts.input_function + " -o " + meta_out

    print(norm_cmd)
    process = system_call_check(norm_cmd.split(" "))

    print(meta_cmd)
    process = system_call_check(meta_cmd.split(" "))

    if opts.tsv:
        norm_out_tsv = opts.out_prefix + ".norm.biom.tsv"
        meta_out_tsv = opts.out_prefix + ".genefamilies.biom.tsv"

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
