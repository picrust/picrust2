#!/usr/bin/env python

__copyright__ = 'Copyright 2018, The PICRUSt Project'
__license__ = 'GPL'
__version__ = '2.0.0-b.6'

import argparse
from picrust2.util import (check_files_exist, convert_humann2_to_picrust2,
                           convert_picrust2_to_humann2,
                           convert_picrust2_to_humann2_merged)

CONVERSION_CHOICES = ['humann2_strat_to_picrust2',
                      'humann2_unstrat_to_picrust2', 
                      'picrust2_unstrat_to_humann2_split',
                      'picrust2_strat_to_humann2_split',
                      'picrust2_to_humann2_merged']

parser = argparse.ArgumentParser(

    description=
    'Converts to and from PICRUSt2 function abundance table. Currently '
    'only converting to and from HUMAnN2 function tables are supported. Both '
    'stratified and unstratified tables are supported. Note that the '
    'categories like \"UNMAPPED\" in the HUMAnN2 tables will be removed '
    'if they have values of 0.',

    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('input_files', metavar='INPUT', type=str, nargs='+',
                    help='Input table to convert. If there are multiple input '
                    'files (e.g. if multiple HUMAnN2 gene tables for different '
                    'samples should be converted to a single PICRUSt2 table) '
                    'then specify them all: file1 file2 file3...')

parser.add_argument('-o', '--output', metavar='OUTPUT', required=True,
                    type=str,
                    help='Path to output. Corresponds to folder name if '
                         'multiple files are output, otherwise it will be a '
                         'filename.')

parser.add_argument('-c', '--conversion', metavar='CONVERSION', required=True,
                    choices=CONVERSION_CHOICES,
                    help='Type of conversion to perform '
                         '(\'humann2_unstrat_to_picrust2\', '
                         '\'humann2_strat_to_picrust2\', '
                         '\'picrust2_unstrat_to_humann2_split\', '
                         '\'picrust2_strat_to_humann2_split\', or '
                         '\'picrust2_to_humann2_merged\').')

parser.add_argument('-v', '--version', default=False, action='version',
                    version='%(prog)s ' + __version__)

def main():

    args = parser.parse_args()

    # Check that input files exist.
    check_files_exist(args.input_files)

    if args.conversion == 'humann2_unstrat_to_picrust2':
        convert_humann2_to_picrust2(args.input_files, args.output, False)

    elif args.conversion == 'humann2_strat_to_picrust2':
        convert_humann2_to_picrust2(args.input_files, args.output, True)

    elif args.conversion == 'picrust2_unstrat_to_humann2_split':

        convert_picrust2_to_humann2(args.input_files, args.output, False)

    elif args.conversion == 'picrust2_strat_to_humann2_split':

        convert_picrust2_to_humann2(args.input_files, args.output, True)

    elif args.conversion == 'picrust2_to_humann2_merged':

        convert_picrust2_to_humann2_merged(args.input_files, args.output)


if __name__ == '__main__':
    main()
