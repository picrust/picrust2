#!/usr/bin/env python

__copyright__ = "Copyright 2018-2020, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.2.1-b"

import argparse
from picrust2.default import default_map
from picrust2.util import add_descrip_col, make_output_dir_for_file
import sys

TRAIT_OPTIONS = ['METACYC', 'COG', 'EC', 'KO', 'PFAM', 'TIGRFAM']

parser = argparse.ArgumentParser(

    description="This script adds a description column to a function " +
                "abundance table and outputs a new file. The user needs " +
                "to specify the input file and what type of functions are " +
                "in the input table. Will throw an error if no ids overlap " +
                "and otherwise will fill in \"not_found\" for the " +
                "description of ids in the function table not in the " +
                "mapfile.",
    epilog='''

Usage:
add_descriptions.py -i IN_TABLE -m KO -o OUT_TABLE

''', formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-i', '--input', metavar='PATH', required=True, type=str,
                    help='Input function abundance table.')

parser.add_argument('-o', '--output', metavar='PATH', type=str, required=True,
                    help='Output function abundance table with added '
                         'description column. If the extension \".gz\" '
                         'is added the table will automatically be gzipped.')

parser.add_argument('-m', '--map_type', type=str.upper, choices=TRAIT_OPTIONS,
                    help='Specifies which default mapping table should be '
                          'used. Use the --custom_map_table option '
                          'to input a non-default mapping table.')

parser.add_argument('--custom_map_table', metavar='PATH', type=str,
                    help='An input map table linking function ids to '
                         'descriptions for each function. ')

parser.add_argument('-v', '--version', default=False, action='version',
                    version="%(prog)s " + __version__)


def main():

    args = parser.parse_args()

    # Determine which mapping type was specified. If neither a default
    # or custom mapping was specified then throw an error.
    if args.map_type and args.custom_map_table:
        sys.exit("Only one of \"--map_type\" or \"--custom_map_table\" can be "
                 "set. Please re-run the command with only one of these "
                 "options.")
    elif args.map_type:
        mapfile = default_map[args.map_type]
    elif args.custom_map_table:
        mapfile = args.custom_map_table
    else:
        sys.exit("A default mapping table needs to be specified with the "
                 "--map_type option, or alternatively a custom mapfile can "
                 "be specified with the --custom_map_table option")

    tab_w_descrip = add_descrip_col(inputfile=args.input, mapfile=mapfile)

    # Output the table to file.
    make_output_dir_for_file(args.output)
    tab_w_descrip.to_csv(path_or_buf=args.output, sep="\t",
                         index=False, compression="infer")


if __name__ == "__main__":
    main()
