#!/usr/bin/env python

from __future__ import division

__author__ = "Gavin Douglas"
__copyright__ = "Copyright 2018, The PICRUSt Project"
__credits__ = ["Gavin Douglas", "Morgan Langille"]
__license__ = "GPL"
__version__ = "2-alpha.6"

import argparse
from biom.table import Table
from biom import load_table
import tempfile
from picrust.util import system_call_check, make_output_dir
from picrust.run_minpath import pathway_counts, minpath_wrapper
from joblib import Parallel, delayed

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

    # Create temporary folder for intermediate files.
    if args.tmp_dir:
        tmp_dir = args.tmp_dir
    else:
        tmp_dir = "minpath_tmp_" + next(tempfile._get_candidate_names())

    make_output_dir(tmp_dir)

    if args.print_cmds:
      print("Creating tmp directory: " + tmp_dir)

    biom_in = load_table(args.input)

    # Remove all empty rows and columns.
    biom_in.remove_empty(axis='whole', inplace=True)

    samples = biom_in.ids()
    functions = biom_in.ids(axis="observation")

    # Initialize set of all pathways.
    all_pathways = set()

    functions_map = {}

    for func in functions:
      functions_map[func] = func.replace('EC:', '')

    biom_in.update_ids(functions_map,
                       axis='observation',
                       strict=True,
                       inplace=True)

    functions = list(functions_map.values())

    sample_path_abun_raw = Parallel(n_jobs=args.threads)(delayed(
                                minpath_wrapper)(sample_id, biom_in,
                                args.map, tmp_dir, functions, args.print_cmds)
                                for sample_id in samples)

    # Figure out what all unique pathway names are.
    all_pathways = []
    for sample_d in sample_path_abun_raw:
      all_pathways = list(set(all_pathways + list(sample_d.keys())))
    all_pathways = set(all_pathways)

    # Loop through all samples and make dictionary of these return pathway
    # abundances.
    sample_path_abun = {}
    for i, sample_id in enumerate(samples):
      sample_path_abun[sample_id] = sample_path_abun_raw[i]

    # Write output file of pathway abundances.
    outfile = open(args.output, "w")

    # Write header-line.
    outfile.write("\t".join(["pathway"] + list(samples)) + "\n")

    # Loop through pathways and write out abundances per sample.
    for pathway in all_pathways:
      out_row = [pathway]
      for sample_id in samples:
        out_row += [str(sample_path_abun[sample_id][pathway])]

      outfile.write("\t".join(out_row) + "\n")

    outfile.close()

    # Remove intermediate files unless "--keep_tmp" option specified.
    if not args.keep_tmp:
        system_call_check("rm -r " + tmp_dir, print_out=args.print_cmds)

if __name__ == "__main__":
    main()
