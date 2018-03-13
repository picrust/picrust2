#!/usr/bin/env python

from __future__ import division

__author__ = "Gavin Douglas"
__copyright__ = "Copyright 2011-2017, The PICRUSt Project"
__credits__ = ["Gavin Douglas", "Morgan Langille"]
__license__ = "GPL"
__version__ = "2-alpha.1"
__maintainer__ = "Gavin Douglas"
__email__ = "gavinmdouglas@gmail.com"
__status__ = "Development"


from cogent.util.option_parsing import parse_command_line_parameters, make_option
from biom.table import Table
from biom import load_table
import tempfile
from picrust.util import system_call_check, make_output_dir
from picrust.run_minpath import pathway_counts, minpath_wrapper
from joblib import Parallel, delayed
from pprint import pprint

script_info = {}

script_info['script_description'] = "Runs MinPath on a table of E.C. numbers"

# Define command-line interface
script_info['output_description'] = "Output is a tab-delimited table of " +\
                                   "predicted pathway abundances"

script_info['required_options'] = [

  make_option('-i', '--input', type="existing_filepath",
              help='the input biom table of gene family abundances'),

  make_option('-m', '--map', type="existing_filepath",
              help='path to MinPath map'),

  make_option('-o', '--output', type="new_filepath",

              help='the output filepath for pathway abundances')
]

script_info['optional_options'] = [

  make_option('--keep_tmp', default=False, action="store_true",
              help='if specified, keep temporary folder ' +
                   '[default: %default]'),

  make_option('-t', '--threads', default=1, type="int",
              help='Number of threads [default: %default]'),

  make_option('--tmp_dir', type="new_filepath",
              help='temporary directory for running MinPath'),

  make_option('--print_cmds', default=False, action="store_true",
              help='if specified, print out wrapped commands to screen ' +
                   '[default: %default]')
]

script_info['version'] = __version__


def main():

    option_parser, opts, args = parse_command_line_parameters(**script_info)

    # Create temporary folder for intermediate files.
    if opts.tmp_dir:
        tmp_dir = opts.tmp_dir
    else:
        tmp_dir = "minpath_tmp_" + next(tempfile._get_candidate_names())

    make_output_dir(tmp_dir)

    if opts.print_cmds:
      print("Creating tmp directory: " + tmp_dir)

    biom_in = load_table(opts.input)

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

    sample_path_abun_raw = Parallel(n_jobs=opts.threads)(delayed(
                                minpath_wrapper)(sample_id, biom_in,
                                opts.map, tmp_dir, functions, opts.print_cmds)
                                for sample_id in samples)

    # Figure out what all unique pathway names are.
    all_pathways = []
    for sample_d in sample_path_abun_raw:
      print(sample_d)
      all_pathways = list(set(all_pathways + list(sample_d.keys())))
    all_pathways = set(all_pathways)

    # Loop through all samples and make dictionary of these return pathway
    # abundances.
    sample_path_abun = {}
    for i, sample_id in enumerate(samples):
      sample_path_abun[sample_id] = sample_path_abun_raw[i]

    pprint(sample_path_abun)

    # Write output file of pathway abundances.
    outfile = open(opts.output, "w")

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
    if not opts.keep_tmp:
        system_call_check("rm -r " + tmp_dir, print_out=opts.print_cmds)

if __name__ == "__main__":
    main()
