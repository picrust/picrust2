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
from picrust.run_minpath import harmonic_mean, pathway_counts

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

    # Initialize dictionary which will contain pathway abundances per sample.
    sample_path_abun = {}

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

    for sample_id in samples:

        # Define MinPath input and outout filenames.
        minpath_in = str(tmp_dir + "/" + sample_id + "_minpath_in.txt")
        minpath_report = str(tmp_dir + "/" + sample_id + "_minpath_report.txt")
        minpath_details = str(tmp_dir + "/" + sample_id + "_minpath_details.txt")

        id_minpath_fh = open(minpath_in, "w")

        # Counter to give each "read" in MinPath input a different id.
        func_num = 0

        for func_id in functions:

          # Get count of each sequence in sample and write that sequence out
          # each time.
          func_count = int(biom_in.get_value_by_ids(obs_id=func_id,
                                                          samp_id=sample_id))
          # If 0 then skip.
          if func_count == 0:
            continue

          for i in range(func_count):
            func_num += 1
            func_name = "seq_" + str(func_num)
            id_minpath_fh.write(func_name + "\t" + func_id + "\n")

        id_minpath_fh.close()

        # Run MinPath on this sample.
        minpath_cmd = "MinPath1.4.py -any " + minpath_in + " -map " +\
                      opts.map + " -report " + minpath_report +\
                      " -details " + minpath_details

        system_call_check(minpath_cmd, print_out=opts.print_cmds)

        # Read through MinPath report and keep track of pathways identified
        # to be present.
        path_present = set()

        with open(minpath_report, "r") as minpath_report_in:
          for line in minpath_report_in:
            line_split = line.split()

            if int(line_split[7]) == 1:
              path_present.add(line_split[-1])

        # Now read in details file and take abundance of pathway to be
        # harmonic mean of gene families in pathway to be abundance of pathway.
        # Abundances of 0 will be added in for gene families not found.

        # Initialize pathway_counts instance for sample.
        sample_path_abun[sample_id] = pathway_counts(sample_id)

        # Initialize dictionary that will contain gene family abundance per
        # pathway.
        gf_abundances = {}

        # Boolean specifying that pathway in details file was called as
        # present by MinPath.
        present = False

        with open(minpath_details, "r") as minpath_details_in:
            for line in minpath_details_in:
              line_split = line.split()

              # If line starts with "path" then keep track of pathway name if
              # it was called as present in report file.
              if line_split[0] == "path":
                if line_split[-1] not in path_present:
                  present = False
                  continue
            
                present = True

                current_pathway = line_split[-1]

                # Initialize list containing gene family abundances.
                gf_abundances[current_pathway] = []

                # Add in abundances of 1 for missing genes.
                for i in range(int(line_split[3]) - int(line_split[5])):
                  gf_abundances[current_pathway] += [1]

            # If line does not start with "path" then only proceed if current
            # pathway is present.
              elif present:
                gf_abundances[current_pathway] += [int(line_split[2])]

        # Loop through all pathways present and get harmonic mean.
        for pathway in gf_abundances.keys():
          all_pathways.add(pathway)

          sample_path_abun[sample_id].set_pathway_abun(pathway,
                                                       harmonic_mean(gf_abundances[pathway]))

    # Write output file of pathway abundances.
    outfile = open(opts.output, "w")

    # Write header-line.
    outfile.write("\t".join(["pathway"] + list(samples)) + "\n")

    # Loop through pathways and write out abundances per sample.
    for pathway in all_pathways:
      out_row = [pathway]
      for sample_id in samples:
        out_row += [str(sample_path_abun[sample_id].return_pathway_abun(pathway))]

      outfile.write("\t".join(out_row) + "\n")

    outfile.close()

    # Remove intermediate files unless "--keep_tmp" option specified.
    if not opts.keep_tmp:
        system_call_check("rm -r " + tmp_dir, print_out=opts.print_cmds)

if __name__ == "__main__":
    main()
