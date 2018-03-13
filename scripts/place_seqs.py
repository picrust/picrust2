#!/usr/bin/env python

from __future__ import division

__author__ = "Gavin Douglas"
__copyright__ = "Copyright 2011-2017, The PICRUSt Project"
__credits__ = ["Gavin Douglas", "Morgan Langille"]
__license__ = "GPL"
__version__ = "2-alpha.2"
__maintainer__ = "Gavin Douglas"
__email__ = "gavinmdouglas@gmail.com"
__status__ = "Development"


from cogent.util.option_parsing import parse_command_line_parameters, make_option
import tempfile
from picrust.util import system_call_check, make_output_dir, read_fasta, read_phylip, write_fasta, write_phylip

script_info = {}
script_info['brief_description'] = "Will place study seqs in phylogeny with " +\
                                   "reference sequences."

script_info['script_description'] = "Wrapper to prep tree before HSP steps. " +\
                                    "Requires tree of reference sequences, " +\
                                    "FASTA MSA of reference seqs, and " +\
                                    "unaligned FASTA of study sequences."

# Define command-line interface.
script_info['output_description'] = "Output is a tree of the reference and " +\
                                    "study sequences"
script_info['required_options'] = [

  make_option('-r', '--ref_msa', type="existing_filepath",
              help='FASTA of aligned reference sequences'),

  make_option('-s', '--study_fasta', type="existing_filepath",
              help='FASTA of unaligned study sequences'),

  make_option('-t', '--tree', type="existing_filepath",
              help='Input tree based on aligned reference sequences.'),

  make_option('-o', '--out_tree', type="new_filepath",
              help='Name of final output tree')

]

script_info['optional_options'] = [

  make_option('--keep_tmp', default=False, action="store_true",
              help='if specified, keep temporary folder ' +
                   '[default: %default]'),

  make_option('--threads', type="int", default=1,
              help='Number of threads to use when possible. ' +
                   '[default: %default]'),

  make_option('--papara_output', type="existing_filepath", default=None,
              help='Path to papara output in Phylip format ' +
                   '(will skip papara step)'),

  make_option('--tmp_dir', type="new_filepath", default=None,
              help='Temporary folder for intermediate files.'),

  make_option('--chunk_size', type="int", default=5000,
              help='Number of query seqs to read in at once for epa-ng. ' +
                   '[default: %default]'),

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
      tmp_dir = "place_seqs_tmp_" + next(tempfile._get_candidate_names())

    make_output_dir(tmp_dir)

    # Define temporary filenames.
    ref_phylip_out = tmp_dir + "/ref_seqs.phylip"

    papara_suffix = next(tempfile._get_candidate_names())

    papara_filename = "papara_alignment." + papara_suffix

    papara_fasta_out = tmp_dir + "/" + papara_filename + ".fna"

    study_fasta_aligned = tmp_dir + "/study_seqs_papara.fasta"

    ref_fasta_aligned = tmp_dir + "/ref_seqs_papara.fasta"

    epa_out_dir = tmp_dir + "/epa_out"

    # Read in ref seqs FASTA.
    ref_msa = read_fasta(opts.ref_msa)

    # Either read in papara output or run it.
    if opts.papara_output:
      # Read in papara output if already done.
      papara_out = read_phylip(opts.papara_output, check_input=True)

    else:
      # Convert ref sequences from MSA FASTA to phylip.
      write_phylip(ref_msa, ref_phylip_out)

      # Run papara to align query seqs to ref sequences.
      system_call_check("papara -t " + opts.tree + " -s " + ref_phylip_out +
                        " -q " + opts.study_fasta + " -j " + str(opts.threads) +
                        " -n " + papara_suffix, print_out=opts.print_cmds)

      # Read in papara phylip output.
      papara_out = read_phylip(papara_filename, check_input=True)

      # Move papara Phylip file, quality, and log to tmp folder.
      system_call_check("mv " + papara_filename + " " + tmp_dir,
                        print_out=opts.print_cmds)
      system_call_check("mv papara_log." + papara_suffix + " " + tmp_dir,
                        print_out=opts.print_cmds)
      system_call_check("mv papara_quality." + papara_suffix + " " + tmp_dir,
                        print_out=opts.print_cmds)

    # Split papara phylip output into FASTA MSA files of study sequences and
    # reference sequences separately.
    all_seq_names = set(list(papara_out.keys()))
    ref_seq_names = set(list(ref_msa.keys()))
    study_seq_names = all_seq_names.difference(ref_seq_names)

    ref_papara_subset = {seq: papara_out[seq] for seq in ref_seq_names}
    study_papara_subset = {seq: papara_out[seq] for seq in study_seq_names}

    write_fasta(ref_papara_subset, ref_fasta_aligned)
    write_fasta(study_papara_subset, study_fasta_aligned)

    # Run EPA (output needs to be to a directory, which will be created).
    make_output_dir(epa_out_dir)

    system_call_check("epa-ng --tree " + opts.tree + " --ref-msa " +
                      ref_fasta_aligned + " --query " + study_fasta_aligned +
                      " --chunk-size " + str(opts.chunk_size) + " -T " +
                      str(opts.threads) + " -w " + epa_out_dir,
                      print_out=opts.print_cmds)

    system_call_check("guppy tog " + epa_out_dir + "/epa_result.jplace" +
                      " -o " + opts.out_tree, print_out=opts.print_cmds)

    # Remove intermediate files unless "--keep_tmp" option specified.
    if not opts.keep_tmp:
        system_call_check("rm -r " + tmp_dir, print_out=opts.print_cmds)


if __name__ == "__main__":
    main()
