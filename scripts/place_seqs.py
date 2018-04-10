#!/usr/bin/env python

__license__ = "GPL"
__version__ = "2-alpha.6"

import argparse
import tempfile
from picrust2.util import (system_call_check, make_output_dir, read_fasta,
                           read_phylip, write_fasta, write_phylip,
                           get_picrust_project_dir)

project_dir = get_picrust_project_dir()

default_fasta = project_dir + \
                "/precalculated/prokaryotic/img_centroid_16S_aligned.fna"
default_tree = project_dir + \
               "/precalculated/prokaryotic/img_centroid_16S_aligned.tree"

parser = argparse.ArgumentParser(

    description="Wrapper to prep tree before HSP steps. Requires unaligned " +
                "FASTA of study sequences. Users can specify a non-default " +
                "reference fasta and treefile if needed.",

    formatter_class=argparse.RawDescriptionHelpFormatter)


parser.add_argument('-s', '--study_fasta', metavar='PATH', required=True,
                    type=str, help='FASTA of unaligned study sequences')

parser.add_argument('-r', '--ref_msa', metavar='PATH', type=str,
                    default=default_fasta,
                    help='FASTA of aligned reference sequences')

parser.add_argument('-t', '--tree', metavar='PATH', type=str,
                    default=default_tree,
                    help='Input tree based on aligned reference sequences.')

parser.add_argument('-o', '--out_tree', metavar='PATH', required=True,
                    type=str, help='Name of final output tree')

parser.add_argument('--keep_tmp', default=False, action='store_true',
                    help='If specified, keep temporary folder')

parser.add_argument('--threads', type=int, default=1,
                    help='Number of threads to use when possible')

parser.add_argument('--papara_output', metavar='PATH', type=str, default=None,
                    help='Path to PaPaRa output in Phylip format (will skip ' +
                         'PaPaRa step)')

parser.add_argument('--tmp_dir', metavar='PATH', type=str, default=None,
                    help='Temporary folder for intermediate files')

parser.add_argument('--chunk_size', type=int, default=5000,
                    help='Number of query seqs to read in at once for epa-ng')

parser.add_argument('--print_cmds', default=False, action='store_true',
                    help='If specified, print out wrapped commands to screen')


def main():

    args = parser.parse_args()

    # Create temporary folder for intermediate files.
    if args.tmp_dir:
      tmp_dir = args.tmp_dir
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
    ref_msa = read_fasta(args.ref_msa)

    # Either read in papara output or run it.
    if args.papara_output:
      # Read in papara output if already done.
      papara_out = read_phylip(args.papara_output, check_input=True)

    else:
      # Convert ref sequences from MSA FASTA to phylip.
      write_phylip(ref_msa, ref_phylip_out)

      # Run papara to align query seqs to ref sequences.
      system_call_check("papara -t " + args.tree + " -s " + ref_phylip_out +
                        " -q " + args.study_fasta + " -j " + str(args.threads) +
                        " -n " + papara_suffix, print_out=args.print_cmds)

      # Read in papara phylip output.
      papara_out = read_phylip(papara_filename, check_input=True)

      # Move papara Phylip file, quality, and log to tmp folder.
      system_call_check("mv " + papara_filename + " " + tmp_dir,
                        print_out=args.print_cmds)
      system_call_check("mv papara_log." + papara_suffix + " " + tmp_dir,
                        print_out=args.print_cmds)
      system_call_check("mv papara_quality." + papara_suffix + " " + tmp_dir,
                        print_out=args.print_cmds)

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

    system_call_check("epa-ng --tree " + args.tree + " --ref-msa " +
                      ref_fasta_aligned + " --query " + study_fasta_aligned +
                      " --chunk-size " + str(args.chunk_size) + " -T " +
                      str(args.threads) + " -w " + epa_out_dir,
                      print_out=args.print_cmds)

    system_call_check("guppy tog " + epa_out_dir + "/epa_result.jplace" +
                      " -o " + args.out_tree, print_out=args.print_cmds)

    # Remove intermediate files unless "--keep_tmp" option specified.
    if not args.keep_tmp:
        system_call_check("rm -r " + tmp_dir, print_out=args.print_cmds)


if __name__ == "__main__":
    main()
