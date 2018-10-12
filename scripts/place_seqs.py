#!/usr/bin/env python

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.0.2-b"

import argparse
from tempfile import TemporaryDirectory
from picrust2.place_seqs import place_seqs_pipeline
from picrust2.default import default_fasta, default_tree, default_hmm
from picrust2.util import make_output_dir, check_files_exist

ALIGN_CHOICES = ["hmmalign", "papara"]

parser = argparse.ArgumentParser(

    description="Wrapper to prep tree before HSP steps. Requires unaligned " +
                "FASTA of study sequences. Users can specify a non-default " +
                "reference fasta and treefile if needed. Note that typically " +
                "the input study sequences are representive sequences of " +
                "operational taxonomic units or amplicon sequence variants.",

    formatter_class=argparse.RawDescriptionHelpFormatter)


parser.add_argument('-s', '--study_fasta', metavar='PATH', required=True,
                    type=str, help='FASTA of unaligned study sequences')

parser.add_argument('-r', '--ref_msa', metavar='PATH', type=str,
                    default=default_fasta,
                    help='FASTA of aligned reference sequences (default: %(default)s).')

parser.add_argument('--hmm', metavar='PATH', type=str,
                    default=default_hmm,
                    help='Hidden markov model of reference MSA (default: %(default)s).')

parser.add_argument('-t', '--tree', metavar='PATH', type=str,
                    default=default_tree,
                    help='Input tree based on aligned reference sequences. ' +
                         '(default: %(default)s).')

parser.add_argument('-o', '--out_tree', metavar='PATH', required=True,
                    type=str, help='Name of final output tree')

parser.add_argument('-a', '--alignment_tool', type=str.lower,
                    default="hmmalign", choices=ALIGN_CHOICES,
                    help='Which program to use for aligning query sequences ' +
                         'to reference MSA prior to EPA-NG step (default: ' +
                         '%(default)s).')

parser.add_argument('--threads', type=int, default=1,
                    help='Number of threads to use (default: %(default)d).')

parser.add_argument('--intermediate', metavar='PATH', type=str, default=None,
                    help='Output folder for intermediate files (wont be ' +
                         'kept unless this option is set.')

parser.add_argument('--chunk_size', type=int, default=5000,
                    help='Number of query seqs to read in at once for epa-ng ' +
                         '(default: %(default)d).')

parser.add_argument('--print_cmds', default=False, action='store_true',
                    help='If specified, print out wrapped commands to screen')

parser.add_argument('-v', '--version', default=False, action='version',
                    version="%(prog)s " + __version__)

def main():

    args = parser.parse_args()

    # Check that input filenames exist.
    check_files_exist([args.study_fasta, args.ref_msa, args.tree])	

    # If intermediate output directory set then create and output there.
    # Otherwise make a temporary directory for the intermediate files.
    if args.intermediate:

        make_output_dir(args.intermediate)

        place_seqs_pipeline(study_fasta=args.study_fasta,
                            ref_msa=args.ref_msa,
                            tree=args.tree,
                            hmm=args.hmm,
                            out_tree=args.out_tree,
                            alignment_tool=args.alignment_tool,
                            threads=args.threads,
                            out_dir=args.intermediate,
                            chunk_size=args.chunk_size,
                            print_cmds=args.print_cmds)

    else:
        with TemporaryDirectory() as temp_dir:
                place_seqs_pipeline(study_fasta=args.study_fasta,
                                    ref_msa=args.ref_msa,
                                    tree=args.tree,
                                    hmm=args.hmm,
                                    out_tree=args.out_tree,
                                    alignment_tool=args.alignment_tool,
                                    threads=args.threads,
                                    out_dir=temp_dir,
                                    chunk_size=args.chunk_size,
                                    print_cmds=args.print_cmds)


if __name__ == "__main__":
    main()
