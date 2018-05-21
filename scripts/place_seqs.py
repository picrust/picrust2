#!/usr/bin/env python

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2-alpha.10"

import argparse
from tempfile import TemporaryDirectory
from picrust2.place_seqs import place_seqs_pipeline
from picrust2.precalc import default_fasta, default_tree
from picrust2.util import make_output_dir, check_files_exist

parser = argparse.ArgumentParser(

    description="Wrapper to prep tree before HSP steps. Requires unaligned " +
                "FASTA of study sequences. Users can specify a non-default " +
                "reference fasta and treefile if needed.",

    formatter_class=argparse.RawDescriptionHelpFormatter)


parser.add_argument('-s', '--study_fasta', metavar='PATH', required=True,
                    type=str, help='FASTA of unaligned study sequences')

parser.add_argument('-r', '--ref_msa', metavar='PATH', type=str,
                    default=default_fasta,
                    help='FASTA of aligned reference sequences (default: %(default)s).')

parser.add_argument('-t', '--tree', metavar='PATH', type=str,
                    default=default_tree,
                    help='Input tree based on aligned reference sequences. ' +
                         '(default: %(default)s).')

parser.add_argument('-o', '--out_tree', metavar='PATH', required=True,
                    type=str, help='Name of final output tree')

parser.add_argument('--threads', type=int, default=1,
                    help='Number of threads to use (default: %(default)d).')

parser.add_argument('--papara_output', metavar='PATH', type=str, default=None,
                    help='Path to PaPaRa output in Phylip format (will skip ' +
                         'PaPaRa step)')

parser.add_argument('--out_dir', metavar='PATH', type=str, default=None,
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
    if args.out_dir:

        make_output_dir(args.out_dir)

        place_seqs_pipeline(study_fasta=args.study_fasta,
                            ref_msa=args.ref_msa,
                            tree=args.tree,
                            out_tree=args.out_tree,
                            threads=args.threads,
                            papara_output=args.papara_output,
                            out_dir=args.out_dir,
                            chunk_size=args.chunk_size,
                            print_cmds=args.print_cmds)

    else:
        with TemporaryDirectory() as temp_dir:
                place_seqs_pipeline(study_fasta=args.study_fasta,
                                    ref_msa=args.ref_msa,
                                    tree=args.tree,
                                    out_tree=args.out_tree,
                                    threads=args.threads,
                                    papara_output=args.papara_output,
                                    out_dir=temp_dir,
                                    chunk_size=args.chunk_size,
                                    print_cmds=args.print_cmds)


if __name__ == "__main__":
    main()
