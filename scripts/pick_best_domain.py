#!/usr/bin/env python

import argparse
from importlib.metadata import version
from picrust2.split_domains import get_lowest_nsti
from picrust2.util import get_tree_nodes, prune_tree, check_files_exist, read_fasta_ids
from picrust2.default import default_ref_dir_bac, default_ref_dir_arc
import sys
import pandas as pd
from os import path

parser = argparse.ArgumentParser(

    description="This script chooses the best domain to use for predictions for each "
                "input sequence. Typically this script would be called after predicting "
                "the copy number of gene families and the NSTI values, but before predicting "
                "any other traits. Given two NSTI tables and trees (typically from placing sequences in "
                "both the bacterial and archaeal trees), this script will output a table that "
                "contains the best domain for each sequence, and the input trees will be filtered to "
                "reflect this.",
    epilog='''
Usage example:
pick_best_domain.py \ \n
              -n1 bacteria \ \n
              -n2 archaea \ \n
              --tree_dom1 bac_out.tre \ \n
              --tree_dom2 arc_out.tre \ \n
              --tree_out_dom1 bac_reduced_out.tre \ \n
              --tree_out_dom2 arc_reduced_out.tre \ \n
              --nsti_table_dom1 bac_marker_nsti_predicted.tsv.gz \ \n
              --nsti_table_dom2 arc_marker_nsti_predicted.tsv.gz \ \n
              --nsti_table_out_dom1 bac_best_marker_nsti_predicted.tsv.gz \ \n
              --nsti_table_out_dom2 arc_best_marker_nsti_predicted.tsv.gz \ \n
              --nsti_table_out_combined combined_marker_nsti_predicted.tsv.gz \

''', formatter_class=argparse.RawDescriptionHelpFormatter)
                         
parser.add_argument('--tree_dom1', metavar='PATH', type=str, default=None,
                    help='The full tree for the first domain in newick format containing '
                         'both study sequences (i.e. ASVs or OTUs) and reference sequences.')

parser.add_argument('--tree_dom2', metavar='PATH', type=str, default=None,
                    help='The full tree for the second domain in newick format containing '
                         'both study sequences (i.e. ASVs or OTUs) and reference sequences.')

parser.add_argument('--tree_out_dom1', metavar='PATH', required=True, type=str,
                    help='Output reference tree for the first domain in newick format containing '
                         'both study sequences (i.e. ASVs or OTUs) and reference sequences. Only '
                         'study sequences that match the first domain will remain in this tree.')

parser.add_argument('--tree_out_dom2', metavar='PATH', required=True, type=str,
                    help='Output reference tree for the second domain in newick format containing '
                         'both study sequences (i.e. ASVs or OTUs) and reference sequences. Only '
                         'study sequences that match the second domain will remain in this tree.')
                         
parser.add_argument('--nsti_table_dom1', metavar='PATH', required=True, type=str,
                    help='The input NSTI table for the first domain in tab-delimited format.')

parser.add_argument('--nsti_table_dom2', metavar='PATH', required=True, type=str,
                    help='The input NSTI table for the second domain in tab-delimited format.')
                    
parser.add_argument('--nsti_table_out_dom1', metavar='PATH', required=True, type=str,
                    help='Output NSTI table for the first domain in tab-delimited format. '
                    'Only study sequences that matched the first domain will remain in this table. '
                    'If the extension \".gz\" is added the table will automatically be gzipped.')

parser.add_argument('--nsti_table_out_dom2', metavar='PATH', required=True, type=str,
                    help='Output NSTI table for the second domain in tab-delimited format. '
                    'Only study sequences that matched the second domain will remain in this table. '
                    'If the extension \".gz\" is added the table will automatically be gzipped.')
                    
parser.add_argument('--nsti_table_out_combined', metavar='PATH', required=True, type=str,
                    help='Output table with predicted NSTI, marker gene copy numbers and best '
                         'domain match per study sequence in input tree. If the extension \".gz\" '
                         'is added the table will automatically be gzipped.')
                         
parser.add_argument('--ref_fasta_dom1', metavar='PATH', type=str, default=None,
                    help='The reference FASTA file for the first domain.')
                    
parser.add_argument('--ref_fasta_dom2', metavar='PATH', type=str, default=None,
                    help='The reference FASTA file for the second domain.')

parser.add_argument('-n1', '--name_dom1', type=str, default="bacteria",
                    help='Name of the domain associated with the first tree and NSTI table. '
                    '(default: %(default)s).')

parser.add_argument('-n2', '--name_dom2', type=str, default="archaea",
                    help='Name of the domain associated with the second tree and NSTI table. '
                    '(default: %(default)s).')
                    
parser.add_argument('--check', default=True, action='store_true',
                    help='Check that the input trees match the input NSTI files.')
                    
parser.add_argument('--verbose', default=False, action='store_true',
                    help='If specified, print out wrapped commands and other '
                         'details to screen.')


def main():

    args = parser.parse_args()
    
    if args.name_dom1 in ['bac', 'bacteria', 'arc', 'archaea'] and args.name_dom2 in ['bac', 'bacteria', 'arc', 'archaea']:
        if args.name_dom1 == args.name_dom2:
            sys.exit("Stopping - gave the same name for --name_dom1 and --name_dom2.")
        if args.ref_fasta_dom1 is None and args.ref_fasta_dom2 is None:
            if args.name_dom1 in ['bac', 'bacteria']:
                ref_dir_dom1 = default_ref_dir_bac
            elif args.name_dom1 in ['arc', 'archaea']:
                ref_dir_dom1 = default_ref_dir_bac
            if args.name_dom2 in ['bac', 'bacteria']:
                ref_dir_dom2 = default_ref_dir_bac
            elif args.name_dom2 in ['arc', 'archaea']:
                ref_dir_dom2 = default_ref_dir_arc
                
            in_dir_dom1, in_dir_dom2 = ref_dir_dom1.rstrip('/'), ref_dir_dom2.rstrip('/')
            base_path_dom1, base_path_dom2 = path.basename(in_dir_dom1), path.basename(in_dir_dom2)
                
            args.ref_fasta_dom1, args.ref_fasta_dom2 = path.join(in_dir_dom1, base_path_dom1 + ".fna"), path.join(in_dir_dom2, base_path_dom2 + ".fna")
        
        else:
            sys.exit("You gave a reference fasta file for only one of the domains. Please give one for both or neither.")
    
    # Check that input filenames exist.
    check_files_exist([args.tree_dom1, args.tree_dom2, args.nsti_table_dom1, args.nsti_table_dom2])
    
    # If check == True, check whether the trees that have been given match the sequences in the NSTI files.
    # If they don't match, exit with an error message.
    if args.check:
        nodes_tree1 = get_tree_nodes(args.tree_dom1)
        nodes_tree2 = get_tree_nodes(args.tree_dom2)
        seqs_dom1 = list(pd.read_csv(filepath_or_buffer=args.nsti_table_dom1, index_col=0, header=0, sep='\t').index.values)
        seqs_dom2 = list(pd.read_csv(filepath_or_buffer=args.nsti_table_dom2, index_col=0, header=0, sep='\t').index.values)
        seqs_in_tree_dom1 = False
        seqs_in_tree_dom2 = False
        for seq in seqs_dom1:
            if seq in nodes_tree1:
                seqs_in_tree_dom1 = True
        
        for seq in seqs_dom2:
            if seq in nodes_tree2:
                seqs_in_tree_dom2 = True
        
        if not seqs_in_tree_dom1 or not seqs_in_tree_dom2:
            sys.exit("Stopping - one or both of the trees and NSTI files don't match.\n"
                     "For domain 1 (" + args.name_dom1 + "), you gave the tree file as: " + str(args.tree_dom1) + "\n"
                     "and the NSTI file as: " + str(args.nsti_table_dom1) + "\n"
                     "For domain 2 (" + args.name_dom2 + "), you gave the tree file as: " + str(args.tree_dom2) + "\n"
                     "and the NSTI file as: " + str(args.nsti_table_dom2) + "\n")
   
    # Get the best matching domain (lowest NSTI) for each of the study sequences
    nsti_lowest_df, nsti_dom1_df, nsti_dom2_df = get_lowest_nsti(args.nsti_table_dom1, 
                                                                args.nsti_table_dom2, 
                                                                args.name_dom1, 
                                                                args.name_dom2,
                                                                verbose=args.verbose)
    
    nsti_dom1_df.to_csv(path_or_buf=args.nsti_table_out_dom1, sep="\t")
    nsti_dom2_df.to_csv(path_or_buf=args.nsti_table_out_dom2, sep="\t")
    nsti_lowest_df.to_csv(path_or_buf=args.nsti_table_out_combined, sep="\t")
    
    # Now prune the trees to contain only the sequences of interest for each domain
    ref_d1_seqs = read_fasta_ids(args.ref_fasta_dom1)
    ref_d2_seqs = read_fasta_ids(args.ref_fasta_dom2)
    prune_tree(list(nsti_dom1_df.index.values)+ref_d1_seqs, args.tree_dom1, args.tree_out_dom1)
    prune_tree(list(nsti_dom2_df.index.values)+ref_d2_seqs, args.tree_dom2, args.tree_out_dom2)
    
    


if __name__ == "__main__":
    main()
