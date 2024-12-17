#!/usr/bin/env python

from os import path
import pandas as pd
import sys
from picrust2.util import (prune_tree, read_fasta_ids)

def get_lowest_nsti(nsti_table_dom1_path,
                    nsti_table_dom2_path,
                    dom1_name,
                    dom2_name,
                    verbose=False):
  
    '''Function that chooses which of the two domains used are most appropriate
    for the predictions for each sequence. 
    Descriptions of all of these input arguments/options are given in the 
    scripts/pick_best_domain.py script.'''
    
    if verbose:
        print("Picking the best domains for all sequences.", file=sys.stderr)
    
    # First read in the data frames containing the NSTI values for each domain.
    nsti_dom1 = pd.read_csv(filepath_or_buffer=nsti_table_dom1_path, sep="\t", header=0, index_col=0)
    nsti_dom2 = pd.read_csv(filepath_or_buffer=nsti_table_dom2_path, sep="\t", header=0, index_col=0)
    
    # Check that each file has an NSTI column
    if 'metadata_NSTI' not in nsti_dom1.columns or 'metadata_NSTI' not in nsti_dom2.columns:
        sys.exit("Stopping - one or both of the files don't contain an column named metadata_NSTI.\n")
    
    # Get the names of all sequences in both trees
    seqs = set(list(nsti_dom1.index.values)+list(nsti_dom2.index.values))
    
    # Make an empty dataframe for the lowest NSTI values, the marker gene copy number, and which of the domains was used for each sequence
    nsti_lowest = pd.DataFrame.from_dict({seq:['', '', '', ''] for seq in seqs}, orient='index', columns=['16S_rRNA_Count', 'metadata_NSTI', 'best_domain', 'closest_reference_genome'])
    nsti_lowest.index.name = 'sequence'
    
    # For each of the sequences, go through and choose the best domain.
    # Note that if the NSTI values are the same for both, dom1 is taken by default.
    for seq in seqs:
        if seq not in nsti_dom1.index.values:
          best, nsti_val, marker_copy_number, closest = dom2_name, nsti_dom2.loc[seq, 'metadata_NSTI'], nsti_dom2.loc[seq, '16S_rRNA_Count'], nsti_dom2.loc[seq, 'closest_reference_genome']
        elif seq not in nsti_dom2.index.values:
          best, nsti_val, marker_copy_number, closest = dom1_name, nsti_dom1.loc[seq, 'metadata_NSTI'], nsti_dom1.loc[seq, '16S_rRNA_Count'], nsti_dom1.loc[seq, 'closest_reference_genome']
        elif nsti_dom1.loc[seq, 'metadata_NSTI'] == nsti_dom2.loc[seq, 'metadata_NSTI']:
          print("Both domains had exactly the same NSTI value for this sequence: "
              + seq + ". Using the first domain: " + dom1_name, file=sys.stderr)
          best, nsti_val, marker_copy_number, closest = dom1_name, nsti_dom1.loc[seq, 'metadata_NSTI'], nsti_dom1.loc[seq, '16S_rRNA_Count'], nsti_dom1.loc[seq, 'closest_reference_genome']
        elif nsti_dom1.loc[seq, 'metadata_NSTI'] < nsti_dom2.loc[seq, 'metadata_NSTI']:
          best, nsti_val, marker_copy_number, closest = dom1_name, nsti_dom1.loc[seq, 'metadata_NSTI'], nsti_dom1.loc[seq, '16S_rRNA_Count'], nsti_dom1.loc[seq, 'closest_reference_genome']
        elif nsti_dom2.loc[seq, 'metadata_NSTI'] < nsti_dom1.loc[seq, 'metadata_NSTI']:
          best, nsti_val, marker_copy_number, closest = dom2_name, nsti_dom2.loc[seq, 'metadata_NSTI'], nsti_dom2.loc[seq, '16S_rRNA_Count'], nsti_dom2.loc[seq, 'closest_reference_genome']
        
        nsti_lowest.loc[seq, :] = marker_copy_number, nsti_val, best, closest
    
    # Filter the original NSTI tables to only contain the sequences that fit best with each domain.
    nsti_dom1 = nsti_dom1.loc[nsti_lowest[nsti_lowest['best_domain'] == dom1_name].index.values, :]
    nsti_dom2 = nsti_dom2.loc[nsti_lowest[nsti_lowest['best_domain'] == dom2_name].index.values, :]
    
    if verbose:
        print(str(nsti_dom1.shape[0]) + " sequences were kept for " + dom1_name +
              "\n"+str(nsti_dom2.shape[0]) + " sequences were kept for " + dom2_name, file=sys.stderr)
    
    return(nsti_lowest, nsti_dom1, nsti_dom2)
  
def combine_domain_predictions(in_path_dom1,
                              in_path_dom2,
                              out_path,
                              verbose=False):
    '''Function for adding together the functional predictions from two separate domains.
    Descriptions of all of these input arguments/options are given in the 
    scripts/combine_domains.py script.'''
    
    # First read in the data frames containing the predictions for each domain.
    in_dom1 = pd.read_csv(filepath_or_buffer=in_path_dom1, sep="\t", header=0, index_col=0)
    in_dom2 = pd.read_csv(filepath_or_buffer=in_path_dom2, sep="\t", header=0, index_col=0)
    
    # Check for duplicated sequence names
    overlap = [seq for seq in in_dom1.index.values if seq in in_dom2.index.values]
    if len(overlap) > 0:
        sys.exit("Stopping - Can't combine domains when there are sequences that are present in both input tables.\n"
                "Your two input files are: " + str(in_path_dom1) + " and " + str(in_path_dom2) + "\n"
                "The following sequences are present in both files: " + ",".join(overlap) + "\n"
                "Check that there are no sequences present in both files before continuing.")
    
    # Check to see whether any functions overlap between the tables.
    # It's not essential for there to be overlap, but if there is none then it probably indicates something weird going on.
    overlap_functions = [func for func in list(in_dom1.columns) if func in list(in_dom2.columns)]
    if len(overlap_functions) == 0:
        print("None of the functions overlapped between the two input files. It's not essential for there to be overlap,"
              "so this will continue running, but it is weird. You might want to check the input files. These are: "
              + str(in_path_dom1) + " and " + str(in_path_dom2), file=sys.stderr)
    
    # Combine the data frames, filling any functions that are not present in both with 0
    combined = pd.concat([in_dom1, in_dom2]).fillna(value=0)
    
    # Save the combined data frame
    combined.to_csv(path_or_buf=out_path, sep='\t')
    
    return
