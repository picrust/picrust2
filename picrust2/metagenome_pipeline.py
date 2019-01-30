#!/usr/bin/env python

from __future__ import division

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.0.4-b"

import sys
import biom
import pandas as pd
import numpy as np
from os import path
from picrust2.util import (biom_to_pandas_df, make_output_dir,
                           three_df_index_overlap_sort)


def run_metagenome_pipeline(input_biom,
                            function,
                            marker,
                            max_nsti,
                            min_reads=1,
                            min_samples=1,
                            strat_out=False,
                            out_dir='metagenome_out'):
    '''Main function to run full metagenome pipeline. Meant to run modular
    functions largely listed below. Will return predicted metagenomes
    straitifed and unstratified by contributing genomes (i.e. taxa).'''

    # Read in input table of sequence abundances and convert to pandas df.
    study_seq_counts = biom_to_pandas_df(biom.load_table(input_biom))

    # Read in predicted function and marker gene abundances.
    pred_function = pd.read_table(function, sep="\t", index_col="sequence")
    pred_marker = pd.read_table(marker, sep="\t", index_col="sequence")

    pred_function.index = pred_function.index.astype(str)
    pred_marker.index = pred_marker.index.astype(str)

    # Initialize empty pandas dataframe to contain NSTI values.
    nsti_val = pd.DataFrame()

    # If NSTI column present then remove all rows with value above specified
    # max value. Also, remove NSTI column (in both dataframes).
    if 'metadata_NSTI' in pred_function.columns:
        pred_function, nsti_val = drop_tips_by_nsti(tab=pred_function,
                                                    nsti_col='metadata_NSTI',
                                                    max_nsti=max_nsti)

    if 'metadata_NSTI' in pred_marker.columns:
        pred_marker, nsti_val = drop_tips_by_nsti(tab=pred_marker,
                                                  nsti_col='metadata_NSTI',
                                                  max_nsti=max_nsti)

    # Re-order predicted abundance tables to be in same order as study seqs.
    # Also, drop any sequence ids that don't overlap across all dataframes.
    study_seq_counts, pred_function, pred_marker = three_df_index_overlap_sort(study_seq_counts,
                                                                               pred_function,
                                                                               pred_marker)

    # Create output directory if it does not already exist.
    make_output_dir(out_dir)

    # Create normalized sequence abundance filename.
    norm_output = path.join(out_dir, "seqtab_norm.tsv")

    # Normalize input study sequence abundances by predicted abundance of
    # marker genes and output normalized table if specified.
    study_seq_counts = norm_by_marker_copies(input_seq_counts=study_seq_counts,
                                             input_marker_num=pred_marker,
                                             norm_filename=norm_output)

    # If NSTI column input then output weighted NSTI values.
    if not nsti_val.empty:
        weighted_nsti_out = path.join(out_dir, "weighted_nsti.tsv")
        calc_weighted_nsti(seq_counts=study_seq_counts,
                           nsti_input=nsti_val,
                           outfile=weighted_nsti_out)

    # Generate tables of functions by sample and return (either stratified or
    # not).
    if strat_out:
        rare_seqs = []

        # Determine which sequences should be in the "RARE" category.
        if min_reads != 1 or min_samples != 1:
            rare_seqs = id_rare_seqs(in_counts=study_seq_counts,
                                     min_reads=min_reads,
                                     min_samples=min_samples)

        return(strat_funcs_by_samples(pred_function, study_seq_counts,
                                      rare_seqs))
    else:
        return(None, unstrat_funcs_only_by_samples(pred_function,
                                                   study_seq_counts))


def strat_funcs_by_samples(func_abun, sample_abun, rare_seqs=[],
                           return_unstrat=True):
    '''Take in function table and study sequence abundance table. Returns
    stratified table of function abundances per sequence (rows) per samples
    (columns). Will also return unstratified format (function by sample).
    Will collapse rare sequences into a single sequence category if
    given a non-empty list of ids.'''

    # Create combined stratified dataframe.
    strat_func = pd.concat((sample_abun.multiply(func_abun[func], axis=0) for func in func_abun.columns),
                           ignore_index=True, axis=0)

    # Set multi-index labels to be sequence and function ids.
    strat_func.index = pd.MultiIndex.from_product((func_abun.columns,
                                                   func_abun.index))
    strat_func.index.names = ['function', 'sequence']

    if len(rare_seqs) > 0:
        # Sum rare rows and make a new row named "RARE" in each df.
        # Remove all the original rare rows.

        # First identify rows to remove, slice them out, and remove from original.
        rare_seqs_rows = strat_func.index.get_level_values('sequence').isin(rare_seqs)
        raw_seqs_slice = strat_func.iloc[rare_seqs_rows]
        strat_func = strat_func.iloc[~rare_seqs_rows]

        # Next remove sequence index and sum based on function id. Re-add "RARE"
        # to be the sequence id.
        raw_seqs_slice.index = raw_seqs_slice.index.droplevel('sequence')
        raw_seqs_slice = raw_seqs_slice.sum(level=['function'])
        raw_seqs_slice['sequence'] = 'RARE'
        raw_seqs_slice.set_index('sequence', append=True, inplace=True)

        # Concat the RARE seqs to the full df.
        strat_func = pd.concat([strat_func, raw_seqs_slice], axis=0)

    # Remove rows that are all 0.
    strat_func = strat_func.loc[~(strat_func == 0).all(axis=1)]

    # Return dataframe and also unstratified dataframe if specified.
    if return_unstrat:
        return(strat_func, strat_func.sum(level='function', axis=0))
    else:
        return(strat_func)

def unstrat_funcs_only_by_samples(func_abun, sample_abun):
    '''Take in function table and study sequence abundance table. Returns
    unstratified table of function abundances by samples.'''

    # Sample ids are taken from sequence abundance table.
    sample_ids = sample_abun.columns.values

    # List that will contain a series for each sample's unstratified
    # function abundnaces.
    sample_funcs = []

    for sample in sample_ids:
        sample_funcs.append(func_abun.mul(sample_abun[sample], axis=0).sum(axis=0))

    # Build dataframe from these series.
    unstrat_func = pd.concat(sample_funcs, axis=1)

    unstrat_func = unstrat_func.loc[~(unstrat_func == 0).all(axis=1)]

    unstrat_func.columns = sample_ids

    unstrat_func.index.name = 'function'

    # Remove rows that are all 0s and return.
    return(unstrat_func)


def drop_tips_by_nsti(tab, nsti_col, max_nsti):
    '''Takes in (1) function table (columns are functions, rows are predicted
    genomes for ASVs, (2) which column of this function table corresponds to
    NSTI values, and (3) the maximum NSTI value for ASVs to be retained. Will
    return table with ASVs above the max NSTI cut-off excluded and a dataframe
    which just contains the NSTI values for ASVs that passed the cut-off. Will
    report how many ASVs were removed and will throw error if all ASVs are.'''

    orig_num_rows = tab.shape[0]

    tab = tab[tab[nsti_col] <= max_nsti]

    filt_num_rows = tab.shape[0]

    if filt_num_rows == 0:
        sys.exit("Stopping - all ASVs filtered from table when max NSTI "
                 "cut-off of " + str(max_nsti) + " used.")
    else:
        num_removed = orig_num_rows - filt_num_rows
        print(str(num_removed) + " of " + str(orig_num_rows) + " ASVs were "
              "above the max NSTI cut-off of " + str(max_nsti) + " and were "
              "removed.", file=sys.stderr)

    # Keep track of NSTI column as separate dataframe and remove this column
    # from the main dataframe.
    nsti_val = tab[[nsti_col]]

    return(tab.drop(nsti_col, axis=1, inplace=False), nsti_val)


def calc_weighted_nsti(seq_counts, nsti_input, outfile=None, return_df=False):
    '''Will calculate weighted NSTI values given sequence count table and NSTI
    value for each sequence. Will output these weighted values to a file if
    output file is specified. Will only return a df if specified.'''

    nsti_mult = seq_counts.mul(nsti_input.metadata_NSTI, axis=0)

    # Get column sums divided by total abundance per sample.
    weighted_nsti = pd.DataFrame(nsti_mult.sum(axis=0)/seq_counts.sum(axis=0))

    weighted_nsti.columns = ["weighted_NSTI"]

    # Write to outfile if specified.
    if outfile:
        weighted_nsti.to_csv(path_or_buf=outfile, sep="\t", header=True,
                             index_label="sample")
    if return_df:
        return(weighted_nsti)


def norm_by_marker_copies(input_seq_counts,
                          input_marker_num,
                          norm_filename=None,
                          round_decimal=2):

    '''Divides sequence counts (which correspond to amplicon sequence
    variants) by the predicted marker gene copies for each sequence. Will write
    out the normalized table if option specified.'''

    input_seq_counts = input_seq_counts.div(input_marker_num.loc[
                                                input_seq_counts.index.values,
                                                input_marker_num.columns.values[0]],
                                            axis="index")

    input_seq_counts = input_seq_counts.round(decimals=round_decimal)

    # Output normalized table if specified.
    if norm_filename:
        input_seq_counts.to_csv(path_or_buf=norm_filename,
                                index_label="normalized",
                                sep="\t")

    return(input_seq_counts)


def id_rare_seqs(in_counts, min_reads, min_samples):
    '''Determine which rows of a sequence countfile are below either the
    cut-offs of min read counts or min samples present.'''

    # Check if "RARE" is the name of a sequence in this table.
    if "RARE" in in_counts.index:
        sys.exit("Stopping: the sequence called \"RARE\" in the sequence " +
                 "abundance table should be re-named.")

    low_freq_seq = set(in_counts[in_counts.sum(axis=1) < min_reads].index)
    few_samples_seq = set(in_counts[(in_counts != 0).astype(int).sum(axis=1) < min_samples].index)

    return(list(low_freq_seq.union(few_samples_seq)))


# ### DEPRECATED ###
# def funcs_by_sample(input_seq_counts, input_function_num, rare_seqs=[],
#                     strat_out=False, proc=1):
#     '''Function that reads in study sequence abundances and predicted
#     number of gene families per study sequence's predicted genome. Will
#     return abundance of functions in each sample (unstratified format). If
#     specified then will also return abundances stratified by contributing
#     sequences.'''

#     # Sample ids are taken from sequence abundance table.
#     sample_ids = input_seq_counts.columns.values

#     # Loop through all samples and get predicted functional abundances
#     # after multiplying each contributing sequence by the abundance in
#     # the sequence abundance dataframe.
#     if proc > 1:

#         sample_funcs = Parallel(n_jobs=proc)(delayed(
#                             func_by_seq_abun)(
#                             input_seq_counts[sample],
#                             input_function_num,
#                             rare_seqs,
#                             strat_out)
#                             for sample in sample_ids)
#     else:
#         # Run in basic loop if only 1 processor specified.
#         sample_funcs = []

#         for sample in sample_ids:
#             sample_funcs += [func_by_seq_abun(input_seq_counts[sample],
#                                               input_function_num, rare_seqs,
#                                               strat_out)]

#     # Build dataframe from list of series or dataframes (one per sample).
#     sample_funcs_df = pd.concat(sample_funcs, axis=1)

#     # Remove rows that are all 0s.
#     sample_funcs_df = sample_funcs_df.loc[~(sample_funcs_df==0).all(axis=1)]

#     # Set column names to be sample ids.
#     sample_funcs_df.columns = sample_ids

#     if strat_out:

#         # Sum rows by function id for unstratified output and set index labels
#         # equal to function ids (remove "sequence" column first) if unstratified
#         # returned above.
#         unstrat_out_df = sample_funcs_df.copy()

#         unstrat_out_df = pd.pivot_table(unstrat_out_df.reset_index().drop("sequence", axis=1),
#                                         index="function", aggfunc=np.sum)

#         return(sample_funcs_df, unstrat_out_df)

#     else:

#         return(None, sample_funcs_df)

# ### DEPRECATED ###
# def func_by_seq_abun(sample_seq_counts, func_abun, rare_seqs=[],
#                      calc_strat=False):
#     '''Given the abundances of sequences in a sample (as a pandas series) and
#     the predicted functions of those sequences (as a pandas dataframe), this
#     function will return the functional abundances after multiplying the
#     abundances of functions contributed by a sequence by that sequence's
#     abundance summed over all sequences (unstratified). If calc_strat=True then
#     instead will return a long-form dataframe with 3 columns: function,
#     sequence, and count. Note that this function assumes that the order of the
#     sequences is the same in both the input series and the dataframe of
#     function abundances.'''

#     func_abun_depth = func_abun.mul(sample_seq_counts, axis=0)

#     # Identify rows corresponding to rare seqs, remove them and add them back
#     # in as sum of all those rows with new name "RARE".
#     if rare_seqs:
#         rare_subset_sum = func_abun_depth.loc[rare_seqs].sum(axis=0)
#         func_abun_depth.drop(labels=rare_seqs, axis=0, inplace=True)
#         func_abun_depth.loc["RARE"] = rare_subset_sum

#     # Generate stratitifed table if option set, otherwise get sum per function.
#     if calc_strat:
#         # Set index labels (sequence ids) to be new column.
#         func_abun_depth["sequence"] = func_abun_depth.index.values

#         # Convert from wide to long table format (only columns for sequence,
#         # function, and count).
#         func_abun_depth_long = pd.melt(func_abun_depth, id_vars=["sequence"],
#                                        var_name="function", value_name="count")

#         func_abun_depth_long.transpose()

#         # Set index labels to be sequence and function columns.
#         func_abun_depth_long = func_abun_depth_long.set_index(keys=["function",
#                                                                     "sequence"])

#         # Convert long-form pandas dataframe to series and return.
#         return(func_abun_depth_long["count"])

#     else:
#         return(func_abun_depth.sum(axis=0))
