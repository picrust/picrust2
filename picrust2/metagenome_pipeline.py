#!/usr/bin/env python

__copyright__ = "Copyright 2018-2019, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.1.4-b"

import sys
import biom
import pandas as pd
import numpy as np
from os import path
from picrust2.util import (read_seqabun, make_output_dir,
                           three_df_index_overlap_sort)


def run_metagenome_pipeline(input_seqabun,
                            function,
                            marker,
                            max_nsti,
                            min_reads=1,
                            min_samples=1,
                            metagenome_contrib=False,
                            strat_out=False,
                            out_dir='metagenome_out'):
    '''Main function to run full metagenome pipeline. Meant to run modular
    functions largely listed below. Will return predicted metagenomes
    straitifed and unstratified by contributing genomes (i.e. taxa).'''

    # Read in input table of sequence abundances.

    study_seq_counts = read_seqabun(input_seqabun)

    # Read in predicted function and marker gene abundances.
    pred_function = pd.read_csv(function, sep="\t", index_col="sequence")
    pred_marker = pd.read_csv(marker, sep="\t", index_col="sequence")

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

    # Determine which sequences should be in the "RARE" category if either
    # the metagenome contributions table or stratified table is requested.
    if metagenome_contrib or strat_out:
        rare_seqs = []

        if min_reads != 1 or min_samples != 1:
            rare_seqs = id_rare_seqs(in_counts=study_seq_counts,
                                     min_reads=min_reads,
                                     min_samples=min_samples)

    # Output metagenome contributions table if specified.
    if metagenome_contrib:
        metagenome_contib_output = metagenome_contrib(func_abun=pred_function,
                                                      sample_abun=study_seq_counts,
                                                      rare_seqs=rare_seqs)

        metagenome_contib_outfile = path.join(out_dir,
                                              "metagenome_contrib.tsv.gz")

        metagenome_contib_output.to_csv(path_or_buf=metagenome_contib_outfile,
                                        sep="\t", index=False,
                                        compression="gzip")

    # Generate tables of functions by sample and return (either stratified or
    # not).
    if strat_out:

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
                           ignore_index=True, axis=0, sort=True)

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
        strat_func = pd.concat([strat_func, raw_seqs_slice], axis=0, sort=True)

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
    unstrat_func = pd.concat(sample_funcs, axis=1, sort=True)

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


def metagenome_contrib(func_abun, sample_abun, rare_seqs=[]):
    '''Take in function table and study sequence abundance table. Returns
    long-form table of how each sequence contributes functions in each
    sample. Note that the old format of columns (such as calling the sequences
    "OTUs" is retained here for backwards compatability. A subset of input
    sequences will be collapsed to a single category called "RARE" if a
    non-empty list is input for the rare_seqs option.'''

    # Convert sample abundance table to relative abundance since this is
    # expected by many tools that use the metagenome contributions table.
    sample_abun = sample_abun.div(sample_abun.sum(axis=0), axis=1) * 100

    # Counter used to identify the first sample.
    s_i = 0

    for sample in sample_abun.columns:
        single_abun = sample_abun[sample]
        single_abun = single_abun.iloc[single_abun.to_numpy().nonzero()]

        func_abun_subset = func_abun.loc[single_abun.index, :]

        # Melt function table sto be long format.
        func_abun_subset['OTU'] = func_abun_subset.index.to_list()

        func_abun_subset_melt = pd.melt(func_abun_subset, id_vars='OTU',
                                        value_name='GeneCountPerGenome',
                                        var_name='Gene')

        # Remove rows where gene count is 0.
        func_abun_subset_melt = func_abun_subset_melt[func_abun_subset_melt.GeneCountPerGenome != 0]

        func_abun_subset_melt['OTUAbundanceInSample'] = single_abun.loc[func_abun_subset_melt['OTU'].to_list()].to_list()

        func_abun_subset_melt['CountContributedByOTU'] = func_abun_subset_melt['GeneCountPerGenome'] * func_abun_subset_melt['OTUAbundanceInSample']

        # Collapse sequences identified as "rare" to the same category.
        rare_seqs = [r for r in rare_seqs if r in func_abun_subset_melt['OTU']]

        if len(rare_seqs) > 0:
            func_abun_subset_melt.loc[func_abun_subset_melt['OTU'].isin(rare_seqs), 'OTU'] = 'RARE' 
            func_abun_subset_melt = func_abun_subset_melt.groupby(['Gene', 'OTU'], as_index=False).sum()

        func_abun_subset_melt['Sample'] = sample

        # Order column names.
        func_abun_subset_melt = func_abun_subset_melt[['Sample',
                                                       'Gene',
                                                       'OTU',
                                                       'OTUAbundanceInSample',
                                                       'GeneCountPerGenome',
                                                       'CountContributedByOTU']]

        if s_i > 0:
            metagenome_contrib_out = pd.concat([metagenome_contrib_out,
                                                func_abun_subset_melt],
                                                axis=0)
        else:
            metagenome_contrib_out = func_abun_subset_melt.copy()
            s_i += 1

    return(metagenome_contrib_out)

