#!/usr/bin/env python

from __future__ import division

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2-alpha.8"

import biom
import itertools
import pandas as pd
import numpy as np
from os import path
from tempfile import TemporaryDirectory
from joblib import Parallel, delayed
from picrust2.util import (biom_to_pandas_df, make_output_dir,
                           three_df_index_overlap_sort)

def run_metagenome_pipeline(input_biom,
                            function,
                            marker,
                            out_dir='metagenome_out',
                            proc=1,
                            output_normfile=False):
    '''Main function to run full metagenome pipeline. Meant to run modular
    functions largely listed below. Will return predicted metagenomes
    straitifed and unstratified by contributing genomes (i.e. taxa).'''

    # Read in input table of sequence abundances and convert to pandas df.
    study_seq_counts = biom_to_pandas_df(biom.load_table(input_biom))

    # Read in predicted function and marker gene abundances.
    pred_function = pd.read_table(function, sep="\t", index_col="tips")
    pred_marker = pd.read_table(marker, sep="\t", index_col="tips")

    # Re-order predicted abundance tables to be in same order as study seqs.
    # Also, drop any sequence ids that don't overlap across all dataframes.
    study_seq_counts, pred_function, pred_marker = three_df_index_overlap_sort(study_seq_counts, 
                                                                               pred_function,
                                                                               pred_marker)

    # Create output directory if it does not already exist.
    make_output_dir(out_dir)

    # Create normalized sequence abundance filename if outfile specified.
    if output_normfile:
        norm_output = path.join(out_dir, "seqtab_norm.tsv")
    else:
        norm_output = None

    # Normalize input study sequence abundances by predicted abundance of
    # marker genes and output normalized table if specified.
    study_seq_counts = norm_by_marker_copies(input_seq_counts=study_seq_counts,
                                             input_marker_num=pred_marker,
                                             norm_filename=norm_output)

    # Get predicted function counts by sample, stratified by contributing
    # genomes and also separately unstratified.
    return(funcs_by_sample(input_seq_counts=study_seq_counts,
                           input_function_num=pred_function))

def norm_by_marker_copies(input_seq_counts,
                          input_marker_num,
                          norm_filename=None,
                          round_norm=True):

    '''Divides sequence counts (which correspond to amplicon sequence
    variants) by the predicted marker gene copies for each sequence. Will write
    out the normalized table if option specified.'''

    input_seq_counts = input_seq_counts.div(input_marker_num.loc[
                                                input_seq_counts.index.values,
                                                input_marker_num.columns.values[0]],
                                            axis="index")

    if round_norm:
        input_seq_counts = input_seq_counts.round(decimals=0)

    # Output normalized table if specified.
    if norm_filename:
        input_seq_counts.to_csv(path_or_buf=norm_filename,
                                index_label="sequence",
                                sep="\t")

    return(input_seq_counts)

def funcs_by_sample(input_seq_counts, input_function_num, proc=1):
    '''Function that reads in study sequence abundances and predicted
    number of gene families per study sequence's predicted genome. Will
    return abundance of functions contributed by each study sequence per
    sample in pandas dataframe and another dataframe in unstratified format.'''

    # Sample ids are taken from sequence abundance table.
    sample_ids = input_seq_counts.columns.values

    # Loop through all samples and get predicted functional abundances
    # after multiplying each contributing sequence by the abundance in
    # the sequence abundance dataframe.
    if proc > 1:
        strat_out = Parallel(n_jobs=proc)(delayed(
                            func_by_seq_abun)(
                            input_seq_counts[sample],
                            input_function_num)
                            for sample in sample_ids)
    else:
        # Run in basic loop if only 1 processor specified.
        strat_out = []

        for sample in sample_ids:
            strat_out += [func_by_seq_abun(input_seq_counts[sample],
                                           input_function_num)]

    # Build dataframe from list of dictionaries per sample.
    strat_out_df = pd.DataFrame(strat_out).transpose()

    # Set column names to be sample ids.
    strat_out_df.columns = sample_ids

    # Remove rows that are all 0s.
    strat_out_df = strat_out_df.loc[~(strat_out_df==0).all(axis=1)]

    # Add function names as column to dataframe.
    strat_out_df["function"] = list(map(lambda x: list(x)[0],
                                    strat_out_df.index.values))

    # Sum rows by function id for unstratified output and set index labels
    # equal to function ids.
    unstrat_out_df = strat_out_df.copy()

    unstrat_out_df = pd.pivot_table(unstrat_out_df,
                                    index="function",
                                    aggfunc=np.sum)

    # Also add sequence as column to stratified output.
    strat_out_df["sequence"] = list(map(lambda x: list(x)[1],
                                    strat_out_df.index.values))

    # Re-order columns.
    strat_out_df = strat_out_df[["function", "sequence"] + list(sample_ids)]

    return(strat_out_df, unstrat_out_df)

def func_by_seq_abun(sample_seq_counts, func_abun):
    '''Given the abundances of sequences in a sample (as a pandas series) and 
    the predicted functions of those sequences (as a pandas dataframe), this 
    function will return the functional abundances after multiplying the 
    abundances of functions contributed by a sequence by that sequence's
    abundance. Will return a long-form dataframe with 3 columns: function,
    sequence, and count. Note that this function assumes that the order of the
    sequences is the same in both the input series and the dataframe of
    function abundances.'''

    func_abun_depth = func_abun.mul(sample_seq_counts, axis=0)

    # Set index labels (sequence ids) to be new column.
    func_abun_depth["sequence"] = func_abun_depth.index.values

    # Convert from wide to long table format (only columns for
    # sequence, function, and count)
    func_abun_depth_long = pd.melt(func_abun_depth, id_vars=["sequence"],
                                   var_name="function", value_name="count")

    # Convert long-form pandas dataframe to dictionary.
    func_dict = {}

    # Loop through all rows of pandas dataframe and add values as nested
    # dictionaries.
    for index, row in func_abun_depth_long.iterrows():
        func_dict[tuple([row["function"], row["sequence"]])] = row["count"]

    return(func_dict)
