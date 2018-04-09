#!/usr/bin/env python

from __future__ import division

__license__ = "GPL"
__version__ = "2-alpha.6"

import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
import numpy as np
import pandas as pd
import itertools


def read_in_rds(rds_path):
    '''Function to read in and return R object in RDS format.'''

    readRDS = robjects.r['readRDS']
    return(readRDS(rds_path))


def fix_dup_counts(poss_counts, count_probs):
    '''Dereplicate same count values and sum their probabilities.'''

    unique_poss_counts = set(poss_counts)
    for count in unique_poss_counts:
        count_i = np.where(poss_counts == count)[0]

        # If this count value is duplicated then sum probs for all
        # instances.
        if len(count_i) > 1:
            summed_prob = 0
            for match_i in sorted(count_i, reverse=True):
                del poss_counts[match_i]

                summed_prob += count_probs[match_i]
                del count_probs[match_i]

            # Add count back in once at end along with summed probability.
            poss_counts += [count]
            count_probs += [summed_prob]

    return poss_counts, count_probs


def all_possible_counts(probs_in, abun_in, poss_counts):
    '''Loop through all rows of input dataframes and return possible counts
    and the probability of each count'''

    count_probs = [1]

    # Loop through all sequences.
    for seqname in probs_in.index.values:

        # Identify non-zero indices.
        nonzero_seq_col = list(probs_in.loc[seqname, ] > 0)

        # Subset to non-zero indices for count and prob rows.
        seq_poss_counts = abun_in.loc[seqname, nonzero_seq_col] *\
                          abun_in.columns.values[nonzero_seq_col]

        seq_count_probs = probs_in.loc[seqname, nonzero_seq_col]

        # Add these seq counts to possible counts.
        new_poss_count = []
        for seq_poss_count in seq_poss_counts:
            new_poss_count += [seq_poss_count + x for x in poss_counts]

        poss_counts = new_poss_count

        # Do the same procedure, but for count probabilities.
        new_count_probs = []
        for count_prob in seq_count_probs:
            new_count_probs += [count_prob * x for x in count_probs]

        count_probs = new_count_probs

    # Check if any duplicates in possible counts and correct if so.
    if len(set(poss_counts)) < len(poss_counts):
        poss_counts, count_probs = fix_dup_counts(poss_counts, count_probs)

    # Convert to numpy arrays and return them sorted by counts.
    poss_counts = np.array(poss_counts)
    count_probs = np.array(count_probs)

    sorted_count_i = poss_counts.argsort()

    return poss_counts[sorted_count_i], count_probs[sorted_count_i]


def sample_func_prob_dist(func_probs, seq_counts, sample_id):
    '''Return probability distribution for specified function in sample.'''

    # Multiply per-sequence func abundances by sequence abundances.
    func_abun = (func_probs > 0).multiply(seq_counts[sample_id],
                                          axis='index')

    # Make sure that both dataframes are ordered the same.
    func_probs = func_probs.reindex(func_abun.index)

    # Identify rows that are > 0 in columns besides the first.    
    non_zero_rows = list(func_abun.iloc[:, 1:-1].sum(axis=1) > 0)

    # If there are no rows that match this description then return 0.
    if sum(non_zero_rows) == 0:
        return([0], [1])

    # Otherwise subset to only the non-zero rows.
    func_probs = func_probs.iloc[non_zero_rows, :]
    func_abun = func_abun.iloc[non_zero_rows, :]

    # Parse out rows that are single counts.
    nonambig_rows = list((func_abun > 0).sum(axis=1) == 1)

    ambig_rows = [not x for x in nonambig_rows]

    # If there are some nonambig, then get how many counts based on these rows.
    if sum(nonambig_rows) > 0:
        func_abun_nonambig = func_abun.iloc[nonambig_rows, :]

        possible_counts = [(func_abun_nonambig * func_abun_nonambig.columns.values).values.sum()]

    # If there are no ambig rows then set possible counts to 0.
    else:
        possible_counts = [0]

    # If there are no ambiguous rows then return starting count calculated above.
    if sum(ambig_rows) == 0:
        return(possible_counts, [1])

    # Otherwise subset to ambiguous rows and get all possible counts.
    return(all_possible_counts(probs_in=func_probs.iloc[ambig_rows, :],
                               abun_in=func_abun.iloc[ambig_rows, :],
                               poss_counts=possible_counts))


def norm_by_marker_copies(input_seq_counts,
                          input_marker_num,
                          output_normfile=False,
                          norm_filename="norm_seq_counts.tsv"):

    '''Divides sequence counts (which correspond to amplicon sequence
    variants) by the predicted marker gene copies for each sequence. Will write
    out the normalized table if option specified.'''

    # Check that all rownames match between the two input dataframes.
    if sorted(input_seq_counts.index) != sorted(input_marker_num.index):
        raise ValueError("Sequence names do not match between input " +
                         "sequence abundances and marker gene copy number " +
                         "table. These names are: ",
                         list(input_seq_counts.index.values),
                         list(input_marker_num.index.values))

    input_seq_counts = input_seq_counts.div(input_marker_num.loc[:, "16S_rRNA_Count"],
                                            axis="index")

    # Output normalized table if specified.
    if output_normfile:
        input_seq_counts.to_csv(path_or_buf=norm_filename,
                                index_label="sequence",
                                sep="\t")

    return(input_seq_counts)


def expectation_and_ci_val(poss_counts, count_probs, rounded=True):
    '''Given lists of possible counts and probabilities of these counts, return
    list of expected count and lower and upper 95% CIs.'''

    # If only 1 possible count given then return it as all 3 values.
    if len(poss_counts) == 1:
        return([poss_counts[0], poss_counts[0], poss_counts[0]])

    # Get expected value.
    exp_count = sum(poss_counts*count_probs)

    # Get counts at probabilities of 0.05 and 0.95.
    count_probs_cumsum = np.cumsum(count_probs)
    lower_ci = poss_counts[np.amin(np.where(count_probs_cumsum >= 0.05))]
    upper_ci = poss_counts[np.amin(np.where(count_probs_cumsum >= 0.95))]

    if rounded:
        exp_count = np.around(exp_count, decimals=2)
        lower_ci = np.around(lower_ci, decimals=2)
        upper_ci = np.around(upper_ci, decimals=2)

    return([exp_count, lower_ci, upper_ci])


def process_func_count_prob(predict_func_probs_raw,
                            func,
                            study_seq_counts):

    func_probs = pandas2ri.ri2py_dataframe(predict_func_probs_raw)

    func_rownames = pandas2ri.ri2py(predict_func_probs_raw.rownames)

    func_probs.set_index(func_rownames, inplace=True)

    func_by_sample = {}

    for sample_id in study_seq_counts.columns.values:

        poss_counts, count_probs = sample_func_prob_dist(func_probs,
                                                         study_seq_counts,
                                                         sample_id)

        func_by_sample_out = expectation_and_ci_val(poss_counts, count_probs)

        func_by_sample[sample_id] = func_by_sample_out[0]
        func_by_sample[sample_id + "_ci5"] = func_by_sample_out[1]
        func_by_sample[sample_id + "_ci95"] = func_by_sample_out[2]

    # Add func name to dictionary so this can later be index labels in df.
    func_by_sample["function"] = func 

    return(func_by_sample)
