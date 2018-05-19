#!/usr/bin/env python

from __future__ import division

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2-alpha.10"

import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
import numpy as np
import pandas as pd
import itertools


def read_in_rds(rds_path):
    '''Function to read in and return R object in RDS format.'''

    readRDS = robjects.r['readRDS']
    return(readRDS(rds_path))


def sample_func_prob_dist(func_probs, seq_counts, sample_id):
    '''Return probability distribution for specified function in sample.'''

    # Multiply per-sequence func abundances by sequence abundances.
    func_abun = (func_probs > 0).multiply(seq_counts[sample_id],
                                          axis='index')

    # Convert type from float to int.
    func_abun = func_abun.astype(int)

    # Make sure that both dataframes are ordered the same.
    func_probs = func_probs.reindex(func_abun.index)

    # Identify rows that are > 0 in columns besides the first.    
    non_zero_rows = list(func_abun.iloc[:, 1:].sum(axis=1) > 0)

    # If there are no rows that match this description then return 0.
    if sum(non_zero_rows) == 0:
        return(np.array([0]), np.array([1]))

    # Otherwise subset to only the non-zero rows.
    func_probs = func_probs.iloc[non_zero_rows, :]
    func_abun = func_abun.iloc[non_zero_rows, :]

    # Parse out rows that are single counts.
    nonambig_rows = list((func_abun > 0).sum(axis=1) == 1)

    ambig_rows = [not x for x in nonambig_rows]

    starting_count = 0

    # If there are some nonambig, then get how many counts based on these rows.
    if sum(nonambig_rows) > 0:

        func_abun_nonambig = func_abun.iloc[nonambig_rows, :]
        starting_count = (func_abun_nonambig * func_abun_nonambig.columns.values).values.sum()
        func_probs = func_probs.iloc[ambig_rows, :]
        func_abun = func_abun.iloc[ambig_rows, :]

    # If there are no ambiguous rows then return starting count calculated above.
    if sum(ambig_rows) == 0:
        return(np.array([starting_count]), np.array([1]))

    # Multiply func abundance by column values (and add starting counts from 
    # unambiguous rows)
    func_abun = (func_abun * func_abun.columns.values) + starting_count

    # Set first column to be 0.
    if 0 in func_abun.columns:
        func_abun[0] = 0

    return(convolve_mult_prob_dist(func_probs, func_abun))


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


def process_func_count_prob(predict_func_probs_raw, func, study_seq_counts):

    '''Determines expected function count and confidence intervals per sample
    and returns dictionary of these values'''

    func_probs = pandas2ri.ri2py_dataframe(predict_func_probs_raw)

    func_rownames = pandas2ri.ri2py(predict_func_probs_raw.rownames)

    # Set index and column names to be same as in original R matrix.
    func_probs.columns = [int(i) for i in predict_func_probs_raw.colnames] 
    func_probs.set_index(func_rownames, inplace=True)

    # Subset to sequences for which there is sequence abundance information.
    func_probs = func_probs.loc[study_seq_counts.index.values]

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


def convolve_mult_prob_dist(prob_input, abun_input):
    """
    Convolve probability distributions together using FFTs. Two input pandas
    dataframes are required. "prob_input" contains the probabilities of each
    abundance for each distribution (rows). "abun_input" contains the abundance 
    that each probability in prob_input refers to. The exception is that there 
    is no abundance of 0 values, which are the first column of the probability 
    dataframe. Returns numpy arrays of the new probabilities and abundances.
    """

    # Get maximum count possible based on abundances for each distribution.
    max_count = sum(np.amax(abun_input, axis=1))

    # Initialize empty array of size (# distributions) x (max count + 1)
    prob_array = np.zeros((abun_input.shape[0], max_count + 1))

    # Add in probabilities for abundance of 0 (if the first column corresponds)
    # to 0.
    if prob_input.columns.values[0] == 0:
        prob_array[: , 0] = prob_input.iloc[ : , 0]

    # Add in probabilities from each distribution for the corresponding
    # abundance values.
    i = 0
    for rowname, abun_row in abun_input.iterrows():

        nonzero_rows = list(abun_row > 0)

        prob_array[i, list(abun_row.loc[nonzero_rows])] = \
                                               prob_input.iloc[i, nonzero_rows]
        i += 1

    # Transform, take the product, and do the inverse transform
    # to get the convolution.
    fft_dists = np.fft.fft(prob_array)
    fft_convolution = fft_dists.prod(axis=0)
    prob_convolution = np.fft.ifft(fft_convolution)

    # Drop absolute values less than 1e-10
    prob_convolution[abs(prob_convolution) < 1e-10] = 0.0

    # Get real values only and identify non-zero indices.
    prob_convolution = prob_convolution.real

    # Raise error if summed probabilty is not equal to 1.
    prob_convolution_sum = sum(prob_convolution)
    if not np.isclose(1.0, prob_convolution_sum):
        raise RuntimeError("Stopping program, sum of probabilities after " +
                           "convolution was " + str(prob_convolution_sum) +
                           " instead of 1.0")

    nonzero_abundances = np.argwhere(prob_convolution > 0).flatten()

    # Get non-zero, real values of prob distribution as list.
    prob_convolution = prob_convolution[nonzero_abundances]

    return nonzero_abundances, prob_convolution
