#!/usr/bin/env python

from __future__ import division

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.0.0-b.2"

from collections import defaultdict
from joblib import Parallel, delayed
from os import path
import pandas as pd
import numpy as np
from picrust2.util import (system_call_check, get_picrust_project_dir)


def run_minpath_pipeline(inputfile,
                         mapfile,
                         proc=1,
                         out_dir=None,
                         print_cmds=False):
    '''Pipeline containing full pipeline for reading input files, making
    calls to functions to run MinPath and to return an output table of
    predicted pathway abundances that can be written to a file.'''

    # Read in table of gene family abundances stratified by contributing
    # seqeunces.
    strat_in = read_strat_genes(inputfile)

    # Get list of sample ids.
    samples = [col for col in strat_in.columns
               if col not in ["function", "sequence"]]

    # Run minpath wrapper on all samples.
    # Note that input stratified table is subsetted to required columns only.
    sample_path_abun_raw = Parallel(n_jobs=proc)(delayed(
                                    minpath_wrapper)(sample_id,
                                    strat_in[["function", "sequence", sample_id]],
                                    mapfile, out_dir, print_cmds)
                                    for sample_id in samples)

    # Split the output into unstratified and stratified.
    sample_path_abun_raw_unstrat = []
    sample_path_abun_raw_strat = []

    for sample_output in sample_path_abun_raw:
        sample_path_abun_raw_unstrat += [sample_output[0]]
        sample_path_abun_raw_strat += [sample_output[1]]

    # Convert these returned lists of series into pandas dataframes.
    sample_path_abun_unstrat = pd.DataFrame(sample_path_abun_raw_unstrat)
    sample_path_abun_strat = pd.DataFrame(sample_path_abun_raw_strat)

    # Set index labels of unstratified dataframe to be sample names.
    sample_path_abun_unstrat.index = samples

    # Replace all missing values (NaN) with 0s (i.e. pathway was missing in
    # that sample) and transpose.
    sample_path_abun_unstrat = sample_path_abun_unstrat.fillna(0).transpose()
    sample_path_abun_strat = sample_path_abun_strat.fillna(0).transpose()

    # Round each value to 2 decimal places.
    sample_path_abun_unstrat = sample_path_abun_unstrat.round(decimals=2)
    sample_path_abun_strat = sample_path_abun_strat.round(decimals=2)

    # Add pathway and sequence as columns of stratified table.
    sample_path_abun_strat.reset_index(inplace=True)

    return(sample_path_abun_unstrat, sample_path_abun_strat)


def read_strat_genes(filename):
    '''Reads in gene abundancy table stratified by contributing sequences
    (output of metagenome_pipeline.py). If an unstratified file is input
    it will return an error.'''

    # Read in input file as pandas dataframe.
    input_df = pd.read_table(filename, sep="\t")

    # Check that expected columns are in table.
    if "function" not in input_df.columns or "sequence" not in input_df.columns:
        raise ValueError("Did not find at least one of the expected " +
                         "in input file (\"function\" and \"sequence\". " +
                         "Make sure the stratified metagenome predictions " +
                         "were input.")

    return(input_df)


def strat_to_unstrat_counts(strat_df, func_col="function"):
    '''Given a pandas dataframe with the columns "sequence", "function" (by
    default), and at least 1 sample column, will return the dataframe after
    removing sequence column and summing all functions per sample. Functions
    will be new index labels.'''

    # Drop column containing sequence ids.
    strat_df = strat_df.drop(["sequence"], axis=1)

    return(pd.pivot_table(data=strat_df, index=func_col, aggfunc=np.sum))


def identify_minpath_present(report_file):
    '''Parse MinPath report output file and returns set containing all pathways
    ids that were called as present.'''

    path_present = set()

    with open(report_file, "r") as minpath_report_in:
        for line in minpath_report_in:
            line_split = line.split()

            if int(line_split[7]) == 1:
                path_present.add(line_split[-1])

    return(path_present)


def parse_minpath_details(details_file: str, path_present: set):
    '''Parse MinPath details output file and returns dictionaries containing
    the abundances of gene families within each pathway and the ids of these
    gene families. Note that the pathways that were called as present in the
    MinPath report file need to be given as an input argument as a set.'''

    # Initialize dictionary that will contain gene family abundance per
    # pathway and one that contains gene family ids in each pathway.
    gf_abund = {}
    gf_names = {}

    # Boolean specifying that pathway in details file was called as
    # present by MinPath.
    present = False

    with open(details_file, "r") as minpath_details_in:
        for line in minpath_details_in:
            line_split = line.split()

            # If line starts with "path" then keep track of pathway name if
            # it was called as present in report file.
            if line_split[0] == "path":
                if line_split[-1] not in path_present:
                    present = False
                    continue

                present = True
                current_pathway = line_split[-1]

                # Initialize list containing gene family abundances and list
                # containing the matching gene family ids in the same order.
                gf_abund[current_pathway] = []
                gf_names[current_pathway] = []

                # Add in abundances of 0 for missing genes (and None as
                # placeholder for missing gene family ids).
                for i in range(int(line_split[3]) - int(line_split[5])):
                    gf_abund[current_pathway] += [0]
                    gf_names[current_pathway] += [None]

            # If line does not start with "path" then only proceed if current
            # pathway is present.
            elif present:
                gf_abund[current_pathway] += [int(float(line_split[2]))]
                gf_names[current_pathway] += [str(line_split[0])]

    return(gf_abund, gf_names)


def path_abun_by_seq(gene_abun, gene_ids, total_sum, path_abun):
    '''Takes in a stratified dataframe, subsets functions to those of interest,
    and pivots by sequence column (takes sum over other columns). Also takes
    in total sum of gene families that went into calculating pathway abundance
    and the calculated pathway abundance. Will return the weighted pathway
    abundance contributed by each sequence.'''

    # Subset to genes in pathway.
    gene_abun = gene_abun.loc[gene_abun['function'].isin(gene_ids)]

    # Drop function column.
    gene_abun = gene_abun.drop(["function"], axis=1)

    # Return dataframe with sum of all genes per sequence.
    seq_path_abun = pd.pivot_table(data=gene_abun, index="sequence",
                                   aggfunc=np.sum)

    # Return weighted pathway abundance (rounded).
    return(np.around((seq_path_abun/total_sum)*path_abun, decimals=2))


def minpath_wrapper(sample_id, strat_input, minpath_map, out_dir,
                    print_opt=False):
    '''Read in sample_id, gene family table, and out_dir, and run MinPath based
    on the gene family abundances. Returns both unstratified and stratified
    pathway abundances as dictionaries in a list.'''

    # Get gene family abundances summed over all sequences for this sample.
    unstrat_input = strat_to_unstrat_counts(strat_input)

    # Define MinPath input and outout filenames.
    minpath_in = path.join(out_dir, sample_id + "_minpath_in.txt")
    minpath_report = path.join(out_dir, sample_id + "_minpath_report.txt")
    minpath_details = path.join(out_dir, sample_id + "_minpath_details.txt")
    minpath_mps = path.join(out_dir, sample_id + "_minpath.mps")

    minpath_output = open(path.join(out_dir, sample_id + "_minpath_out.txt"),
                          "w")

    id_minpath_fh = open(minpath_in, "w")

    # Loop over all functions (which are the index labels in unstrat table).
    for func_id in unstrat_input.index.values:
        # Get count of each sequence in sample and write that sequence out
        # along with count if non-zero abundance.
        func_count = unstrat_input.loc[func_id, sample_id]

        # If 0 then skip.
        if func_count == 0:
            continue

        id_minpath_fh.write(func_id + "\t" + str(func_count) + "\n")

    id_minpath_fh.close()

    # Run MinPath on this sample.
    path2minpath = path.join(get_picrust_project_dir(), 'MinPath',
                                 'MinPath12hmp.py')

    minpath_cmd = path2minpath + " -any " + minpath_in + " -map " +\
                  minpath_map + " -report " + minpath_report +\
                  " -details " + minpath_details + " -mps " + minpath_mps

    system_call_check(minpath_cmd, print_out=print_opt,
                      stdout=minpath_output)

    # Read through MinPath report and keep track of pathways identified
    # to be present.
    path_present = identify_minpath_present(minpath_report)

    # Now read in details file and take abundance of pathway to be
    # mean of top 1/2 most abundant gene families.
    # Abundances of 0 will be added in for gene families not found.
    gf_abundances, gf_ids = parse_minpath_details(minpath_details, path_present)

    # Initialize series and dataframe that will contain pathway abundances.
    unstrat_abun = pd.Series()
    strat_abun = pd.DataFrame(columns=["pathway", "sequence", sample_id])
    strat_abun = strat_abun.set_index(["pathway", "sequence"])

    # Loop through all pathways present and get mean of 1/2 most abundant.
    for pathway in gf_abundances.keys():

        # Like HUMAnN2, sort enzyme reactions, take second half, and get
        # their mean abundance.

        # First get indices of sorted list.
        sorted_index = list(np.argsort(gf_abundances[pathway]))
        sorted_gf_abundances = [gf_abundances[pathway][i] for i in sorted_index]
        sorted_gf_ids = [gf_ids[pathway][i] for i in sorted_index]

        # Take second half of gene family abundances and ids lists.
        half_i = int(len(sorted_gf_abundances) / 2)
        gf_abundances_subset = sorted_gf_abundances[half_i:]
        gf_ids_subset = sorted_gf_ids[half_i:]

        # Take mean for unstratified pathway abundance.
        unstrat_abun[pathway] = sum(gf_abundances_subset)/len(gf_abundances_subset)

        # Get stratified pathway abundances by sequences.
        strat_path_abun = path_abun_by_seq(strat_input,
                                           gf_ids_subset,
                                           sum(gf_abundances_subset),
                                           unstrat_abun[pathway])
        # Remove rows that are all 0.
        strat_path_abun[strat_path_abun[sample_id] > 0]

        # Add pathway as new column.
        strat_path_abun["pathway"] = [pathway]*strat_path_abun.shape[0]

        strat_path_abun.set_index("pathway", append=True, inplace=True)

        # Changes levels of index labels.
        strat_path_abun = strat_path_abun.reorder_levels(["pathway",
                                                          "sequence"])

        strat_abun = pd.concat([strat_abun, strat_path_abun], levels=["pathway", "sequence"])

    # Return unstratified and stratified abundances.
    # Note that the stratified abundances are converted to a series.
    return([unstrat_abun, strat_abun[sample_id]])
