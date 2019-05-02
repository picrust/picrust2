#!/usr/bin/env python

from __future__ import division

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.1.3-b"

from os import path
import pandas as pd
from math import ceil
from joblib import Parallel, delayed
from picrust2.util import system_call_check, TemporaryDirectory

def castor_hsp_workflow(tree_path,
                        trait_table_path,
                        hsp_method,
                        chunk_size=500,
                        calc_nsti=False,
                        calc_ci=False,
                        check_input=False,
                        num_proc=1,
                        ran_seed=None):
    '''Runs full HSP workflow. Main purpose is to read in trait table and run
    HSP on subsets of column at a time to be more memory efficient. Will return
    a single table of predictions and also a table of CIs (if specified).'''

    # Read in trait table as pandas dataframe.
    trait_tab = pd.read_table(trait_table_path, sep="\t", index_col="assembly",
                              dtype={'assembly' : str})

    # Calculate NSTI values if option set.
    if calc_nsti:
        nsti_values = castor_nsti(tree_path, trait_tab.index.values)

    # Create output directory for writing trait table subsets.
    with TemporaryDirectory() as temp_dir:

        num_chunks = int(trait_tab.shape[1]) / (chunk_size + 1)

        # Get all table subsets and write to file.
        # Also keep list of all temporary filenames.
        file_subsets = []

        for i in range(ceil(num_chunks)):

            subset_file = path.join(temp_dir, "subset_tab_" + str(i))

            subset_tab = trait_tab.iloc[:, i * chunk_size:(i + 1) * chunk_size]

            subset_tab.to_csv(path_or_buf=subset_file,
                              index_label="assembly",
                              sep="\t")

            file_subsets.append(subset_file)

        castor_out_raw = Parallel(n_jobs=num_proc)(delayed(
                                    castor_hsp_wrapper)(tree_path,
                                                        trait_in,
                                                        hsp_method,
                                                        calc_ci,
                                                        check_input,
                                                        ran_seed)
                                    for trait_in in file_subsets)

    # Get lists of predictions and CIs for all chunks.
    predict_out_chunks = []
    ci_out_chunks = []

    for i in range(len(castor_out_raw)):
        predict_out_chunks.append(castor_out_raw[i][0])
        ci_out_chunks.append(castor_out_raw[i][1])

    predict_out_combined = pd.concat(predict_out_chunks, axis=1, sort=True)

    # Add NSTI as column as well if option specified.
    if calc_nsti:
        predict_out_combined = pd.concat([predict_out_combined, nsti_values],
                                         axis=1, sort=True)

    ci_out_combined = None

    if calc_ci:
        ci_out_combined = pd.concat(ci_out_chunks, axis=1, sort=True)

    return(predict_out_combined, ci_out_combined)
    

def castor_hsp_wrapper(tree_path, trait_tab, hsp_method, calc_ci=False,
                       check_input=False, ran_seed=None):
    '''Wrapper for making system calls to castor_hsp.py Rscript.'''

    castor_hsp_script = path.join(path.dirname(path.abspath(__file__)),
                                  'Rscripts', 'castor_hsp.R')

    # Need to format boolean setting as string for R to read in as argument.
    if calc_ci:
        calc_ci_setting = "TRUE"
    else:
        calc_ci_setting = "FALSE"

    if check_input:
        check_input_setting = "TRUE"
    else:
        check_input_setting = "FALSE"

    # Create temporary directory for writing output files of castor_hsp.R
    with TemporaryDirectory() as temp_dir:

        output_count_path = path.join(temp_dir, "predicted_counts.txt")
        output_ci_path = path.join(temp_dir, "predicted_ci.txt")

        hsp_cmd = " ".join(["Rscript",
                            castor_hsp_script,
                            tree_path,
                            trait_tab,
                            hsp_method,
                            calc_ci_setting,
                            check_input_setting,
                            output_count_path,
                            output_ci_path,
                            str(ran_seed)])

        # Run castor_hsp.R
        system_call_check(hsp_cmd)

        # Load the output into Table objects
        try:
            asr_table = pd.read_table(filepath_or_buffer=output_count_path,
                                  sep="\t", index_col="sequence")
        except IOError:
            raise ValueError("Cannot read in expected output file" +
                            output_ci_path)

        if calc_ci:
            asr_ci_table = pd.read_table(filepath_or_buffer=output_ci_path,
                                  sep="\t", index_col="sequence")
        else:
            asr_ci_table = None

    # Return list with predicted counts and CIs.
    return [asr_table, asr_ci_table]


def castor_nsti(tree_path,
                known_tips):
    '''Will calculate distance from each study sequence to the closest
    reference sequence. Takes in the path to treefile and the known tips
    (i.e. the rownames in the trait table - the reference genome ids).'''
    castor_nsti_script = path.join(path.dirname(path.abspath(__file__)),
                                   'Rscripts', 'castor_nsti.R')

    # Create temporary directory for working in.
    with TemporaryDirectory() as temp_dir:

        # Output known tip names to temp file
        # (note this object is a numpy.ndarray)
        known_tips_out = path.join(temp_dir, "known_tips.txt")
        known_tips.tofile(known_tips_out, sep="\n")

        nsti_tmp_out = path.join(temp_dir, "nsti_out.txt")

        # Run Rscript.
        system_call_check(" ".join(["Rscript",
                                    castor_nsti_script,
                                    tree_path,
                                    known_tips_out,
                                    nsti_tmp_out]))

        # Read in calculated NSTI values.
        nsti_out = pd.read_table(nsti_tmp_out, sep="\t", index_col="sequence")

    # Make sure that the table has the correct number of rows.
    if len(known_tips) != nsti_out.shape[0]:
        ValueError("Number of rows in returned NSTI table is incorrect.")

    return(nsti_out)

