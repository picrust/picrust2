#!/usr/bin/env python

from __future__ import division

__license__ = "GPL"
__version__ = "2-alpha.8"

from os import remove, path
import pandas as pd
from picrust2.util import (system_call_check, get_picrust_project_dir,
                           generate_temp_filename)

def castor_hsp_wrapper(tree_path,
                       trait_table_path,
                       hsp_method,
                       calc_nsti=False,
                       calc_ci=False,
                       check_input=False,
                       num_cores=1,
                       rds_outfile=None,
                       ran_seed=None,
                       HALT_EXEC=False):
                       
    '''Runs the castor_hsp.R Rscript and returns predicted counts and
    confidence interval (if applicable) tables as pandas dataframes'''

    castor_hsp_script_fp = path.join(get_picrust_project_dir(), 'picrust2', 
                                     'Rscripts', 'castor_hsp.R')

    tmp_output_count_path = generate_temp_filename()
    tmp_output_ci_path = generate_temp_filename()

    # Need to format boolean setting as string for R to read in as argument.
    if calc_nsti:
        calc_nsti_setting = "TRUE"
    else:
        calc_nsti_setting = "FALSE"

    if calc_ci:
        calc_ci_setting = "TRUE"
    else:
        calc_ci_setting = "FALSE"

    if check_input:
        check_input_setting = "TRUE"
    else:
        check_input_setting = "FALSE"

    if rds_outfile:
        write_rds = "TRUE"
    else:
        write_rds = "FALSE"

    hsp_cmd = " ".join(["Rscript",
                        castor_hsp_script_fp,
                        tree_path,
                        trait_table_path,
                        hsp_method,
                        calc_nsti_setting,
                        calc_ci_setting,
                        check_input_setting,
                        str(num_cores),
                        tmp_output_count_path,
                        tmp_output_ci_path,
                        str(rds_outfile),
                        str(ran_seed),
                        write_rds])

    # Run castor_hsp.R here
    result = system_call_check(hsp_cmd)

    # Load the output into Table objects
    try:
        asr_table = pd.read_table(filepath_or_buffer=tmp_output_count_path,
                                  sep="\t", index_col="sequence")
    except IOError:
        raise ValueError("Cannot read in expected output file" +
                         tmp_output_count_path)

    remove(tmp_output_count_path)

    if calc_ci:
        asr_ci_table = pd.read_table(filepath_or_buffer=tmp_output_ci_path,
                                  sep="\t", index_col="sequence")
        remove(tmp_output_ci_path)
    else:
        asr_ci_table = None

    return asr_table, asr_ci_table


def castor_hsp_loocv_wrapper(tree_path,
                             trait_table_path,
                             tips_path,
                             hsp_method,
                             expected_out_path,
                             predicted_out_path,
                             metrics_out_path,
                             num_cores=1,
                             HALT_EXEC=False):
                       
    '''Runs the castor_hsp_loocv.R Rscript and writes out result tables'''
    castor_loocv_hsp_script_fp = path.join(get_picrust_project_dir(),
                                           'picrust2', 'Rscripts',
                                           'castor_hsp_loocv.R')

    loocv_cmd = " ".join(["Rscript",
                          castor_loocv_hsp_script_fp,
                          tree_path,
                          trait_table_path,
                          tips_path,
                          hsp_method,
                          expected_out_path,
                          predicted_out_path,
                          metrics_out_path,
                          str(num_cores)])

    # Run castor_hsp_loocv.R here
    result = system_call_check(loocv_cmd)
