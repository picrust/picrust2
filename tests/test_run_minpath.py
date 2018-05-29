#!/usr/bin/env python

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.0.0-b.2"

import unittest
import pandas as pd
from os import path
from tempfile import TemporaryDirectory
from picrust2.run_minpath import (minpath_wrapper, run_minpath_pipeline,
                                  read_strat_genes)
from picrust2.util import get_picrust_project_dir

# Path to test directory.
test_dir_path = path.join(get_picrust_project_dir(), "tests")

in_metagenome_abun = path.join(test_dir_path, "test_data", "run_minpath",
                               "test_metagenome_out.tsv")

exp_minpath_out_strat = path.join(test_dir_path, "test_data", "run_minpath",
                                  "expected_out_strat_path.tsv")

exp_minpath_out_unstrat = path.join(test_dir_path, "test_data", "run_minpath",
                                    "expected_out_unstrat_path.tsv")

map_ec2path_prokaryotic = path.join(get_picrust_project_dir(), "MinPath",
                                    "ec2metacyc_picrust_prokaryotic.txt")


class minpath_wrapper_tests(unittest.TestCase):
    """Tests for minpath_wrapper function."""

    def test_minpath_wrapper_single_sample(self):
        '''Test running minpath_wrapper on single sample.'''

        strat_in = read_strat_genes(in_metagenome_abun)

        with TemporaryDirectory() as temp_dir:
            path_abun = minpath_wrapper(sample_id="sample2",
                                        strat_input=strat_in[["function",
                                                              "sequence",
                                                              "sample2"]],
                                        minpath_map=map_ec2path_prokaryotic,
                                        out_dir=temp_dir)

        # Convert to pandas dataframe.
        unstrat_path_abun_df = pd.DataFrame([path_abun[0]]).transpose()

        #  Set index title and column name.
        unstrat_path_abun_df.index.name = "pathway"
        unstrat_path_abun_df.columns = ["sample2"]

        # Compare this predicted column to expected (after removing rows that
        # are 0).
        exp_path_abun = pd.read_csv(exp_minpath_out_unstrat, sep="\t",
                                    index_col="pathway")

        exp_path_abun_s2 = exp_path_abun[["sample2"]]

        pd.testing.assert_frame_equal(exp_path_abun_s2, unstrat_path_abun_df,
                                      check_like=True)


class run_minpath_pipeline_tests(unittest.TestCase):
    """Tests for run_minpath_pipeline function."""

    def test_basic_pipeline_2_proc(self):
        '''Test running full pipeline over 2 processes.'''

        with TemporaryDirectory() as temp_dir:
            test_unstrat, test_strat = run_minpath_pipeline(inputfile=in_metagenome_abun,
                                                            mapfile=map_ec2path_prokaryotic,
                                                            proc=2,
                                                            out_dir=temp_dir)

        # Compare to expected pathway abundances.
        exp_path_abun_strat = pd.read_csv(exp_minpath_out_strat, sep="\t")

        exp_path_abun_unstrat = pd.read_csv(exp_minpath_out_unstrat, sep="\t",
                                            index_col="pathway")

        test_unstrat.index.name = "pathway"

        # Sort stratified files (different versions can sort the output slightly differently).
        test_strat.sort_values(['pathway', 'sequence'], inplace=True)
        exp_path_abun_strat.sort_values(['pathway', 'sequence'], inplace=True)

        # Reset index labels.
        test_strat.reset_index(inplace=True, drop=True)
        exp_path_abun_strat.reset_index(inplace=True, drop=True)

        pd.testing.assert_frame_equal(exp_path_abun_unstrat, test_unstrat,
                                      check_like=True)
        pd.testing.assert_frame_equal(exp_path_abun_strat, test_strat,
                                      check_like=True)


if __name__ == '__main__':
    unittest.main()
