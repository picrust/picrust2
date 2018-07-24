#!/usr/bin/env python

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.0.0-b.4"

import unittest
import pandas as pd
from os import path
from tempfile import TemporaryDirectory
from picrust2.run_minpath import (minpath_wrapper, run_minpath_pipeline,
                                  read_metagenome_input, unstrat_minpath,
                                  strat_minpath, strat_to_unstrat_counts)
from picrust2.util import get_picrust_project_dir

# Path to test directory.
test_dir_path = path.join(get_picrust_project_dir(), "tests")

in_metagenome_strat = path.join(test_dir_path, "test_data", "run_minpath",
                                "test_metagenome_out.tsv")

in_metagenome_unstrat = path.join(test_dir_path, "test_data", "run_minpath",
                                  "test_metagenome_unstrat_out.tsv")

exp_minpath_out_strat = path.join(test_dir_path, "test_data", "run_minpath",
                                  "expected_out_strat_path.tsv")

exp_minpath_out_unstrat = path.join(test_dir_path, "test_data", "run_minpath",
                                    "expected_out_unstrat_path.tsv")

map_ec2path_prokaryotic = path.join(get_picrust_project_dir(), "MinPath",
                                    "ec2metacyc_picrust_prokaryotic.txt")


class minpath_wrapper_tests(unittest.TestCase):
    """Tests for running MinPath on stratified and unstratified tables."""

    def test_strat_minpath_single_sample(self):
        '''Test running strat_minpath on single sample.'''

        strat_in, strat_tab = read_metagenome_input(in_metagenome_strat)

        with TemporaryDirectory() as temp_dir:
            path_abun = strat_minpath(sample_id="sample2",
                                      strat_input=strat_in[["function",
                                                            "sequence",
                                                            "sample2"]],
                                      minpath_map=map_ec2path_prokaryotic,
                                      out_dir=temp_dir)

        # Convert to pandas dataframe.
        unstrat_path_abun_df = pd.DataFrame([path_abun[0]]).transpose()
        strat_path_abun_df = pd.DataFrame([path_abun[1]]).transpose()

        # Set index title and column name.
        unstrat_path_abun_df.index.name = "pathway"
        unstrat_path_abun_df.columns = ["sample2"]

        strat_path_abun_df.reset_index(inplace=True)

        # Compare this predicted column to expected (after removing rows that
        # are 0).
        exp_path_unstrat = pd.read_csv(exp_minpath_out_unstrat, sep="\t",
                                    index_col="pathway")

        exp_path_unstrat_s2 = exp_path_unstrat[["sample2"]]

        exp_path_strat = pd.read_csv(exp_minpath_out_strat, sep="\t")

        # Sort stratified files (different versions can sort the output
        # slightly differently).
        strat_path_abun_df.sort_values(['pathway', 'sequence'], inplace=True)
        exp_path_strat.sort_values(['pathway', 'sequence'], inplace=True)

        # Reset index labels.
        strat_path_abun_df.reset_index(inplace=True, drop=True)
        exp_path_strat.reset_index(inplace=True, drop=True)

        exp_path_strat_s2 = exp_path_strat[["pathway", "sequence", "sample2"]]

        pd.testing.assert_frame_equal(exp_path_unstrat_s2, unstrat_path_abun_df,
                                      check_like=True)

        pd.testing.assert_frame_equal(exp_path_strat_s2, strat_path_abun_df,
                                      check_like=True)


    def test_unstrat_minpath_3_samples(self):
        '''Test running unstrat_minpath on single sample.'''

        strat_in, strat_tab = read_metagenome_input(in_metagenome_strat)

        unstrat_in = strat_to_unstrat_counts(strat_in)
        unstrat_in["function"] = unstrat_in.index

        with TemporaryDirectory() as temp_dir:

            path_abun = []

            for samp in ["sample1", "sample2", "sample3"]:
                path_abun += [unstrat_minpath(sample_id=samp,
                                              unstrat_input=unstrat_in[["function",
                                                                  samp]],
                                              minpath_map=map_ec2path_prokaryotic,
                                              out_dir=temp_dir)]

        # Convert to pandas dataframe.
        unstrat_path_abun_df = pd.DataFrame(path_abun).fillna(0).transpose()
        unstrat_path_abun_df.index.name = "pathway"
        unstrat_path_abun_df.columns = ["sample1", "sample2", "sample3"]

        # Compare this predicted column to expected (after removing rows that
        # are 0).
        exp_path_abun_unstrat = pd.read_csv(exp_minpath_out_unstrat, sep="\t",
                                    index_col="pathway")

        pd.testing.assert_frame_equal(exp_path_abun_unstrat,
                                      unstrat_path_abun_df,
                                      check_like=True)


class run_minpath_pipeline_tests(unittest.TestCase):
    """Tests for run_minpath_pipeline function."""

    def test_strat_pipeline_2_proc(self):
        '''Test running full pipeline over 2 processes.'''

        with TemporaryDirectory() as temp_dir:
            test_unstrat, test_strat = run_minpath_pipeline(inputfile=in_metagenome_strat,
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

    def test_unstrat_pipeline_2_proc(self):
        '''Test running full pipeline with unstratified table over 2
        processes.'''

        with TemporaryDirectory() as temp_dir:
            test_unstrat, test_strat = run_minpath_pipeline(inputfile=in_metagenome_unstrat,
                                                            mapfile=map_ec2path_prokaryotic,
                                                            proc=2,
                                                            out_dir=temp_dir)

        # Compare to expected pathway abundances.
        exp_path_abun_unstrat = pd.read_csv(exp_minpath_out_unstrat, sep="\t",
                                            index_col="pathway")

        test_unstrat.index.name = "pathway"

        pd.testing.assert_frame_equal(exp_path_abun_unstrat, test_unstrat,
                                      check_like=True)


if __name__ == '__main__':
    unittest.main()
