#!/usr/bin/env python

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.0.0-b.9"

import unittest
import pandas as pd
from os import path
from tempfile import TemporaryDirectory
from picrust2.default import default_regroup_map, default_pathway_map
from picrust2.run_minpath import (minpath_wrapper, run_minpath_pipeline,
                                  read_metagenome_input, unstrat_minpath,
                                  strat_minpath, strat_to_unstrat_counts,
                                  PathwaysDatabase, regroup_func_ids)
from picrust2.util import get_picrust_project_dir

# Path to test directory.
test_dir_path = path.join(get_picrust_project_dir(), "tests", "test_data",
                          "run_minpath")

# Paths to input files and expected outputs.
in_metagenome_strat = path.join(test_dir_path, "test_input_ec_strat.tsv")
in_metagenome_strat2 = path.join(test_dir_path, "test_input_ec_strat2.tsv")
in_metagenome_unstrat = path.join(test_dir_path, "test_input_ec_unstrat.tsv")

exp_minpath_abun_strat = path.join(test_dir_path, "exp_path_abun_strat.tsv")
exp_minpath_abun_strat_per_genome = path.join(test_dir_path, "exp_path_abun_strat_per-genome.tsv")
exp_minpath_abun_unstrat = path.join(test_dir_path, "exp_path_abun_unstrat.tsv")
exp_minpath_cov_strat_per_genome = path.join(test_dir_path, "exp_path_cov_strat_per-genome.tsv")
exp_minpath_cov_unstrat = path.join(test_dir_path, "exp_path_cov_unstrat.tsv")

class run_minpath_tests(unittest.TestCase):
    """Tests for running MinPath on stratified and unstratified tables."""

    def test_strat_default_pipeline(self):
        '''Test running strat_minpath default pipeline. Make sure that
        community wide stratified abundances are calculated correctly and
        that unstratified abundances are right.'''

        with TemporaryDirectory() as temp_dir:
            unstrat_path_abun_df, unstrat_path_cov_df, strat_path_abun_df, strat_cov = run_minpath_pipeline(in_metagenome_strat2,
                                                                                                 default_pathway_map,
                                                                                                 proc=1,
                                                                                                 out_dir=temp_dir,
                                                                                                 regroup_mapfile=default_regroup_map,
                                                                                                 gap_fill=True,
                                                                                                 per_sequence_contrib=False,
                                                                                                 print_cmds=False)


        # Compare these predicted tables to expected tables.
        exp_abun_unstrat = pd.read_csv(exp_minpath_abun_unstrat, sep="\t",
                                       index_col="pathway")

        exp_cov_unstrat = pd.read_csv(exp_minpath_cov_unstrat, sep="\t",
                                       index_col="pathway")


        exp_abun_strat = pd.read_csv(exp_minpath_abun_strat, sep="\t")

        # Sort stratified files (different versions can sort the output
        # slightly differently).
        strat_path_abun_df.sort_values(['pathway', 'sequence'], inplace=True)
        exp_abun_strat.sort_values(['pathway', 'sequence'], inplace=True)

        # Reset index labels.
        exp_abun_strat.reset_index(drop=True, inplace=True)
        strat_path_abun_df.reset_index(drop=True, inplace=True)

        pd.testing.assert_frame_equal(exp_abun_unstrat, unstrat_path_abun_df,
                                      check_like=True, check_less_precise=True)

        pd.testing.assert_frame_equal(exp_cov_unstrat, unstrat_path_cov_df,
                                      check_like=True, check_less_precise=True)

        # Check with less precision here since the HUMAnN2 output that is used
        # as expected abundances are not rounded.
        pd.testing.assert_frame_equal(exp_abun_strat, strat_path_abun_df,
                                      check_like=True, check_less_precise=True)


    def test_unstrat_default_pipeline(self):
        '''Test running default pipeline on unstratified input table.'''

        with TemporaryDirectory() as temp_dir:
            unstrat_path_abun_df, unstrat_path_cov_df, strat_abun, strat_cov = run_minpath_pipeline(in_metagenome_unstrat,
                                                                                                 default_pathway_map,
                                                                                                 proc=1,
                                                                                                 out_dir=temp_dir,
                                                                                                 regroup_mapfile=default_regroup_map,
                                                                                                 gap_fill=True,
                                                                                                 per_sequence_contrib=False,
                                                                                                 print_cmds=False)

        # Compare these predicted tables to expected tables.
        exp_abun_unstrat = pd.read_csv(exp_minpath_abun_unstrat, sep="\t",
                                       index_col="pathway")

        exp_cov_unstrat = pd.read_csv(exp_minpath_cov_unstrat, sep="\t",
                                       index_col="pathway")

        pd.testing.assert_frame_equal(exp_abun_unstrat, unstrat_path_abun_df,
                                      check_like=True, check_less_precise=True)

        pd.testing.assert_frame_equal(exp_cov_unstrat, unstrat_path_cov_df,
                                      check_like=True, check_less_precise=True)


    def test_strat_per_genome_pipeline(self):
        '''Test running strat_minpath default pipeline. Make sure that
        per genome contributions are correct (per_sequence_contrib set).'''

        with TemporaryDirectory() as temp_dir:
            unstrat_path_abun_df, unstrat_path_cov_df, strat_path_abun_df, strat_path_cov_df = run_minpath_pipeline(in_metagenome_strat,
                                                                                                 default_pathway_map,
                                                                                                 proc=1,
                                                                                                 out_dir=temp_dir,
                                                                                                 regroup_mapfile=default_regroup_map,
                                                                                                 gap_fill=True,
                                                                                                 per_sequence_contrib=True,
                                                                                                 print_cmds=False)


        # Compare these predicted tables to expected tables.
        exp_abun_unstrat = pd.read_csv(exp_minpath_abun_unstrat, sep="\t",
                                       index_col="pathway")

        exp_cov_unstrat = pd.read_csv(exp_minpath_cov_unstrat, sep="\t",
                                       index_col="pathway")
        exp_abun_strat = pd.read_csv(exp_minpath_abun_strat_per_genome, sep="\t")
        exp_cov_strat = pd.read_csv(exp_minpath_cov_strat_per_genome, sep="\t")

        # Sort stratified files (different versions can sort the output
        # slightly differently).
        strat_path_abun_df.sort_values(['pathway', 'sequence'], inplace=True)
        exp_abun_strat.sort_values(['pathway', 'sequence'], inplace=True)
        strat_path_cov_df.sort_values(['pathway', 'sequence'], inplace=True)
        exp_cov_strat.sort_values(['pathway', 'sequence'], inplace=True)

        # Reset index labels.
        exp_abun_strat.reset_index(drop=True, inplace=True)
        strat_path_abun_df.reset_index(drop=True, inplace=True)

        exp_cov_strat.reset_index(drop=True, inplace=True)
        strat_path_cov_df.reset_index(drop=True, inplace=True)

        pd.testing.assert_frame_equal(exp_abun_unstrat, unstrat_path_abun_df,
                                      check_like=True, check_less_precise=True)

        pd.testing.assert_frame_equal(exp_cov_unstrat, unstrat_path_cov_df,
                                      check_like=True, check_less_precise=True)

        pd.testing.assert_frame_equal(exp_abun_strat, strat_path_abun_df,
                                      check_like=True, check_less_precise=True)

        pd.testing.assert_frame_equal(exp_cov_strat, strat_path_cov_df,
                                      check_like=True, check_less_precise=True)


if __name__ == '__main__':
    unittest.main()
