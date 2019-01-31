#!/usr/bin/env python

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.0.4-b"

import unittest
import pandas as pd
from os import path
from tempfile import TemporaryDirectory
from picrust2.default import default_regroup_map, default_pathway_map
from picrust2.pathway_pipeline import (minpath_wrapper,
                                       pathway_pipeline,
                                       read_metagenome_input,
                                       unstrat_pathway_levels,
                                       basic_strat_pathway_levels,
                                       strat_to_unstrat_counts,
                                       PathwaysDatabase,
                                       regroup_func_ids)

# Path to test directory.
test_dir_path = path.join(path.dirname(path.abspath(__file__)), "test_data",
                          "pathway_pipeline")

# Paths to input files and expected outputs.
in_metagenome_strat = path.join(test_dir_path, "test_input_ec_strat.tsv")
in_metagenome_unstrat = path.join(test_dir_path, "test_input_ec_unstrat.tsv")

exp_abun_strat_file = path.join(test_dir_path, "exp_path_abun_strat.tsv")
exp_abun_unstrat_file = path.join(test_dir_path, "exp_path_abun_unstrat.tsv")
exp_cov_unstrat_file = path.join(test_dir_path, "exp_path_cov_unstrat.tsv")


in_metagenome_unstrat_per_seq = path.join(test_dir_path,
                                          "per_seq_contrib_input",
                                          "pred_metagenome_unstrat.tsv")

in_metagenome_strat_per_seq = path.join(test_dir_path,
                                        "per_seq_contrib_input",
                                        "pred_metagenome_strat.tsv")

in_per_seq_abun = path.join(test_dir_path, "per_seq_contrib_input",
                            "seqtab_norm.tsv")

in_per_seq_func = path.join(test_dir_path, "per_seq_contrib_input",
                            "per_seq_func.txt")

exp_abun_unstrat_per_genome_file = path.join(test_dir_path, "per_seq_contrib_input", "humann2_run", "humann2_unstrat_pathabun_output.tsv")
exp_cov_unstrat_per_genome_file = path.join(test_dir_path, "per_seq_contrib_input", "humann2_run", "humann2_unstrat_pathcov_output.tsv")
exp_abun_strat_per_genome_file = path.join(test_dir_path, "per_seq_contrib_input", "humann2_run", "humann2_strat_pathabun_output_seq_norm.tsv")
exp_cov_strat_per_genome_file = path.join(test_dir_path, "per_seq_contrib_input", "humann2_run", "humann2_strat_pathcov_output_seq_norm.tsv")

class run_minpath_tests(unittest.TestCase):
    """Tests for running MinPath on stratified and unstratified tables."""

    def test_strat_default_pipeline(self):
        '''Test running strat_minpath default pipeline. Make sure that
        community wide stratified abundances are calculated correctly and
        that unstratified abundances are right.'''

        with TemporaryDirectory() as temp_dir:
            unstrat_path_abun_df, unstrat_path_cov_df, strat_path_abun_df, strat_cov = pathway_pipeline(in_metagenome_strat,
                                                                                                        default_pathway_map,
                                                                                                        proc=1,
                                                                                                        out_dir=temp_dir,
                                                                                                        run_minpath=True,
                                                                                                        coverage=True,
                                                                                                        regroup_mapfile=default_regroup_map,
                                                                                                        gap_fill_on=True,
                                                                                                        per_sequence_contrib=False,
                                                                                                        print_cmds=False)


        # Compare these predicted tables to expected tables.
        exp_abun_unstrat = pd.read_csv(exp_abun_unstrat_file, sep="\t",
                                       index_col="pathway")

        exp_cov_unstrat = pd.read_csv(exp_cov_unstrat_file, sep="\t",
                                       index_col="pathway")


        exp_abun_strat = pd.read_csv(exp_abun_strat_file, sep="\t")

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
            unstrat_path_abun_df, unstrat_path_cov_df, strat_abun, strat_cov = pathway_pipeline(in_metagenome_unstrat,
                                                                                                default_pathway_map,
                                                                                                proc=1,
                                                                                                out_dir=temp_dir,
                                                                                                run_minpath=True,
                                                                                                coverage=True,
                                                                                                regroup_mapfile=default_regroup_map,
                                                                                                gap_fill_on=True,
                                                                                                per_sequence_contrib=False,
                                                                                                print_cmds=False)

        # Compare these predicted tables to expected tables.
        exp_abun_unstrat = pd.read_csv(exp_abun_unstrat_file, sep="\t",
                                       index_col="pathway")

        exp_cov_unstrat = pd.read_csv(exp_cov_unstrat_file, sep="\t",
                                       index_col="pathway")

        pd.testing.assert_frame_equal(exp_abun_unstrat, unstrat_path_abun_df,
                                      check_like=True, check_less_precise=True)

        pd.testing.assert_frame_equal(exp_cov_unstrat, unstrat_path_cov_df,
                                      check_like=True, check_less_precise=True)


    def test_strat_per_genome_pipeline(self):
        '''Test running strat_minpath default pipeline. Make sure that
        per genome contributions are correct (per_sequence_contrib set).'''

        with TemporaryDirectory() as temp_dir:
            unstrat_path_abun_df, unstrat_path_cov_df, strat_path_abun_df, strat_path_cov_df = pathway_pipeline(in_metagenome_unstrat_per_seq,
                                                                                                                default_pathway_map,
                                                                                                                proc=1,
                                                                                                                out_dir=temp_dir,
                                                                                                                run_minpath=True,
                                                                                                                coverage=True,
                                                                                                                regroup_mapfile=default_regroup_map,
                                                                                                                gap_fill_on=True,
                                                                                                                per_sequence_contrib=True,
                                                                                                                per_sequence_abun=in_per_seq_abun,
                                                                                                                per_sequence_function=in_per_seq_func,
                                                                                                                print_cmds=False)


        # Compare these predicted tables to expected tables.
        exp_abun_unstrat = pd.read_csv(exp_abun_unstrat_per_genome_file,
                                       sep="\t", index_col="pathway")

        exp_cov_unstrat = pd.read_csv(exp_cov_unstrat_per_genome_file,
                                      sep="\t", index_col="pathway")
        exp_abun_strat = pd.read_csv(exp_abun_strat_per_genome_file, sep="\t")
        exp_cov_strat = pd.read_csv(exp_cov_strat_per_genome_file, sep="\t")

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

    def test_strat_per_genome_pipeline_strat_input(self):
        '''Test running strat_minpath default pipeline. Make sure that
        per genome contributions are correct (per_sequence_contrib set).
        In this case the input is a stratified table.'''

        with TemporaryDirectory() as temp_dir:
            unstrat_path_abun_df, unstrat_path_cov_df, strat_path_abun_df, strat_path_cov_df = pathway_pipeline(in_metagenome_strat_per_seq,
                                                                                                                default_pathway_map,
                                                                                                                proc=1,
                                                                                                                out_dir=temp_dir,
                                                                                                                run_minpath=True,
                                                                                                                coverage=True,
                                                                                                                regroup_mapfile=default_regroup_map,
                                                                                                                gap_fill_on=True,
                                                                                                                per_sequence_contrib=True,
                                                                                                                per_sequence_abun=in_per_seq_abun,
                                                                                                                per_sequence_function=in_per_seq_func,
                                                                                                                print_cmds=False)

        # Compare these predicted tables to expected tables.
        exp_abun_unstrat = pd.read_csv(exp_abun_unstrat_per_genome_file,
                                       sep="\t", index_col="pathway")

        exp_cov_unstrat = pd.read_csv(exp_cov_unstrat_per_genome_file,
                                      sep="\t", index_col="pathway")
        exp_abun_strat = pd.read_csv(exp_abun_strat_per_genome_file, sep="\t")
        exp_cov_strat = pd.read_csv(exp_cov_strat_per_genome_file, sep="\t")

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
