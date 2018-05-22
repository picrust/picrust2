#!/usr/bin/env python

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.0.0-b.1"

import unittest
from os import path
import pandas as pd
import biom
from tempfile import TemporaryDirectory
from picrust2.util import get_picrust_project_dir, biom_to_pandas_df
from picrust2.metagenome_pipeline import (run_metagenome_pipeline,
                                          norm_by_marker_copies)

# Set paths to test files.
test_dir_path = path.join(get_picrust_project_dir(), "tests", "test_data",
                          "metagenome_pipeline")

seqtab_biom = path.join(test_dir_path, "test_input_sequence_abun.biom")
seqtab_tsv = path.join(test_dir_path, "test_input_sequence_abun.tsv")
func_predict = path.join(test_dir_path, "test_predicted_func.tsv")
marker_predict = path.join(test_dir_path, "test_predicted_marker.tsv")

exp_strat = path.join(test_dir_path, "metagenome_out",
                      "pred_metagenome_strat.tsv")

exp_unstrat = path.join(test_dir_path, "metagenome_out",
                        "pred_metagenome_unstrat.tsv")

exp_norm = path.join(test_dir_path, "metagenome_out", "seqtab_norm.tsv")

# Read in test inputs and expected files.
func_predict_in = pd.read_table(func_predict, sep="\t", index_col="sequence")
marker_predict_in = pd.read_table(marker_predict, sep="\t",
                                  index_col="sequence")

exp_strat_in = pd.read_table(exp_strat, sep="\t")

exp_unstrat_in = pd.read_table(exp_unstrat, sep="\t", index_col="function")

exp_norm_in = pd.read_table(exp_norm, sep="\t", index_col="sequence")


class metagenome_pipeline_test(unittest.TestCase):

    def test_full_pipeline_strat_tsv(self):
        '''Test that run_metagenome_pipeline works on tsv input seqtab.'''

        with TemporaryDirectory() as temp_dir:
            strat_out, unstrat_out = run_metagenome_pipeline(input_biom=seqtab_tsv,
                                                             function=func_predict,
                                                             marker=marker_predict,
                                                             out_dir=temp_dir)

        # Need to reset index names since these aren't in output files.
        strat_out.index = range(30)

        pd.testing.assert_frame_equal(strat_out, exp_strat_in)

    def test_full_pipeline_unstrat_tsv(self):
        '''Test that run_metagenome_pipeline works on tsv input seqtab.'''

        with TemporaryDirectory() as temp_dir:
            strat_out, unstrat_out = run_metagenome_pipeline(input_biom=seqtab_tsv,
                                                             function=func_predict,
                                                             marker=marker_predict,
                                                             out_dir=temp_dir)

        pd.testing.assert_frame_equal(unstrat_out, exp_unstrat_in)

    def test_full_pipeline_strat_tsv_2proc(self):
        '''Test that run_metagenome_pipeline works on tsv input seqtab and
        running on 2 processes.'''

        with TemporaryDirectory() as temp_dir:
            strat_out, unstrat_out = run_metagenome_pipeline(input_biom=seqtab_tsv,
                                                             function=func_predict,
                                                             marker=marker_predict,
                                                             out_dir=temp_dir,
                                                             proc=2)

        # Need to reset index names since these aren't in output files.
        strat_out.index = range(30)
    
        pd.testing.assert_frame_equal(strat_out, exp_strat_in)

    def test_full_pipeline_strat_biom(self):
        '''Test that run_metagenome_pipeline works on tsv input seqtab.'''

        with TemporaryDirectory() as temp_dir:
            strat_out, unstrat_out = run_metagenome_pipeline(input_biom=seqtab_biom,
                                                             function=func_predict,
                                                             marker=marker_predict,
                                                             out_dir=temp_dir)
        # Need to reset index names since these aren't in output files.
        strat_out.index = range(30)

        pd.testing.assert_frame_equal(strat_out, exp_strat_in)

    def test_full_pipeline_unstrat_biom(self):
        '''Test that run_metagenome_pipeline works on tsv input seqtab.'''

        with TemporaryDirectory() as temp_dir:
            strat_out, unstrat_out = run_metagenome_pipeline(input_biom=seqtab_biom,
                                                             function=func_predict,
                                                             marker=marker_predict,
                                                             out_dir=temp_dir)

        pd.testing.assert_frame_equal(unstrat_out, exp_unstrat_in)

    def test_norm_by_marker_copies(self):
        '''Test that expected normalized sequence abundance table generated.'''

        seqtab_in = biom_to_pandas_df(biom.load_table(seqtab_biom))

        # Get output index labels in same order as expected.
        seqtab_in = seqtab_in.reindex(exp_norm_in.index)

        test_norm = norm_by_marker_copies(input_seq_counts=seqtab_in,
                                          input_marker_num=marker_predict_in,
                                          norm_filename=None)

        # Test whether normalized table matches expected table.
        pd.testing.assert_frame_equal(test_norm, exp_norm_in)


if __name__ == '__main__':
    unittest.main()
