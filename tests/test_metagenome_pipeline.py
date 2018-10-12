#!/usr/bin/env python

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.0.2-b"

import unittest
from os import path
import pandas as pd
import biom
from tempfile import TemporaryDirectory
from picrust2.util import get_picrust_project_dir, biom_to_pandas_df
from picrust2.metagenome_pipeline import (run_metagenome_pipeline,
                                          norm_by_marker_copies,
                                          calc_weighted_nsti,
                                          id_rare_seqs)

# Set paths to test files.
test_dir_path = path.join(get_picrust_project_dir(), "tests", "test_data",
                          "metagenome_pipeline")

seqtab_biom = path.join(test_dir_path, "test_input_sequence_abun.biom")
seqtab_tsv = path.join(test_dir_path, "test_input_sequence_abun.tsv")
func_predict = path.join(test_dir_path, "test_predicted_func.tsv")
marker_predict = path.join(test_dir_path, "test_predicted_marker.tsv")

nsti_in_path = path.join(test_dir_path, "test_nsti_in.tsv")

exp_strat = path.join(test_dir_path, "metagenome_out",
                      "pred_metagenome_strat.tsv")

exp_unstrat = path.join(test_dir_path, "metagenome_out",
                        "pred_metagenome_unstrat.tsv")

exp_strat_rare = path.join(test_dir_path, "metagenome_out",
                           "pred_metagenome_strat_RARE.tsv")

exp_norm = path.join(test_dir_path, "metagenome_out", "seqtab_norm.tsv")

# Read in test inputs and expected files.
func_predict_in = pd.read_table(func_predict, sep="\t", index_col="sequence")
marker_predict_in = pd.read_table(marker_predict, sep="\t",
                                  index_col="sequence")

exp_strat_in = pd.read_table(exp_strat, sep="\t")
exp_strat_in = exp_strat_in.set_index(["function", "sequence"])

exp_strat_in_rare = pd.read_table(exp_strat_rare, sep="\t")
exp_strat_in_rare = exp_strat_in_rare.set_index(["function", "sequence"])


exp_unstrat_in = pd.read_table(exp_unstrat, sep="\t", index_col="function")

exp_norm_in = pd.read_table(exp_norm, sep="\t", index_col="sequence")

nsti_in = pd.read_table(nsti_in_path, sep="\t", index_col="sequence")

class metagenome_pipeline_test(unittest.TestCase):

    def test_full_pipeline_strat_tsv(self):
        '''Test that run_metagenome_pipeline works on tsv input seqtab.'''

        with TemporaryDirectory() as temp_dir:
            strat_out, unstrat_out = run_metagenome_pipeline(input_biom=seqtab_tsv,
                                                             function=func_predict,
                                                             marker=marker_predict,
                                                             max_nsti=2,
                                                             out_dir=temp_dir,
                                                             strat_out=True)

        pd.testing.assert_frame_equal(strat_out, exp_strat_in, check_like=True)

        pd.testing.assert_frame_equal(unstrat_out, exp_unstrat_in,
                                      check_like=True)

    def test_full_pipeline_unstrat_tsv_when_no_strat(self):
        '''Test that run_metagenome_pipeline works on tsv input seqtab when
        strat_out=False.'''

        with TemporaryDirectory() as temp_dir:
            strat_out, unstrat_out = run_metagenome_pipeline(input_biom=seqtab_tsv,
                                                             function=func_predict,
                                                             marker=marker_predict,
                                                             max_nsti=2,
                                                             out_dir=temp_dir,
                                                             strat_out=False)

        pd.testing.assert_frame_equal(unstrat_out, exp_unstrat_in,
                                      check_like=True)

    def test_full_pipeline_strat_tsv_2proc(self):
        '''Test that run_metagenome_pipeline works on tsv input seqtab and
        running on 2 processes.'''

        with TemporaryDirectory() as temp_dir:
            strat_out, unstrat_out = run_metagenome_pipeline(input_biom=seqtab_tsv,
                                                             function=func_predict,
                                                             marker=marker_predict,
                                                             max_nsti=2,
                                                             out_dir=temp_dir,
                                                             proc=2,
                                                             strat_out=True)
    
        pd.testing.assert_frame_equal(strat_out, exp_strat_in, check_like=True)

    def test_full_pipeline_strat_biom(self):
        '''Test that run_metagenome_pipeline create corrected stratified output
        on biom input seqtab.'''

        with TemporaryDirectory() as temp_dir:
            strat_out, unstrat_out = run_metagenome_pipeline(input_biom=seqtab_biom,
                                                             function=func_predict,
                                                             marker=marker_predict,
                                                             max_nsti=2,
                                                             out_dir=temp_dir,
                                                             strat_out=True)

        pd.testing.assert_frame_equal(strat_out, exp_strat_in, check_like=True)

    def test_full_pipeline_unstrat_biom_when_no_strat(self):
        '''Test that run_metagenome_pipeline create corrected unstratified
        output on biom input seqtab when strat_out=False.'''

        with TemporaryDirectory() as temp_dir:
            strat_out, unstrat_out = run_metagenome_pipeline(input_biom=seqtab_biom,
                                                             function=func_predict,
                                                             marker=marker_predict,
                                                             max_nsti=2,
                                                             out_dir=temp_dir,
                                                             strat_out=False)

        pd.testing.assert_frame_equal(unstrat_out, exp_unstrat_in,
                                      check_like=True)

    def test_full_pipeline_unstrat_biom(self):
        '''Test that run_metagenome_pipeline create corrected unstratified
        output on biom input seqtab.'''

        with TemporaryDirectory() as temp_dir:
            strat_out, unstrat_out = run_metagenome_pipeline(input_biom=seqtab_biom,
                                                             function=func_predict,
                                                             marker=marker_predict,
                                                             max_nsti=2,
                                                             out_dir=temp_dir,
                                                             strat_out=False)

        pd.testing.assert_frame_equal(unstrat_out, exp_unstrat_in,
                                      check_like=True)

    def test_norm_by_marker_copies(self):
        '''Test that expected normalized sequence abundance table generated.'''

        seqtab_in = biom_to_pandas_df(biom.load_table(seqtab_biom))

        # Get output index labels in same order as expected.
        seqtab_in = seqtab_in.reindex(exp_norm_in.index)

        test_norm = norm_by_marker_copies(input_seq_counts=seqtab_in,
                                          input_marker_num=marker_predict_in,
                                          norm_filename=None)

        # Test whether normalized table matches expected table.
        pd.testing.assert_frame_equal(test_norm, exp_norm_in, check_like=True)

    def test_weighted_nsti(self):
        '''Test that expected weighted NSTI values are calculated.'''

        weighted_nsti_out = calc_weighted_nsti(exp_norm_in, nsti_in)

        expected_weighted_nsti = { "samples": ["sample1", "sample2",
                                                "sample3"],
                                   "weighted_NSTI": [0.292857143, 0.348275862, 
                                                      0.436111111] }

        expected_weighted_nsti_df = pd.DataFrame.from_dict(expected_weighted_nsti)

        expected_weighted_nsti_df.set_index("samples", inplace=True)

        expected_weighted_nsti_df.index.name = None

        pd.testing.assert_frame_equal(weighted_nsti_out,
                                      expected_weighted_nsti_df,
                                      check_like=True)

class rare_seqs_test(unittest.TestCase):
    '''Checks that \"RARE\" category is being collapsed to correctly.'''

    def test_rare_4_reads(self):
        '''Check that correct sequences are identified as rare when a cut-off
        of 4 reads is used.'''

        seqtab_in = biom_to_pandas_df(biom.load_table(seqtab_biom))

        rare_seqs = id_rare_seqs(seqtab_in, 4, 1)

        self.assertSetEqual(set(rare_seqs), set(["2558860574", "extra"]))

    def test_rare_2_samp(self):
        '''Check that correct sequences are identified as rare when a cut-off
        of 2 samples is used.'''

        seqtab_in = biom_to_pandas_df(biom.load_table(seqtab_biom))

        rare_seqs = id_rare_seqs(seqtab_in, 1, 2)

        self.assertSetEqual(set(rare_seqs), set(["2558860574", "2571042244"]))

    def test_full_pipeline_strat_tsv_rare_category(self):
        '''Test that run_metagenome_pipeline works on tsv input seqtab and when
        rare seqs are collapsed into RARE category'''

        with TemporaryDirectory() as temp_dir:
            strat_out, unstrat_out = run_metagenome_pipeline(input_biom=seqtab_tsv,
                                                             function=func_predict,
                                                             marker=marker_predict,
                                                             max_nsti=2,
                                                             min_reads=4,
                                                             min_samples=2,
                                                             out_dir=temp_dir,
                                                             strat_out=True)

        pd.testing.assert_frame_equal(strat_out, exp_strat_in_rare,
                                      check_like=True)

        pd.testing.assert_frame_equal(unstrat_out, exp_unstrat_in,
                                      check_like=True)   


if __name__ == '__main__':
    unittest.main()
