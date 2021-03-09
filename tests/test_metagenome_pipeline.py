#!/usr/bin/env python

__copyright__ = "Copyright 2018-2020, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.3.0-b"

import unittest
from os import path
import pandas as pd
import biom
from picrust2.util import TemporaryDirectory
from picrust2.metagenome_pipeline import (run_metagenome_pipeline,
                                          norm_by_marker_copies,
                                          calc_weighted_nsti,
                                          id_rare_seqs,
                                          drop_tips_by_nsti)

# Set paths to test files.
test_dir_path = path.join(path.dirname(path.abspath(__file__)), "test_data",
                          "metagenome_pipeline")

seqtab_biom = path.join(test_dir_path, "test_input_sequence_abun.biom")
seqtab_tsv = path.join(test_dir_path, "test_input_sequence_abun.tsv.gz")
seqtab_msf = path.join(test_dir_path, "test_input_sequence_abun.msf")

func_predict = path.join(test_dir_path, "test_predicted_func.tsv.gz")
marker_predict = path.join(test_dir_path, "test_predicted_marker.tsv.gz")

nsti_in_path = path.join(test_dir_path, "test_nsti_in.tsv.gz")

seqtab_tsv_simple = path.join(test_dir_path, "test_sequence_abun_simple.txt.gz")

func_simple_in = path.join(test_dir_path, "test_predicted_func_simple.txt.gz")

marker_simple_in = path.join(test_dir_path, "test_predicted_marker_simple.txt.gz")

exp_strat_simple = path.join(test_dir_path, "metagenome_out",
                             "expected_metagenome_contrib.txt.gz")

exp_strat_wide = path.join(test_dir_path, "metagenome_out",
                           "pred_metagenome_strat.tsv.gz")

exp_unstrat_simple = path.join(test_dir_path, "metagenome_out",
                               "pred_metagenome_unstrat_simple.tsv.gz")

exp_unstrat = path.join(test_dir_path, "metagenome_out",
                        "pred_metagenome_unstrat.tsv.gz")

exp_unstrat = path.join(test_dir_path, "metagenome_out",
                        "pred_metagenome_unstrat.tsv.gz")

exp_strat_simple_rare = path.join(test_dir_path, "metagenome_out",
                                  "expected_metagenome_contrib_rare.txt.gz")

exp_norm = path.join(test_dir_path, "metagenome_out", "seqtab_norm.tsv.gz")

# Read in test inputs and expected files.
func_predict_in = pd.read_csv(func_predict, sep="\t", dtype={'sequence': str})
func_predict_in.set_index('sequence', drop=True, inplace=True)

marker_predict_in = pd.read_csv(marker_predict, sep="\t",
                                dtype={'sequence': str})
marker_predict_in.set_index('sequence', drop=True, inplace=True)

exp_strat_simple_in = pd.read_csv(exp_strat_simple, sep="\t",
                                       dtype={'sequence': str, 'taxon': str,
                                              'function': str, 'sample': str})

exp_strat_simple_rare_in = pd.read_csv(exp_strat_simple_rare, sep="\t",
                                       dtype={'sequence': str, 'taxon': str,
                                              'function': str, 'sample': str})

exp_strat_wide_in = pd.read_csv(exp_strat_wide, sep="\t",
                                dtype={'sequence': str, 'function': str})
exp_strat_wide_in = exp_strat_wide_in.set_index(["function", "sequence"])

exp_unstrat_in = pd.read_csv(exp_unstrat, sep="\t", dtype={'function': str})
exp_unstrat_in.set_index('function', drop=True, inplace=True)

exp_unstrat_simple_in = pd.read_csv(exp_unstrat_simple, sep="\t",
                                    dtype={'function': str})
exp_unstrat_simple_in.set_index('function', drop=True, inplace=True)

exp_norm_in = pd.read_csv(exp_norm, sep="\t", dtype={'normalized': str})
exp_norm_in.set_index('normalized', drop=True, inplace=True)

nsti_in = pd.read_csv(nsti_in_path, sep="\t", dtype={'sequence': str})
nsti_in.set_index('sequence', drop=True, inplace=True)


class metagenome_pipeline_test(unittest.TestCase):

    def test_full_pipeline_strat_tsv(self):
        '''Test that run_metagenome_pipeline works on tsv input seqtab.'''

        with TemporaryDirectory() as temp_dir:
            strat_out, unstrat_out = run_metagenome_pipeline(input_seqabun=seqtab_tsv_simple,
                                                             function=func_simple_in,
                                                             marker=marker_simple_in,
                                                             max_nsti=1.9,
                                                             out_dir=temp_dir,
                                                             strat_out=True,
                                                             wide_table=False)

        pd.testing.assert_frame_equal(strat_out.reset_index(drop=True),
                                      exp_strat_simple_in.reset_index(drop=True),
                                      check_like=True, atol=1e-3)
        pd.testing.assert_frame_equal(unstrat_out, exp_unstrat_simple_in,
                                      check_like=True, atol=1e-3)

    def test_full_pipeline_strat_tsv_skip_norm(self):
        '''Test that run_metagenome_pipeline works on tsv input seqtab and skip
        normalization.'''

        with TemporaryDirectory() as temp_dir:
            strat_out, unstrat_out = run_metagenome_pipeline(input_seqabun=seqtab_tsv_simple,
                                                             function=func_simple_in,
                                                             max_nsti=1.9,
                                                             out_dir=temp_dir,
                                                             strat_out=True,
                                                             wide_table=False,
                                                             skip_norm=True)

        pd.testing.assert_frame_equal(strat_out.reset_index(drop=True),
                                      exp_strat_simple_in.reset_index(drop=True),
                                      check_like=True, atol=1e-3)
        pd.testing.assert_frame_equal(unstrat_out, exp_unstrat_simple_in,
                                      check_like=True, atol=1e-3)

    def test_full_pipeline_strat_wide_tsv(self):
        '''Test that run_metagenome_pipeline works on tsv input seqtab. Compare
        with wide-format table in this case.'''

        with TemporaryDirectory() as temp_dir:
            strat_out, unstrat_out = run_metagenome_pipeline(input_seqabun=seqtab_tsv,
                                                             function=func_predict,
                                                             marker=marker_predict,
                                                             max_nsti=1.9,
                                                             out_dir=temp_dir,
                                                             strat_out=True,
                                                             wide_table=True)

        pd.testing.assert_frame_equal(strat_out, exp_strat_wide_in,
                                      check_like=True)

        pd.testing.assert_frame_equal(unstrat_out, exp_unstrat_in,
                                      check_like=True)

    def test_full_pipeline_unstrat_tsv_when_no_strat(self):
        '''Test that run_metagenome_pipeline works on tsv input seqtab when
        strat_out=False.'''

        with TemporaryDirectory() as temp_dir:
            strat_out, unstrat_out = run_metagenome_pipeline(input_seqabun=seqtab_tsv,
                                                             function=func_predict,
                                                             marker=marker_predict,
                                                             max_nsti=2.1,
                                                             out_dir=temp_dir,
                                                             strat_out=False)

        pd.testing.assert_frame_equal(unstrat_out, exp_unstrat_in,
                                      check_like=True)

    def test_full_pipeline_unstrat_msf_when_no_strat(self):
        '''Test that run_metagenome_pipeline works on mothur shared file input
        seqtab when strat_out=False.'''

        with TemporaryDirectory() as temp_dir:
            strat_out, unstrat_out = run_metagenome_pipeline(input_seqabun=seqtab_msf,
                                                             function=func_predict,
                                                             marker=marker_predict,
                                                             max_nsti=2.1,
                                                             out_dir=temp_dir,
                                                             strat_out=False)

        pd.testing.assert_frame_equal(unstrat_out, exp_unstrat_in,
                                      check_like=True)

    def test_full_pipeline_strat_wide_biom(self):
        '''Test that run_metagenome_pipeline creates correct stratified output
        on biom input seqtab. Compare with wide-format table in this case.'''

        with TemporaryDirectory() as temp_dir:
            strat_out, unstrat_out = run_metagenome_pipeline(input_seqabun=seqtab_biom,
                                                             function=func_predict,
                                                             marker=marker_predict,
                                                             max_nsti=2.0,
                                                             out_dir=temp_dir,
                                                             strat_out=True,
                                                             wide_table=True)

        pd.testing.assert_frame_equal(strat_out, exp_strat_wide_in,
                                      check_like=True)

    def test_full_pipeline_unstrat_biom_when_no_strat(self):
        '''Test that run_metagenome_pipeline create corrected unstratified
        output on biom input seqtab when strat_out=False.'''

        with TemporaryDirectory() as temp_dir:
            strat_out, unstrat_out = run_metagenome_pipeline(input_seqabun=seqtab_biom,
                                                             function=func_predict,
                                                             marker=marker_predict,
                                                             max_nsti=1.8,
                                                             out_dir=temp_dir,
                                                             strat_out=False)

        pd.testing.assert_frame_equal(unstrat_out, exp_unstrat_in,
                                      check_like=True)

    def test_full_pipeline_unstrat_biom(self):
        '''Test that run_metagenome_pipeline create corrected unstratified
        output on biom input seqtab.'''

        with TemporaryDirectory() as temp_dir:
            strat_out, unstrat_out = run_metagenome_pipeline(input_seqabun=seqtab_biom,
                                                             function=func_predict,
                                                             marker=marker_predict,
                                                             max_nsti=2.1,
                                                             out_dir=temp_dir,
                                                             strat_out=False)

        pd.testing.assert_frame_equal(unstrat_out, exp_unstrat_in,
                                      check_like=True)

    def test_norm_by_marker_copies(self):
        '''Test that expected normalized sequence abundance table generated.'''

        seqtab_in = biom.load_table(seqtab_biom).to_dataframe(dense=True)

        # Get output index labels in same order as expected.
        seqtab_in = seqtab_in.reindex(exp_norm_in.index)

        test_norm = norm_by_marker_copies(input_seq_counts=seqtab_in,
                                          input_marker_num=marker_predict_in,
                                          norm_filename=None)

        # Test whether normalized table matches expected table.
        pd.testing.assert_frame_equal(test_norm, exp_norm_in, check_like=True)

    def test_weighted_nsti(self):
        '''Test that expected weighted NSTI values are calculated.'''

        weighted_nsti_out = calc_weighted_nsti(exp_norm_in, nsti_in,
                                               return_df=True)

        expected_weighted_nsti = {"samples": ["sample1", "sample2", "sample3"],
                                  "weighted_NSTI": [0.292857143, 0.348275862,
                                                    0.436111111]}

        expected_weighted_nsti_df = pd.DataFrame.from_dict(expected_weighted_nsti)

        expected_weighted_nsti_df.set_index("samples", inplace=True)

        expected_weighted_nsti_df.index.name = None

        pd.testing.assert_frame_equal(weighted_nsti_out,
                                      expected_weighted_nsti_df,
                                      check_like=True)

    def test_nsti_filtering(self):
        '''Test that NSTI values filtered out correctly by checking for
        expected sequences to be retained..'''

        test_file = path.join(path.dirname(path.abspath(__file__)),
                              "test_data", "hsp", "hsp_output",
                              "mp_pred_out_nsti.tsv.gz")

        pred_test_in = pd.read_csv(test_file, sep="\t", index_col="sequence")

        pred_test_in_filt, nsti_col = drop_tips_by_nsti(tab=pred_test_in,
                                                     nsti_col="metadata_NSTI",
                                                        max_nsti=0.003)

        expected_passing_seqs = set(['2571042249_cluster', '2593339006',
                                 '2571042244', '2571042654', '2568526369',
                                 '2574180429_cluster', '2593338844',
                                 '2568526487_cluster', '2574180282_cluster'])

        self.assertSetEqual(expected_passing_seqs,
                            set(list(pred_test_in_filt.index)))

    def test_nsti_filtering_all_err(self):
        '''Test that error thrown when all ASVs thrown out.'''

        test_file = path.join(path.dirname(path.abspath(__file__)),
                              "test_data", "hsp", "hsp_output",
                              "mp_pred_out_nsti.tsv.gz")

        pred_test_in = pd.read_csv(test_file, sep="\t", index_col="sequence")

        with self.assertRaises(SystemExit):
            pred_test_in_filt, nsti_col = drop_tips_by_nsti(tab=pred_test_in,
                                                      nsti_col="metadata_NSTI",
                                                            max_nsti=0.000001)


class rare_seqs_test(unittest.TestCase):
    '''Checks that \"RARE\" category is being collapsed to correctly.'''

    def test_rare_4_reads(self):
        '''Check that correct sequences are identified as rare when a cut-off
        of 4 reads is used.'''

        seqtab_in = biom.load_table(seqtab_biom).to_dataframe(dense=True)

        rare_seqs = id_rare_seqs(seqtab_in, 4, 1)

        self.assertSetEqual(set(rare_seqs), set(["2558860574", "extra"]))

    def test_rare_2_samp(self):
        '''Check that correct sequences are identified as rare when a cut-off
        of 2 samples is used.'''

        seqtab_in = biom.load_table(seqtab_biom).to_dataframe(dense=True)

        rare_seqs = id_rare_seqs(seqtab_in, 1, 2)

        self.assertSetEqual(set(rare_seqs), set(["2558860574", "2571042244"]))

    def test_full_pipeline_strat_rare_category_tsv(self):
        '''Test that run_metagenome_pipeline works on tsv input seqtab and when
        rare seqs are collapsed into RARE category'''

        with TemporaryDirectory() as temp_dir:
            strat_out, unstrat_out = run_metagenome_pipeline(input_seqabun=seqtab_tsv_simple,
                                                             function=func_simple_in,
                                                             marker=marker_simple_in,
                                                             max_nsti=2.1,
                                                             min_reads=10,
                                                             min_samples=2,
                                                             out_dir=temp_dir,
                                                             strat_out=True,
                                                             wide_table=False)

        pd.testing.assert_frame_equal(strat_out.reset_index(drop=True),
                                      exp_strat_simple_rare_in.reset_index(drop=True),
                                      check_like=True)

        pd.testing.assert_frame_equal(unstrat_out, exp_unstrat_simple_in,
                                      check_like=True)


if __name__ == '__main__':
    unittest.main()
