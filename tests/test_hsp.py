#!/usr/bin/env python

import unittest
from os import path
import pandas as pd
import hashlib
import gzip
from picrust2.default import default_tables_bac, default_tables_arc
from picrust2.wrap_hsp import (castor_hsp_workflow,
                               castor_nsti)

# Read in expected output files.
test_dir_path = path.join(path.dirname(path.abspath(__file__)), "test_data",
                          "hsp")

in_traits1 = path.join(test_dir_path, "known_traits.tsv.gz")
in_tree1 = path.join(test_dir_path, "tree.tre")

hsp_mp_pred = path.join(test_dir_path, "hsp_output", "mp_pred_out.tsv.gz")
hsp_mp_pred_nsti = path.join(test_dir_path, "hsp_output",
                             "mp_pred_out_nsti.tsv.gz")
hsp_emp_prob_pred_ci = path.join(test_dir_path, "hsp_output", "emp_prob_pred_out_ci.tsv.gz")

hsp_emp_prob_pred = path.join(test_dir_path, "hsp_output",
                              "emp_prob_pred_out.tsv.gz")
hsp_pic_pred = path.join(test_dir_path, "hsp_output", "pic_pred_out.tsv.gz")
hsp_scp_pred = path.join(test_dir_path, "hsp_output", "scp_pred_out.tsv.gz")
hsp_subtree_average_pred = path.join(test_dir_path, "hsp_output",
                                     "subtree_average_pred_out.tsv.gz")

hsp_mp_pred_in = pd.read_csv(hsp_mp_pred, sep="\t", index_col="sequence")

hsp_mp_pred_in_nsti = pd.read_csv(hsp_mp_pred_nsti, sep="\t",
                                  index_col="sequence")

hsp_emp_prob_pred_in_ci = pd.read_csv(hsp_emp_prob_pred_ci, sep="\t",
                                      index_col="sequence")

hsp_emp_prob_pred_in = pd.read_csv(hsp_emp_prob_pred, sep="\t",
                                   index_col="sequence")
hsp_pic_pred_in = pd.read_csv(hsp_pic_pred, sep="\t", index_col="sequence")
hsp_scp_pred_in = pd.read_csv(hsp_scp_pred, sep="\t", index_col="sequence")
hsp_subtree_average_pred_in = pd.read_csv(hsp_subtree_average_pred, sep="\t",
                                          index_col="sequence")


class castor_hsp_workflow_tests(unittest.TestCase):
    """Tests for castor hsp workflow."""

    # Each of the below "simple" tests check that predictions match the
    # expected values in "test_data/hsp"
    def test_mp_simple(self):

        predict_out, ci_out = castor_hsp_workflow(tree_path=in_tree1,
                                                  trait_table_path=in_traits1,
                                                  hsp_method="mp",
                                                  ran_seed=10)

        pd.testing.assert_frame_equal(predict_out, hsp_mp_pred_in, check_like=True)

    def test_emp_prob_simple(self):

        predict_out, ci_out = castor_hsp_workflow(tree_path=in_tree1,
                                                  trait_table_path=in_traits1,
                                                  hsp_method="emp_prob",
                                                  ran_seed=10)

        # Since values can differ depending on exact dependency versions, just comparing dimension and names.
        predict_out[:] = 0
        hsp_emp_prob_pred_in[:] = 0

        pd.testing.assert_frame_equal(predict_out, hsp_emp_prob_pred_in, check_like=True)

    def test_pic_simple(self):

        predict_out, ci_out = castor_hsp_workflow(tree_path=in_tree1,
                                                  trait_table_path=in_traits1,
                                                  hsp_method="pic",
                                                  ran_seed=10)

        # Since values can differ depending on exact dependency versions, just comparing dimension and names.
        predict_out[:] = 0
        hsp_pic_pred_in[:] = 0

        pd.testing.assert_frame_equal(predict_out, hsp_pic_pred_in, check_like=True)

    def test_scp_simple(self):

        predict_out, ci_out = castor_hsp_workflow(tree_path=in_tree1,
                                                  trait_table_path=in_traits1,
                                                  hsp_method="scp",
                                                  ran_seed=10)

        # Since values can differ depending on exact dependency versions, just comparing dimension and names.
        predict_out[:] = 0
        hsp_scp_pred_in[:] = 0

        pd.testing.assert_frame_equal(predict_out, hsp_scp_pred_in, check_like=True)

    def test_subtree_average_simple(self):

        predict_out, ci_out = castor_hsp_workflow(tree_path=in_tree1,
                                                  trait_table_path=in_traits1,
                                                  hsp_method="subtree_average",
                                                  ran_seed=10)

        # Since values can differ depending on exact dependency versions, just comparing dimension and names.
        predict_out[:] = 0
        hsp_subtree_average_pred_in[:] = 0

        pd.testing.assert_frame_equal(predict_out, hsp_subtree_average_pred_in, check_like=True)

    def test_emp_prob_ci(self):
        '''Test that Emp Prob confidence intervals calculated correctly.'''
        predict_out, ci_out = castor_hsp_workflow(tree_path=in_tree1,
                                                  trait_table_path=in_traits1,
                                                  hsp_method="emp_prob",
                                                  ran_seed=10,
                                                  calc_ci=True)

        pd.testing.assert_frame_equal(ci_out, hsp_emp_prob_pred_in_ci, check_like=True)

    def test_nsti(self):
        '''Test that calculated NSTI values match expected.'''

        in_traits1_df = pd.read_csv(in_traits1, sep="\t",
                                    index_col="assembly",
                                    dtype={'assembly': str})

        nsti_out = castor_nsti(tree_path=in_tree1,
                               known_tips=in_traits1_df.index.values,
                               verbose=True)
        
        # Drop the best_domain column that was added for PICRUSt2-v2.6.0 so that the file still matches the old one
        nsti_out = nsti_out.drop('closest_reference_genome', axis=1)

        # Only compare NSTI column.
        hsp_mp_pred_in_nsti_subset = hsp_mp_pred_in_nsti.loc[:, ["metadata_NSTI"]]

        pd.testing.assert_frame_equal(nsti_out, hsp_mp_pred_in_nsti_subset, check_like=True)


class table_mdf5sum_tests(unittest.TestCase):
    '''Check that md5sum values of default tables match expected values.'''

    def test_default_table_md5sum_bac(self):

        marker_16S_hash = hashlib.md5()
        ec_hash = hashlib.md5()
        ko_hash = hashlib.md5()
        go_hash = hashlib.md5()
        pfam_hash = hashlib.md5()
        bigg_hash = hashlib.md5()
        cazy_hash = hashlib.md5()
        gene_names_hash = hashlib.md5()

        with gzip.open(default_tables_bac["16S"], 'rt') as marker_16S_in:
            marker_16S_hash.update(marker_16S_in.read().encode())

        with gzip.open(default_tables_bac["EC"], 'rt') as ec_in:
            ec_hash.update(ec_in.read().encode())

        with gzip.open(default_tables_bac["KO"], 'rt') as ko_in:
            ko_hash.update(ko_in.read().encode())

        with gzip.open(default_tables_bac["GO"], 'rt') as go_in:
            go_hash.update(go_in.read().encode())

        with gzip.open(default_tables_bac["PFAM"], 'rt') as pfam_in:
            pfam_hash.update(pfam_in.read().encode())

        with gzip.open(default_tables_bac["BIGG"], 'rt') as bigg_in:
            bigg_hash.update(bigg_in.read().encode())

        with gzip.open(default_tables_bac["CAZY"], 'rt') as cazy_in:
            cazy_hash.update(cazy_in.read().encode())

        with gzip.open(default_tables_bac["GENE_NAMES"], 'rt') as gene_names_in:
            gene_names_hash.update(gene_names_in.read().encode())        

        obs_hash = [marker_16S_hash.hexdigest(), ec_hash.hexdigest(),
                    ko_hash.hexdigest(), go_hash.hexdigest(),
                    pfam_hash.hexdigest(), bigg_hash.hexdigest(),
                    cazy_hash.hexdigest(), gene_names_hash.hexdigest()]

        exp_hash = ["3e870a33486bb7800dbfe497a2d6db11",
                    "96bc131933a661a192b8e0a1e0182254",
                    "af7363ea3d67e96077a37d7d8bdafd89",
                    "4ee72f88f84a3653d1d15abb27cb42b7",
                    "a892738bd3e0cd1454c26f2a78decf6f",
                    "b7e85755cf1016c3e55e0dbcdef56c39",
                    "7f0520b27de9bb85cc7a1c3d37a4fec5",
                    "17e8e23bb1f5b73ecfed2cb721e26449"]

        # Check that md5sum values match expected values.
        self.assertEqual(obs_hash, exp_hash)

    def test_default_table_md5sum_arc(self):

        marker_16S_hash = hashlib.md5()
        ec_hash = hashlib.md5()
        ko_hash = hashlib.md5()
        go_hash = hashlib.md5()
        pfam_hash = hashlib.md5()
        bigg_hash = hashlib.md5()
        cazy_hash = hashlib.md5()
        gene_names_hash = hashlib.md5()

        with gzip.open(default_tables_arc["16S"], 'rt') as marker_16S_in:
            marker_16S_hash.update(marker_16S_in.read().encode())

        with gzip.open(default_tables_arc["EC"], 'rt') as ec_in:
            ec_hash.update(ec_in.read().encode())

        with gzip.open(default_tables_arc["KO"], 'rt') as ko_in:
            ko_hash.update(ko_in.read().encode())

        with gzip.open(default_tables_arc["GO"], 'rt') as go_in:
            go_hash.update(go_in.read().encode())

        with gzip.open(default_tables_arc["PFAM"], 'rt') as pfam_in:
            pfam_hash.update(pfam_in.read().encode())

        with gzip.open(default_tables_arc["BIGG"], 'rt') as bigg_in:
            bigg_hash.update(bigg_in.read().encode())

        with gzip.open(default_tables_arc["CAZY"], 'rt') as cazy_in:
            cazy_hash.update(cazy_in.read().encode())

        with gzip.open(default_tables_arc["GENE_NAMES"], 'rt') as gene_names_in:
            gene_names_hash.update(gene_names_in.read().encode())

        obs_hash = [marker_16S_hash.hexdigest(), ec_hash.hexdigest(),
                    ko_hash.hexdigest(), go_hash.hexdigest(),
                    pfam_hash.hexdigest(), bigg_hash.hexdigest(),
                    cazy_hash.hexdigest(), gene_names_hash.hexdigest()]

        exp_hash = ["2a617a0c15d3aef4f849e12f38f6ef70",
                    "7b13c586c3cb989257e30abeb5271abb",
                    "f9e548b617672ff0028e0430e0c66e49",
                    "b48906876606ac4a89f532156afde80c",
                    "b3d5db008edc070e13f8736a1c25f10c",
                    "4523a5256d2db24d6fb0d9605b9cfc75",
                    "5d8208d7aedec81d054d429c735bf56d",
                    "8c5f8eb084c7cc028174e89744d2b0b3"]

        # Check that md5sum values match expected values.
        self.assertEqual(obs_hash, exp_hash)

class castor_hsp_chunk_tests(unittest.TestCase):
    '''Check all trait columns are predicted even with small chunk size.'''

    def test_all_traits_predicted_any_chunk_size(self):
        exp_cols = set(pd.read_csv(in_traits1, sep="\t", index_col="assembly").columns)

        for chunk_size in [1, 2, len(exp_cols), 500]:
            predict_out, _ = castor_hsp_workflow(tree_path=in_tree1,
                                                 trait_table_path=in_traits1,
                                                 hsp_method="mp",
                                                 chunk_size=chunk_size)

            self.assertEqual(set(predict_out.columns), exp_cols, msg="chunk_size=" + str(chunk_size))



if __name__ == '__main__':
    unittest.main()
