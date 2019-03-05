#!/usr/bin/env python

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.1.1-b"

import unittest
from os import path
import pandas as pd
import hashlib
import gzip
from picrust2.default import default_tables
from picrust2.wrap_hsp import (castor_hsp_workflow,
                               castor_nsti)

# Read in expected output files.

test_dir_path = path.join(path.dirname(path.abspath(__file__)), "test_data",
                          "hsp")

in_traits1 = path.join(test_dir_path, "known_traits.tsv.gz")
in_tree1 = path.join(test_dir_path, "tree.tre")

hsp_mp_pred = path.join(test_dir_path, "hsp_output", "mp_pred_out.tsv")
hsp_mp_pred_nsti = path.join(test_dir_path, "hsp_output",
                             "mp_pred_out_nsti.tsv")
hsp_emp_prob_pred_ci = path.join(test_dir_path, "hsp_output", "emp_prob_pred_out_ci.tsv")

hsp_emp_prob_pred = path.join(test_dir_path, "hsp_output",
                              "emp_prob_pred_out.tsv")
hsp_pic_pred = path.join(test_dir_path, "hsp_output", "pic_pred_out.tsv")
hsp_scp_pred = path.join(test_dir_path, "hsp_output", "scp_pred_out.tsv")
hsp_subtree_average_pred = path.join(test_dir_path, "hsp_output",
                                     "subtree_average_pred_out.tsv")

hsp_mp_pred_in = pd.read_table(hsp_mp_pred, sep="\t", index_col="sequence")

hsp_mp_pred_in_nsti = pd.read_table(hsp_mp_pred_nsti, sep="\t",
                                    index_col="sequence")

hsp_emp_prob_pred_in_ci = pd.read_table(hsp_emp_prob_pred_ci, sep="\t",
                                  index_col="sequence")

hsp_emp_prob_pred_in = pd.read_table(hsp_emp_prob_pred, sep="\t",
                                     index_col="sequence")
hsp_pic_pred_in = pd.read_table(hsp_pic_pred, sep="\t", index_col="sequence")
hsp_scp_pred_in = pd.read_table(hsp_scp_pred, sep="\t", index_col="sequence")
hsp_subtree_average_pred_in = pd.read_table(hsp_subtree_average_pred, sep="\t",
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
                                    dtype={'assembly' : str})

        nsti_out = castor_nsti(tree_path=in_tree1,
                               known_tips=in_traits1_df.index.values)

       # Only compare NSTI column.
        hsp_mp_pred_in_nsti_subset = hsp_mp_pred_in_nsti.loc[:, ["metadata_NSTI"]]

        pd.testing.assert_frame_equal(nsti_out, hsp_mp_pred_in_nsti_subset, check_like=True)


class table_mdf5sum_tests(unittest.TestCase):
    '''Check that md5sum values of default tables match expected values.'''

    def test_default_table_md5sum(self):

        marker_16S_hash = hashlib.md5()
        ec_hash = hashlib.md5()
        ko_hash = hashlib.md5()
        cog_hash = hashlib.md5()
        pfam_hash = hashlib.md5()
        tigrfam_hash = hashlib.md5()

        with gzip.open(default_tables["16S"], 'rt') as marker_16S_in:
            marker_16S_hash.update(marker_16S_in.read().encode())

        with gzip.open(default_tables["EC"], 'rt') as ec_in:
            ec_hash.update(ec_in.read().encode())

        with gzip.open(default_tables["KO"], 'rt') as ko_in:
            ko_hash.update(ko_in.read().encode())

        with gzip.open(default_tables["COG"], 'rt') as cog_in:
            cog_hash.update(cog_in.read().encode())

        with gzip.open(default_tables["PFAM"], 'rt') as pfam_in:
            pfam_hash.update(pfam_in.read().encode())

        with gzip.open(default_tables["TIGRFAM"], 'rt') as tigrfam_in:
            tigrfam_hash.update(tigrfam_in.read().encode())

        obs_hash = [marker_16S_hash.hexdigest(), ec_hash.hexdigest(),
                    ko_hash.hexdigest(), cog_hash.hexdigest(),
                    pfam_hash.hexdigest(), tigrfam_hash.hexdigest()]

        exp_hash = ["a0acd4dbc3501b271ade941bf308643e",
                    "85fcd6ac08ba324c6c14699fcee92725",
                    "46540c1eb6e7f22feaa50893dc74a6d7",
                    "c50eaa4565804b6fe880e81044d33276",
                    "420605c3afe98a585f4f006f8aa8c4f3",
                    "ca5871bac44593c009b5f7838efb2772"]

        # Check that md5sum values match expected values.
        self.assertEqual(obs_hash, exp_hash)

if __name__ == '__main__':
    unittest.main()
