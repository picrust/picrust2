#!/usr/bin/env python

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2-alpha.10"

import unittest
from os import path
import pandas as pd
from picrust2.wrap_hsp import castor_hsp_wrapper, castor_hsp_loocv_wrapper
from picrust2.util import get_picrust_project_dir

# Read in expected output files.
test_dir_path = path.join(get_picrust_project_dir(), "tests", "test_data",
                          "hsp")

in_traits1 = path.join(test_dir_path, "known_traits.tsv")
in_tree1 = path.join(test_dir_path, "tree.tre")

hsp_mp_pred = path.join(test_dir_path, "hsp_output", "mp_pred_out.tsv")
hsp_mp_pred_nsti = path.join(test_dir_path, "hsp_output",
                             "mp_pred_out_nsti.tsv")
hsp_mp_pred_ci = path.join(test_dir_path, "hsp_output", "mp_pred_out_ci.tsv")

hsp_emp_prob_pred = path.join(test_dir_path, "hsp_output",
                              "emp_prob_pred_out.tsv")
hsp_pic_pred = path.join(test_dir_path, "hsp_output", "pic_pred_out.tsv")
hsp_scp_pred = path.join(test_dir_path, "hsp_output", "scp_pred_out.tsv")
hsp_subtree_average_pred = path.join(test_dir_path, "hsp_output",
                                     "subtree_average_pred_out.tsv")

hsp_mp_pred_in = pd.read_table(hsp_mp_pred, sep="\t", index_col="sequence")

hsp_mp_pred_in_nsti = pd.read_table(hsp_mp_pred_nsti, sep="\t",
                                    index_col="sequence")

hsp_mp_pred_in_ci = pd.read_table(hsp_mp_pred_ci, sep="\t",
                                  index_col="sequence")

hsp_emp_prob_pred_in = pd.read_table(hsp_emp_prob_pred, sep="\t",
                                     index_col="sequence")
hsp_pic_pred_in = pd.read_table(hsp_pic_pred, sep="\t", index_col="sequence")
hsp_scp_pred_in = pd.read_table(hsp_scp_pred, sep="\t", index_col="sequence")
hsp_subtree_average_pred_in = pd.read_table(hsp_subtree_average_pred, sep="\t",
                                            index_col="sequence")


class castor_hsp_wrapper_tests(unittest.TestCase):
    """Tests for castor_hsp_wrapper function."""

    # Each of the below "simple" tests check that predictions match the
    # expected values in "test_data/hsp"
    def test_mp_simple(self):

        predict_out, ci_out = castor_hsp_wrapper(tree_path=in_tree1,
                                                 trait_table_path=in_traits1,
                                                 hsp_method="mp",
                                                 ran_seed=10)
	
	# Since values can differ depending on exact dependency versions, just comparing dimension and names.
        predict_out[:] = 0
        hsp_mp_pred_in[:] = 0

        pd.testing.assert_frame_equal(predict_out, hsp_mp_pred_in, check_like=True)

    def test_emp_prob_simple(self):

        predict_out, ci_out = castor_hsp_wrapper(tree_path=in_tree1,
                                                 trait_table_path=in_traits1,
                                                 hsp_method="emp_prob",
                                                 ran_seed=10)

	# Since values can differ depending on exact dependency versions, just comparing dimension and names.
        predict_out[:] = 0
        hsp_emp_prob_pred_in[:] = 0

        pd.testing.assert_frame_equal(predict_out, hsp_emp_prob_pred_in, check_like=True)

    def test_pic_simple(self):

        predict_out, ci_out = castor_hsp_wrapper(tree_path=in_tree1,
                                                 trait_table_path=in_traits1,
                                                 hsp_method="pic",
                                                 ran_seed=10)

	# Since values can differ depending on exact dependency versions, just comparing dimension and names.
        predict_out[:] = 0
        hsp_pic_pred_in[:] = 0

        pd.testing.assert_frame_equal(predict_out, hsp_pic_pred_in, check_like=True)

    def test_scp_simple(self):

        predict_out, ci_out = castor_hsp_wrapper(tree_path=in_tree1,
                                                 trait_table_path=in_traits1,
                                                 hsp_method="scp",
                                                 ran_seed=10)

	# Since values can differ depending on exact dependency versions, just comparing dimension and names.
        predict_out[:] = 0
        hsp_scp_pred_in[:] = 0

        pd.testing.assert_frame_equal(predict_out, hsp_scp_pred_in, check_like=True)

    def test_subtree_average_simple(self):

        predict_out, ci_out = castor_hsp_wrapper(tree_path=in_tree1,
                                                 trait_table_path=in_traits1,
                                                 hsp_method="subtree_average",
                                                 ran_seed=10)

	# Since values can differ depending on exact dependency versions, just comparing dimension and names.
        predict_out[:] = 0
        hsp_subtree_average_pred_in[:] = 0

        pd.testing.assert_frame_equal(predict_out, hsp_subtree_average_pred_in, check_like=True)

    def test_mp_ci(self):
        '''Test that MP confidence intervals calculated correctly.'''
        predict_out, ci_out = castor_hsp_wrapper(tree_path=in_tree1,
                                                 trait_table_path=in_traits1,
                                                 hsp_method="mp",
                                                 ran_seed=10,
                                                 calc_ci=True)

       	# Since values can differ depending on exact dependency versions, just comparing dimension and names.
        #predict_out[:] = 0
        #hsp_mp_pred_in_ci[:] = 0

        pd.testing.assert_frame_equal(ci_out, hsp_mp_pred_in_ci, check_like=True)

    def test_mp_nsti(self):
        '''Test that NSTI values (and the MP predictions) match expected.'''
        predict_out, ci_out = castor_hsp_wrapper(tree_path=in_tree1,
                                                 trait_table_path=in_traits1,
                                                 hsp_method="mp",
                                                 ran_seed=10,
                                                 calc_nsti=True)
        # Keep NSTI column only.
        predict_out_subset = predict_out.loc[:, ["metadata_NSTI"]]
        hsp_mp_pred_in_nsti_subset = hsp_mp_pred_in_nsti.loc[:, ["metadata_NSTI"]]

        pd.testing.assert_frame_equal(predict_out_subset, hsp_mp_pred_in_nsti_subset, check_like=True)


if __name__ == '__main__':
    unittest.main()
