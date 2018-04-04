#!/usr/bin/env python

__author__ = "Gavin Douglas"
__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2-alpha.6"

import sys, os, unittest
from picrust.wrap_hsp import castor_hsp_wrapper, castor_hsp_loocv_wrapper
from cogent import LoadTable

# Path to test directory.
test_dir_path = os.path.dirname(sys.argv[0])

in_traits1 = test_dir_path + "/" + "test_data/hsp/known_traits.tsv"
in_tree1 = test_dir_path + "/" + "test_data/hsp/tree.tre"

hsp_mp_pred = test_dir_path + "/" + "test_data/hsp/hsp_output/mp_pred_out.tsv"
hsp_emp_prob_pred = test_dir_path + "/" + "test_data/hsp/hsp_output/emp_prob_pred_out.tsv"
hsp_pic_pred = test_dir_path + "/" + "test_data/hsp/hsp_output/pic_pred_out.tsv"
hsp_scp_pred = test_dir_path + "/" + "test_data/hsp/hsp_output/scp_pred_out.tsv"
hsp_subtree_average_pred = test_dir_path + "/" + "test_data/hsp/hsp_output/subtree_average_pred_out.tsv"

hsp_mp_pred_in = LoadTable(hsp_mp_pred, header=True, sep='\t')
hsp_emp_prob_pred_in = LoadTable(hsp_emp_prob_pred, header=True, sep='\t')
hsp_pic_pred_in = LoadTable(hsp_pic_pred, header=True, sep='\t')
hsp_scp_pred_in = LoadTable(hsp_scp_pred, header=True, sep='\t')
hsp_subtree_average_pred_in = LoadTable(hsp_subtree_average_pred, header=True, sep='\t')


class castor_hsp_wrapper_tests(unittest.TestCase):
    """Tests for castor_hsp_wrapper function."""

    # Each of the below "simple" tests check that predictions match the
    # expected values in "test_data/hsp"
    def test_mp_simple(self):

        predict_out, ci_out = castor_hsp_wrapper(tree_path=in_tree1,
                                                 trait_table_path=in_traits1,
                                                 hsp_method="mp",
                                                 ran_seed=10)

        self.assertEqual(predict_out, hsp_mp_pred_in)


    def test_emp_prob_simple(self):

        predict_out, ci_out = castor_hsp_wrapper(tree_path=in_tree1,
                                                 trait_table_path=in_traits1,
                                                 hsp_method="emp_prob",
                                                 ran_seed=10)

        self.assertEqual(predict_out, hsp_emp_prob_pred_in)


    def test_pic_simple(self):

        predict_out, ci_out = castor_hsp_wrapper(tree_path=in_tree1,
                                                 trait_table_path=in_traits1,
                                                 hsp_method="pic",
                                                 ran_seed=10)

        self.assertEqual(predict_out, hsp_pic_pred_in)


    def test_scp_simple(self):

        predict_out, ci_out = castor_hsp_wrapper(tree_path=in_tree1,
                                                 trait_table_path=in_traits1,
                                                 hsp_method="scp",
                                                 ran_seed=10)

        self.assertEqual(predict_out, hsp_scp_pred_in)


    def test_subtree_average_simple(self):

        predict_out, ci_out = castor_hsp_wrapper(tree_path=in_tree1,
                                                 trait_table_path=in_traits1,
                                                 hsp_method="subtree_average",
                                                 ran_seed=10)

        self.assertEqual(predict_out, hsp_subtree_average_pred_in)


    # With mp method:
        # Check that can run with and without confidence intervals (+ check that values match exactly)
        # Check that can run with and without NSTI calculation (+ check that values match exactly)
        # Check that can run with and without --check option.
        # Check that can run on 1 or 2 processors.

    # Also try out a couple of different trees (e.g. rooted vs unrooted)
    # And a trait table with a lot of variation in values.


###class castor_hsp_wrapper_tests(unittest.TestCase):

if __name__ == '__main__':
    unittest.main()
