#!/usr/bin/env python

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2-alpha.9"

import unittest
from os import path
import pandas as pd
from tempfile import TemporaryDirectory
from picrust2.util import (write_fasta, read_fasta, write_phylip, read_phylip,
                           three_df_index_overlap_sort)

# Inititalize 3 test pandas dataframes.
test1 = pd.DataFrame.from_dict({"a": [1, 2, 3],
                                "b": [10, 6, 5],
                                "c": [1, 3, 4]},
                               orient='index')

test2 = pd.DataFrame.from_dict({"b": [10, 7, 0],
                                "c": [0, 0, 0],
                                "d": [5, 4, 3]},
                               orient='index')

test3 = pd.DataFrame.from_dict({"e": [0, 5, 3],
                                "c": [9, 0, 7],
                                "b": [1, 0, 2]},
                               orient='index')


class util_test(unittest.TestCase):

    def test_read_write_fasta(self):
        '''Basic test that FASTA files are read and written correctly.'''

        test_seqs_dict = {"seq1": "GNATNGAC",
                          "seq2": "GTCGTGGC",
                          "seq3": "GNCTGAGATTAACC"}

        # Write these sequences temp file and then read them back in again.
        with TemporaryDirectory() as temp_dir:
            outfile = path.join(temp_dir, "test.fna")

            write_fasta(test_seqs_dict, outfile)

            test_seqs_dict_in = read_fasta(outfile)

        self.assertEqual(test_seqs_dict, test_seqs_dict_in)

    def test_read_write_phylip(self):
        '''Basic test that Phylip files are read and written correctly.'''

        test_seqs_dict = {"seq1": "GNATNGAC",
                          "seq2": "GTCGTGGC",
                          "seq3": "GNCTGAGA"}

        # Write these sequences temp file and then read them back in again.
        with TemporaryDirectory() as temp_dir:
            outfile = path.join(temp_dir, "test.phylip")

            write_phylip(test_seqs_dict, outfile)

            test_seqs_dict_in = read_phylip(outfile)

        self.assertEqual(test_seqs_dict, test_seqs_dict_in)

    def test_three_df_index_overlap_sort(self):
        '''Test that function to subset and sort dataframes to overlapping
        index labels is generating expected output.'''

        # Get subset of overlapped index labels ordered for each dataframe.
        test1_out, test2_out, test3_out = three_df_index_overlap_sort(test1,
                                                                      test2,
                                                                      test3)

        # Get expected output with different approach as well.
        exp_out = [test1.loc[["b", "c"], :], test2.loc[["b", "c"], :],
                   test3.loc[["b", "c"], :]]

        # Check that output dataframes match the expected objects.
        pd.testing.assert_frame_equal(test1_out, exp_out[0])
        pd.testing.assert_frame_equal(test2_out, exp_out[1])
        pd.testing.assert_frame_equal(test3_out, exp_out[2])

    def test_three_df_index_overlap_sort_err(self):
        '''Test that function to subset and sort dataframes to overlapping
        index labels will give error if index label don't overlap..'''

        # Define df without overlapping index labels.
        test_no_overlap = pd.DataFrame.from_dict({"e": [0, 5, 3],
                                                  "f": [9, 0, 7],
                                                  "g": [1, 0, 2]},
                                                 orient='index')

        # Check that ValueError assertion raised.
        self.assertRaises(ValueError, three_df_index_overlap_sort, test1,
                          test2, test_no_overlap)


if __name__ == '__main__':
    unittest.main()
