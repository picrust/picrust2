#!/usr/bin/env python

__author__ = "Gavin Douglas"
__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2-alpha.7"

from biom import load_table
import sys, os, unittest
import pandas as pd
from picrust2.run_minpath import minpath_wrapper, run_minpath_pipeline
from picrust2.util import (generate_temp_filename, get_picrust_project_dir,
                           make_tmp_directory, system_call_check)

# Path to test directory.
test_dir_path = get_picrust_project_dir() + "/tests"

in_gene_abun = test_dir_path + "/" + \
               "test_data/run_minpath/test.genefamilies.biom"

in_gene_abun_in = load_table(in_gene_abun)
in_gene_abun_in.remove_empty(axis='whole', inplace=True)

functions = in_gene_abun_in.ids(axis="observation")

exp_minpath_out = test_dir_path + "/" + \
                  "test_data/run_minpath/expected_minpath_out.tsv"

map_ec2path_prokaryotic = get_picrust_project_dir() + "/MinPath/" + \
                          "ec2metacyc_picrust_prokaryotic.txt"


class minpath_wrapper_tests(unittest.TestCase):
    """Tests for minpath_wrapper function."""

    def test_minpath_wrapper_single_sample(self):
        '''Test running minpath_wrapper on single sample.'''

        single_tmp_dir = make_tmp_directory(dir_prefix=test_dir_path + "/tmp/") 

        path_abun = minpath_wrapper(sample_id="B1",
                                    biom_in=in_gene_abun_in,
                                    minpath_map=map_ec2path_prokaryotic,
                                    tmp_dir=single_tmp_dir,
                                    functions=functions)

        # Delete tmp directory.
        system_call_check("rm -r " + test_dir_path + "/tmp/")

        # Convert to pandas dataframe.
        path_abun_df = pd.DataFrame([path_abun])

        path_abun_df = path_abun_df.transpose()

        #  Set index title and column name.
        path_abun_df.index.name = "pathway"
        path_abun_df.columns = ["B1"]

        # Compare this predicted column to expected (after removing rows that
        # are 0).
        exp_path_abun = pd.read_csv(exp_minpath_out,
                                    sep="\t",
                                    index_col="pathway")

        exp_path_abun_B1 = exp_path_abun[["B1"]]

        # Remove rows that are all 0s.
        exp_path_abun_B1 = exp_path_abun_B1.loc[~(exp_path_abun_B1==0).all(axis=1)]

        pd.testing.assert_frame_equal(exp_path_abun_B1, path_abun_df)


class run_minpath_pipeline_tests(unittest.TestCase):
    """Tests for run_minpath_pipeline function."""

    def test_basic_pipeline_2_threads(self):
        '''Test running full pipeline on 2 threads.'''

        test_metacyc_out = run_minpath_pipeline(inputfile=in_gene_abun,
                                                mapfile=map_ec2path_prokaryotic,
                                                keep_tmp=False,
                                                threads=2,
                                                tmp_dir=test_dir_path + "/tmp/")

        # Compare to expected pathway abundances.
        exp_path_abun = pd.read_csv(exp_minpath_out, sep="\t",
                                    index_col="pathway")

        test_metacyc_out.index.name = "pathway"

        pd.testing.assert_frame_equal(exp_path_abun, test_metacyc_out)


if __name__ == '__main__':
    unittest.main()
