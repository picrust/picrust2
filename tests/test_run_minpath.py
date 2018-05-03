#!/usr/bin/env python

__author__ = "Gavin Douglas"
__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2-alpha.8"

from biom import load_table
import unittest
import pandas as pd
from os import path
from tempfile import TemporaryDirectory
from picrust2.run_minpath import minpath_wrapper, run_minpath_pipeline
from picrust2.util import get_picrust_project_dir

# Path to test directory.
test_dir_path = path.join(get_picrust_project_dir(), "tests")

in_gene_abun = path.join(test_dir_path, "test_data", "run_minpath",
                         "test.genefamilies.biom")

in_gene_abun_in = load_table(in_gene_abun)
in_gene_abun_in.remove_empty(axis='whole', inplace=True)

functions = in_gene_abun_in.ids(axis="observation")

exp_minpath_out = path.join(test_dir_path, "test_data", "run_minpath",
                            "expected_minpath_out.tsv")

map_ec2path_prokaryotic = path.join(get_picrust_project_dir(), "MinPath",
                                    "ec2metacyc_picrust_prokaryotic.txt")


class minpath_wrapper_tests(unittest.TestCase):
    """Tests for minpath_wrapper function."""

    def test_minpath_wrapper_single_sample(self):
        '''Test running minpath_wrapper on single sample.'''

        with TemporaryDirectory() as temp_dir:
            path_abun = minpath_wrapper(sample_id="B1",
                                        biom_in=in_gene_abun_in,
                                        minpath_map=map_ec2path_prokaryotic,
                                        out_dir=temp_dir,
                                        functions=functions)

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
        exp_path_abun_B1 = exp_path_abun_B1.loc[~(exp_path_abun_B1 == 0).all(axis=1)]

        pd.testing.assert_frame_equal(exp_path_abun_B1, path_abun_df)


class run_minpath_pipeline_tests(unittest.TestCase):
    """Tests for run_minpath_pipeline function."""

    def test_basic_pipeline_2_threads(self):
        '''Test running full pipeline on 2 threads.'''

        with TemporaryDirectory() as temp_dir:
            test_metacyc_out = run_minpath_pipeline(inputfile=in_gene_abun,
                                                    mapfile=map_ec2path_prokaryotic,
                                                    threads=2,
                                                    out_dir=temp_dir)

        # Compare to expected pathway abundances.
        exp_path_abun = pd.read_csv(exp_minpath_out, sep="\t",
                                    index_col="pathway")

        test_metacyc_out.index.name = "pathway"

        pd.testing.assert_frame_equal(exp_path_abun, test_metacyc_out)


if __name__ == '__main__':
    unittest.main()
