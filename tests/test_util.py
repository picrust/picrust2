#!/usr/bin/env python

__copyright__ = "Copyright 2018-2019, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.0.0-b.3"

import unittest
from os import path
import pandas as pd
import hashlib
import gzip
from picrust2.util import (write_fasta,
                           read_fasta,
                           write_phylip,
                           read_phylip,
                           three_df_index_overlap_sort,
                           add_descrip_col,
                           convert_humann2_to_picrust2,
                           convert_picrust2_to_humann2,
                           convert_picrust2_to_humann2_merged,
                           read_seqabun,
                           TemporaryDirectory)

from picrust2.default import default_map

metagenome_pipeline_test_dir_path = path.join(path.dirname(path.abspath(__file__)),
                                              "test_data",
                                              "metagenome_pipeline")
seqtab_biom = path.join(metagenome_pipeline_test_dir_path,
                       "test_input_sequence_abun.biom")
seqtab_msf = path.join(metagenome_pipeline_test_dir_path,
                       "test_input_sequence_abun.msf")

descrip_test_dir_path = path.join(path.dirname(path.abspath(__file__)),
                                  "test_data",
                                  "add_descriptions")

descrip_test_dir_out_path = path.join(descrip_test_dir_path, "output")

# Set paths to test input and output files for add_descriptions.py tests.
ec_unstrat_in = path.join(descrip_test_dir_path, "ec_unstrat_test.txt")
ec_unstrat_exp = path.join(descrip_test_dir_out_path, "ec_unstrat_exp.txt")

ec_strat_in = path.join(descrip_test_dir_path, "ec_strat_test.txt")
ec_strat_exp = path.join(descrip_test_dir_out_path, "ec_strat_exp.txt")

ec_nomatch_in = path.join(descrip_test_dir_path, "ec_nomatch_test.txt")

metacyc_unstrat_in = path.join(descrip_test_dir_path, "metacyc_unstrat_test.txt")
metacyc_unstrat_exp = path.join(descrip_test_dir_out_path, "metacyc_unstrat_exp.txt")

cog_unstrat_in = path.join(descrip_test_dir_path, "cog_unstrat_test.txt")
cog_unstrat_exp = path.join(descrip_test_dir_out_path, "cog_unstrat_exp.txt")

ko_unstrat_in = path.join(descrip_test_dir_path, "ko_unstrat_test.txt")
ko_unstrat_exp = path.join(descrip_test_dir_out_path, "ko_unstrat_exp.txt")

pfam_unstrat_in = path.join(descrip_test_dir_path, "pfam_unstrat_test.txt")
pfam_unstrat_exp = path.join(descrip_test_dir_out_path, "pfam_unstrat_exp.txt")

tigrfam_unstrat_in = path.join(descrip_test_dir_path, "tigrfam_unstrat_test.txt")
tigrfam_unstrat_exp = path.join(descrip_test_dir_out_path, "tigrfam_unstrat_exp.txt")

# Set paths to input and output files for convert_table.py tests.
convert_test_dir_path = path.join(path.dirname(path.abspath(__file__)), "test_data",
                                  "convert_table")

convert_test_dir_out_path = path.join(convert_test_dir_path, "expected_out")

humann2_strat_in = path.join(convert_test_dir_path, "humann2_strat_example.tsv")
humann2_strat_exp = path.join(convert_test_dir_out_path, "humann2_strat_example_picrust2.tsv")

humann2_unstrat_in = path.join(convert_test_dir_path, "humann2_unstrat_example.tsv")
humann2_unstrat_exp = path.join(convert_test_dir_out_path, "humann2_unstrat_example_picrust2.tsv")

humann2_strat_in_split1 = path.join(convert_test_dir_path, "humann2_path_test_A1.tsv")
humann2_strat_in_split2 = path.join(convert_test_dir_path, "humann2_path_test_D2.tsv")
humann2_strat_in_split3 = path.join(convert_test_dir_path, "humann2_path_test_PD3.tsv")
humann2_strat_split_exp = path.join(convert_test_dir_out_path, "humann2_path_test_picrust2.tsv")

picrust2_strat_in = path.join(convert_test_dir_path, "picrust2_strat_tmp.tsv")

picrust2_strat_exp1 = path.join(convert_test_dir_out_path, "picrust2_strat_to_humann2_split", "s1_split.tsv")
picrust2_strat_exp2 = path.join(convert_test_dir_out_path, "picrust2_strat_to_humann2_split", "s2_split.tsv")
picrust2_strat_exp3 = path.join(convert_test_dir_out_path, "picrust2_strat_to_humann2_split", "s3_split.tsv")

picrust2_unstrat_in1 = path.join(convert_test_dir_path, "picrust2_unstrat_tmp1.tsv")
picrust2_unstrat_in2 = path.join(convert_test_dir_path, "picrust2_unstrat_tmp2.tsv")
picrust2_unstrat_exp = path.join(convert_test_dir_out_path, "picrust2_unstrat_tmp1_tmp2_humann2-format.tsv")

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

    def test_seqabun_reading(self):
        '''Test that mothur shared format and BIOM tables read in
        identically.'''

        seqtab_biom_in = read_seqabun(seqtab_biom)

        seqtab_msf_in = read_seqabun(seqtab_msf)

        pd.testing.assert_frame_equal(seqtab_biom_in, seqtab_msf_in,
                                      check_dtype=False)

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


class add_description_tests(unittest.TestCase):

    def test_ec_unstrat(self):
        '''Test that correct descriptions added to table of EC abundances.'''

        obs_out = add_descrip_col(ec_unstrat_in, default_map["EC"])

        # Read in expected out.
        exp_out = pd.read_csv(ec_unstrat_exp, sep='\t')

        pd.testing.assert_frame_equal(obs_out, exp_out, check_like=True)


    def test_ec_strat(self):
        '''Test that correct descriptions added to table of stratified EC
        abundances.'''

        obs_out = add_descrip_col(ec_strat_in, default_map["EC"])

        # Read in expected out.
        exp_out = pd.read_csv(ec_strat_exp, sep='\t')

        pd.testing.assert_frame_equal(obs_out, exp_out, check_like=True)


    def test_nomatch_fail(self):
        '''Test that a system exit is made when the ids do not overlap.'''

        with self.assertRaises(SystemExit):
            add_descrip_col(ec_nomatch_in, default_map["EC"])

    def test_metacyc_unstrat(self):
        '''Test that correct descriptions added to table of metacyc pathway
        abundances.'''

        obs_out = add_descrip_col(metacyc_unstrat_in, default_map["METACYC"])

        # Read in expected out.
        exp_out = pd.read_csv(metacyc_unstrat_exp, sep='\t')

        pd.testing.assert_frame_equal(obs_out, exp_out, check_like=True)

    def test_cog_unstrat(self):
        '''Test that correct descriptions added to table of COG abundances.'''

        obs_out = add_descrip_col(cog_unstrat_in, default_map["COG"])

        # Read in expected out.
        exp_out = pd.read_csv(cog_unstrat_exp, sep='\t')

        pd.testing.assert_frame_equal(obs_out, exp_out, check_like=True)

    def test_ko_unstrat(self):
        '''Test that correct descriptions added to table of KO abundances.'''

        obs_out = add_descrip_col(ko_unstrat_in, default_map["KO"])

        # Read in expected out.
        exp_out = pd.read_csv(ko_unstrat_exp, sep='\t')

        pd.testing.assert_frame_equal(obs_out, exp_out, check_like=True)

    def test_pfam_unstrat(self):
        '''Test that correct descriptions added to table of PFAM abundances.'''

        obs_out = add_descrip_col(pfam_unstrat_in, default_map["PFAM"])

        # Read in expected out.
        exp_out = pd.read_csv(pfam_unstrat_exp, sep='\t')

        pd.testing.assert_frame_equal(obs_out, exp_out, check_like=True)

    def test_tigrfam_unstrat(self):
        '''Test that correct descriptions added to table of TIGRFAM
        abundances.'''

        obs_out = add_descrip_col(tigrfam_unstrat_in, default_map["TIGRFAM"])

        # Read in expected out.
        exp_out = pd.read_csv(tigrfam_unstrat_exp, sep='\t')

        pd.testing.assert_frame_equal(obs_out, exp_out, check_like=True)

    def test_default_mapping_md5sum(self):
        '''Check that md5sum values of default mapfiles match expected values.'''

        ec_hash = hashlib.md5()
        ko_hash = hashlib.md5()
        cog_hash = hashlib.md5()
        pfam_hash = hashlib.md5()
        tigrfam_hash = hashlib.md5()
        metacyc_hash = hashlib.md5()

        with gzip.open(default_map["EC"], 'rt') as ec_in:
            ec_hash.update(ec_in.read().encode())

        with gzip.open(default_map["KO"], 'rt') as ko_in:
            ko_hash.update(ko_in.read().encode())

        with gzip.open(default_map["COG"], 'rt') as cog_in:
            cog_hash.update(cog_in.read().encode())

        with gzip.open(default_map["PFAM"], 'rt') as pfam_in:
            pfam_hash.update(pfam_in.read().encode())

        with gzip.open(default_map["TIGRFAM"], 'rt') as tigrfam_in:
            tigrfam_hash.update(tigrfam_in.read().encode())

        with gzip.open(default_map["METACYC"], 'rt') as metacyc_in:
            metacyc_hash.update(metacyc_in.read().encode())

        obs_hash = [ec_hash.hexdigest(), ko_hash.hexdigest(),
                    cog_hash.hexdigest(), pfam_hash.hexdigest(),
                    tigrfam_hash.hexdigest(), metacyc_hash.hexdigest()]

        exp_hash = ["61b2fcd300fd53124c4d7f4b8e97b281",
                    "f1cc419051a23bb60d1015762d221deb",
                    "6012eaf8b2f9e336a725cd97af8cf05d",
                    "b24d1f3cae10efd452b964a8589963e1",
                    "8115f710df156c46908d112bd80def4a",
                    "d441bce3c19effa1e474711f6c6cdbeb"]

        # Check that md5sum values match expected values.
        self.assertEqual(obs_hash, exp_hash)

class convert_table_tests(unittest.TestCase):

    def test_picrust2_to_humann2_merged(self):

        with TemporaryDirectory() as temp_dir:
            outfile = path.join(temp_dir, "test_out")
            convert_picrust2_to_humann2_merged([picrust2_unstrat_in1,
                                                picrust2_unstrat_in2],
                                               outfile)
            obs_out = pd.read_csv(outfile, sep="\t", index_col=0)

        exp_out = pd.read_csv(picrust2_unstrat_exp, sep="\t", index_col=0)

        pd.testing.assert_frame_equal(obs_out, exp_out, check_like=True)

    def test_picrust2_strat_to_humann2_split(self):

        with TemporaryDirectory() as temp_dir:
            outfolder = path.join(temp_dir, "outfiles")
            convert_picrust2_to_humann2([picrust2_strat_in,
                                         picrust2_unstrat_in1], outfolder, True)
            obs_out1 = pd.read_csv(path.join(outfolder, "sample1_humann2-format.tsv"), sep="\t", index_col=0)
            obs_out2 = pd.read_csv(path.join(outfolder, "sample2_humann2-format.tsv"), sep="\t", index_col=0)
            obs_out3 = pd.read_csv(path.join(outfolder, "sample3_humann2-format.tsv"), sep="\t", index_col=0)

        exp_out1 = pd.read_csv(picrust2_strat_exp1, sep="\t", index_col=0)
        exp_out2 = pd.read_csv(picrust2_strat_exp2, sep="\t", index_col=0)
        exp_out3 = pd.read_csv(picrust2_strat_exp3, sep="\t", index_col=0)

        pd.testing.assert_frame_equal(obs_out1, exp_out1, check_like=True)
        pd.testing.assert_frame_equal(obs_out2, exp_out2, check_like=True)
        pd.testing.assert_frame_equal(obs_out3, exp_out3, check_like=True)

    def test_humann2_unstrat_to_picrust2(self):

        with TemporaryDirectory() as temp_dir:

            outfile = path.join(temp_dir, "test_out")
            convert_humann2_to_picrust2([humann2_unstrat_in], outfile, False)

            obs_out = pd.read_csv(outfile, sep="\t", index_col=0)

        exp_out = pd.read_csv(humann2_unstrat_exp, sep="\t", index_col=0)

        pd.testing.assert_frame_equal(obs_out, exp_out, check_like=True)

    def test_humann2_strat_to_picrust2(self):

        with TemporaryDirectory() as temp_dir:

            outfile = path.join(temp_dir, "test_out")
            convert_humann2_to_picrust2([humann2_strat_in], outfile, True)

            obs_out = pd.read_csv(outfile, sep="\t", index_col=[0, 1])

        exp_out = pd.read_csv(humann2_strat_exp, sep="\t", index_col=[0, 1])

        pd.testing.assert_frame_equal(obs_out, exp_out, check_like=True)

    def test_humann2_strat_split_to_picrust2(self):

        with TemporaryDirectory() as temp_dir:

            outfile = path.join(temp_dir, "test_out")
            convert_humann2_to_picrust2([humann2_strat_in_split1,
                                         humann2_strat_in_split2,
                                         humann2_strat_in_split3],
                                        outfile, True)

            obs_out = pd.read_csv(outfile, sep="\t", index_col=[0, 1])

        exp_out = pd.read_csv(humann2_strat_split_exp, sep="\t", index_col=[0, 1])

        pd.testing.assert_frame_equal(obs_out, exp_out, check_like=True)

if __name__ == '__main__':
    unittest.main()
