#!/usr/bin/env python

__copyright__ = "Copyright 2018-2020, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.2.1-b"

import unittest
import gzip
import hashlib
from os import path
from picrust2.util import read_phylip, read_fasta, TemporaryDirectory
from picrust2.default import (default_ref_dir, default_fasta, default_tree,
                              default_hmm, default_model)
from picrust2.place_seqs import (place_seqs_pipeline, split_ref_study_papara,
                                 run_epa_ng, gappa_jplace_to_newick,
                                 identify_ref_files, parse_jplace)

# Set paths to test files.
test_dir_path = path.join(path.dirname(path.abspath(__file__)), "test_data",
                          "place_seqs")

test_study_seqs = path.join(test_dir_path, "study_seqs_test.fasta")

test_ref_dir = path.join(test_dir_path, "img_centroid_16S_aligned_head30")
test_tree = path.join(test_ref_dir, "img_centroid_16S_aligned_head30.tre")
test_msa = path.join(test_ref_dir, "img_centroid_16S_aligned_head30.fna.gz")
test_hmm = path.join(test_ref_dir, "img_centroid_16S_aligned_head30.hmm")
test_model = path.join(test_ref_dir, "img_centroid_16S_aligned_head30.model")

exp_study_fasta = path.join(test_dir_path, "place_seqs_output",
                            "place_seqs_working",
                            "study_seqs_papara.fasta")

exp_ref_fasta = path.join(test_dir_path, "place_seqs_output",
                          "place_seqs_working",
                          "ref_seqs_papara.fasta")

exp_newick = path.join(test_dir_path, "place_seqs_output",
                       "img_centroid_16S_aligned_head30_placed.tre")

exp_jplace = path.join(test_dir_path, "place_seqs_output",
                       "place_seqs_working",
                       "epa_out", "epa_result.jplace")

test_jplace1 = path.join(test_dir_path, "jplace_files", "epa_result1.jplace")

test_jplace2 = path.join(test_dir_path, "jplace_files", "epa_result2.jplace")


class place_seqs_tests(unittest.TestCase):
    '''Tests for place seqs pipeline (functions in picrust2/place_seqs.py)'''

    def test_default_md5sum(self):
        '''Test that default files match expected md5sum values.'''

        # Calculate md5sum for fasta and treefile respectively.
        fasta_hash = hashlib.md5()
        tree_hash = hashlib.md5()
        hmm_hash = hashlib.md5()
        model_hash = hashlib.md5()

        with gzip.open(default_fasta, 'rt') as fasta_in:
            fasta_hash.update(fasta_in.read().encode())

        with open(default_tree) as tree_in:
            tree_hash.update(tree_in.read().encode())

        with open(default_hmm) as hmm_in:
            hmm_hash.update(hmm_in.read().encode())

        with open(default_model) as model_in:
            model_hash.update(model_in.read().encode())

        # Check that md5sum values match expected values.
        self.assertEqual([fasta_hash.hexdigest(), tree_hash.hexdigest(),
                          hmm_hash.hexdigest(), model_hash.hexdigest()],
                         ['a1675a5c0a2c1941ff2e71c63049a7b4',
                          'f247071837f74c156dc530736cb6d453',
                          'd50b0dac445b5243e86816dbdeadf898',
                          '478b5011ea8aeb0c720e9bb68774fabd'])


    def test_jplace_parse(self):
        '''Test that jplace files are being parsed consistently.'''

        with TemporaryDirectory() as temp_dir:

            test_jplace1_out = temp_dir + path.basename(test_jplace1) +\
                               "_parsed"

            test_jplace2_out = temp_dir + path.basename(test_jplace2) +\
                               "_parsed"

            parse_jplace(test_jplace1, test_jplace1_out)

            parse_jplace(test_jplace2, test_jplace2_out)

            test_jplace1_hash = hashlib.md5()
            test_jplace2_hash = hashlib.md5()

            with open(test_jplace1_out, 'rt') as out1:
                test_jplace1_hash.update(out1.read().encode())

            with open(test_jplace2_out, 'rt') as out2:
                test_jplace2_hash.update(out2.read().encode())

        self.assertEqual(test_jplace1_hash.hexdigest(),
                         test_jplace2_hash.hexdigest())


    def test_gappa_jplace_to_newick(self):
        '''Basic test for gappa_jplace_to_newick function.'''

        # Read in expected newick output.
        exp_newick_in = open(exp_newick).read()

        with TemporaryDirectory() as temp_dir:
            newick_out = path.join(temp_dir, "out.tre")

            gappa_jplace_to_newick(jplace_file=exp_jplace, outfile=newick_out)

            obs_newick_in = open(newick_out).read()

        self.assertEqual(exp_newick_in, obs_newick_in)

    def test_run_epa_ng(self):
        '''Basic test to check whether EPA-NG wrapper can be run. Exact
        matches to a treefile are not checked since slight differences
        are expected depending on different versions.'''

        with TemporaryDirectory() as temp_dir:
            run_epa_ng(tree=test_tree,
                       model=test_model,
                       ref_msa_fastafile=exp_ref_fasta,
                       study_msa_fastafile=exp_study_fasta,
                       out_dir=temp_dir)

    def test_run_place_seqs_pipeline(self):
        '''Basic test of full place seqs pipeline. As for EPA-NG, exact
        matches to a treefile are not checked since slight differences
        are expected depending on different versions.'''

        with TemporaryDirectory() as temp_dir:
            tmp_tree = path.join(temp_dir, "out.tre")



            place_seqs_pipeline(study_fasta=test_study_seqs,
                                ref_dir=test_ref_dir,
                                out_tree=tmp_tree,
                                threads=1,
                                out_dir=temp_dir,
                                chunk_size=5000,
                                print_cmds=False)

    def test_identify_ref_files(self):
        '''Test for reference files being identified correctly.'''

        expected_files = [default_fasta, default_tree, default_hmm,
                          default_model]

        identified_files = identify_ref_files(default_ref_dir)

        self.assertEqual(expected_files, identified_files)


if __name__ == '__main__':
    unittest.main()
