#!/usr/bin/env python

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.0.2-b"

import unittest
from os import path
from tempfile import TemporaryDirectory
from picrust2.util import get_picrust_project_dir, system_call_check

# Paths to input files.
test_dir_path = path.join(get_picrust_project_dir(), "tests")

test_study_seqs = path.join(test_dir_path, "test_data", "place_seqs",
                            "study_seqs_test.fasta")

test_tree = path.join(test_dir_path, "test_data", "place_seqs",
                      "img_centroid_16S_aligned_head30.tre")

test_msa = path.join(test_dir_path, "test_data", "place_seqs",
                     "img_centroid_16S_aligned_head30.fna")

test_hmm = path.join(test_dir_path, "test_data", "place_seqs",
                     "img_centroid_16S_aligned_head30.hmm")

test_known_marker = path.join(test_dir_path, "test_data", "workflow",
                              "workflow_known_marker.tsv")

test_known_traits = path.join(test_dir_path, "test_data", "workflow",
                              "workflow_known_traits.tsv")

test_seq_abun_tsv = path.join(test_dir_path, "test_data", "workflow",
                             "workflow_seq_abun.tsv")

test_seq_abun_biom = path.join(test_dir_path, "test_data", "workflow",
                               "workflow_seq_abun.biom")

class workflow_test(unittest.TestCase):
    '''Test for whether full pipeline runs without error. This is intended to
    catch incompatabilities between the different scripts that the other tests
    might not catch.'''

    def test_full_pipeline_tsv(self):
        '''Test that full pipeline can be run without error with
        TSV sequence abundance table.'''

        with TemporaryDirectory() as temp_dir:

            out_tree = path.join(temp_dir, "out.tre")

            system_call_check("place_seqs.py -s " + test_study_seqs + " -r " +
                              test_msa + " -t " + test_tree + " --hmm " +
                              test_hmm + " -o " + out_tree)

            hsp_out_prefix = path.join(temp_dir, "hsp_out")
            hsp_out_prefix_marker = path.join(temp_dir, "hsp_out_marker")

            system_call_check("hsp.py -t " + out_tree +
                " --observed_trait_table " + test_known_traits + " -n -c " +
                "-o " + hsp_out_prefix)

            system_call_check("hsp.py -t " + out_tree +
                " --observed_trait_table " + test_known_marker + " -n -c " +
                "-o " + hsp_out_prefix_marker)

            traits_predict = path.join(temp_dir, hsp_out_prefix +
                                       ".tsv")

            marker_predict = path.join(temp_dir, hsp_out_prefix_marker +
                                       ".tsv")

            metagenome_out = path.join(temp_dir, "meta_out")

            system_call_check("metagenome_pipeline.py -i " + test_seq_abun_tsv +
                              " -f " + traits_predict + " --strat_out -m " +
                              marker_predict + " -o " + metagenome_out)

            metagenome_outfile = path.join(metagenome_out,
                                           "pred_metagenome_strat.tsv")

            system_call_check("run_minpath.py -i " + metagenome_outfile +
                              " -o " + temp_dir)

    def test_picrust2_pipeline_script(self):
        '''Test that full pipeline can be run successfully with
        picrust2_pipeline.py'''

        with TemporaryDirectory() as temp_dir:

            out_dir = path.join(temp_dir, "pipeline_out")

            system_call_check("picrust2_pipeline.py -s " + test_study_seqs +
                              " -i " + test_seq_abun_tsv + " -o " + out_dir +
                              " -r " + test_msa + " -t " + test_tree +
                              " --hmm " + test_hmm +
                              " --custom_trait_tables " + test_known_traits +
                              " --marker_gene_table " + test_known_marker)

if __name__ == '__main__':
    unittest.main()
