#!/usr/bin/env python

import unittest
from os import path
from picrust2.util import system_call_check, TemporaryDirectory

# Paths to input files.
test_dir_path = path.join(path.dirname(path.abspath(__file__)))

test_study_seqs = path.join(test_dir_path, "test_data", "place_seqs",
                            "study_seqs_test.fasta")

test_ref_dir = path.join(test_dir_path, "test_data", "place_seqs",
                         "img_centroid_16S_aligned_head30")

test_known_marker = path.join(test_dir_path, "test_data", "workflow",
                              "workflow_known_marker.tsv.gz")

test_known_traits = path.join(test_dir_path, "test_data", "workflow",
                              "workflow_known_traits.tsv.gz")

test_seq_abun_tsv = path.join(test_dir_path, "test_data", "workflow",
                              "workflow_seq_abun.tsv.gz")

test_seq_abun_biom = path.join(test_dir_path, "test_data", "workflow",
                               "workflow_seq_abun.biom")


class workflow_test(unittest.TestCase):
    '''Test for whether full pipeline runs without error. This is intended to
    catch incompatabilities between the different scripts that the other tests
    might not catch.'''

    def test_help_calls(self):
        '''Quick test that help commands can be run for all scripts, which
        is a check that there are no obvious argparse syntax errors.'''

        system_call_check("add_descriptions_oldIMG.py -h")
        system_call_check("convert_table.py -h")
        system_call_check("hsp.py -h")
        system_call_check("metagenome_pipeline.py -h")
        system_call_check("picrust2_pipeline_oldIMG.py -h")
        system_call_check("place_seqs.py -h")
        system_call_check("pathway_pipeline.py -h")
        system_call_check("shuffle_predictions.py -h")

    def test_full_pipeline_tsv(self):
        '''Test that full pipeline can be run without error with
        TSV sequence abundance table.'''

        with TemporaryDirectory() as temp_dir:

            out_tree = path.join(temp_dir, "out.tre")

            system_call_check("place_seqs.py -s " + test_study_seqs + " -r " +
                              test_ref_dir + " -o " + out_tree)

            traits_predict = path.join(temp_dir, "hsp_out.tsv.gz")
            marker_predict = path.join(temp_dir, "hsp_out_marker.tsv.gz")

            system_call_check("hsp.py -t " + out_tree +
                              " --observed_trait_table " + test_known_traits +
                              " -n -o " + traits_predict)

            system_call_check("hsp.py -t " + out_tree +
                              " --observed_trait_table " + test_known_marker +
                              " -n -o " + marker_predict)

            shuffled_out = path.join(temp_dir, "shuffled_out")

            system_call_check("shuffle_predictions.py -i " + traits_predict +
                              " -o " + shuffled_out +
                              " -s " + str(14141) +
                              " -r " + str(3))

            metagenome_out = path.join(temp_dir, "meta_out")

            system_call_check("metagenome_pipeline.py -i " + test_seq_abun_tsv +
                              " -f " + traits_predict +
                              " --strat_out " +
                              " -m " + marker_predict +
                              " -o " + metagenome_out)

            metagenome_outfile = path.join(metagenome_out,
                                           "pred_metagenome_unstrat.tsv.gz")

            system_call_check("pathway_pipeline.py -i " + metagenome_outfile +
                              " -o " + temp_dir)

    def test_picrust2_pipeline_script(self):
        '''Test that full pipeline can be run successfully with
        picrust2_pipeline_oldIMG.py'''

        with TemporaryDirectory() as temp_dir:

            out_dir = path.join(temp_dir, "pipeline_out")

            system_call_check("picrust2_pipeline_oldIMG.py -s " + test_study_seqs +
                              " -i " + test_seq_abun_tsv +
                              " -o " + out_dir +
                              " -r " + test_ref_dir +
                              " -p 1" +
                              " --custom_trait_tables " + test_known_traits +
                              " --marker_gene_table " + test_known_marker +
                              " --reaction_func " + test_known_traits +
                              " --max_nsti 1.9" +
                              " --min_reads 2" +
                              " --min_samples 2" +
                              " --skip_minpath" +
                              " --no_gap_fill" +
                              " --coverage" +
                              " --remove_intermediate" +
                              " --verbose")

    def test_picrust2_pipeline_script_per_seq_contrib_strat(self):
        '''Test that full pipeline can be run successfully with
        picrust2_pipeline_oldIMG.py with the --per_sequence_contrib and --stratified
        options.'''

        with TemporaryDirectory() as temp_dir:

            out_dir = path.join(temp_dir, "pipeline_out")

            system_call_check("picrust2_pipeline_oldIMG.py -s " + test_study_seqs +
                              " -i " + test_seq_abun_tsv +
                              " -o " + out_dir +
                              " -r " + test_ref_dir +
                              " -p 1" +
                              " --custom_trait_tables " + test_known_traits +
                              " --marker_gene_table " + test_known_marker +
                              " --reaction_func " + test_known_traits +
                              " --max_nsti 1.9" +
                              " --min_reads 2" +
                              " --min_samples 2" +
                              " --skip_minpath" +
                              " --no_gap_fill" +
                              " --coverage" +
                              " --remove_intermediate" +
                              " --stratified" +
                              " --per_sequence_contrib"
                              " --verbose")

    def test_picrust2_pipeline_script_per_seq_contrib_strat_skip_norm(self):
        '''Test that full pipeline can be run successfully with
        picrust2_pipeline_oldIMG.py with the --per_sequence_contrib and --stratified
        options as well as --skip_norm.'''

        with TemporaryDirectory() as temp_dir:

            out_dir = path.join(temp_dir, "pipeline_out")

            system_call_check("picrust2_pipeline_oldIMG.py -s " + test_study_seqs +
                              " -i " + test_seq_abun_tsv +
                              " -o " + out_dir +
                              " -r " + test_ref_dir +
                              " -p 1" +
                              " --custom_trait_tables " + test_known_traits +
                              " --marker_gene_table " + test_known_marker +
                              " --reaction_func " + test_known_traits +
                              " --max_nsti 1.9" +
                              " --min_reads 2" +
                              " --min_samples 2" +
                              " --skip_minpath" +
                              " --no_gap_fill" +
                              " --coverage" +
                              " --remove_intermediate" +
                              " --stratified" +
                              " --per_sequence_contrib" +
                              " --skip_norm" +
                              " --verbose")


if __name__ == '__main__':
    unittest.main()
