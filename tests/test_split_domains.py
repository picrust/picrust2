#!/usr/bin/env python

import unittest
import pandas as pd
from os import path
import os
from picrust2.util import make_output_dir, TemporaryDirectory
from picrust2.default import (default_ref_dir_bac, default_ref_dir_arc,
                              default_tables_bac, default_tables_arc)
from picrust2.place_seqs import place_seqs_pipeline
from picrust2.wrap_hsp import castor_hsp_workflow
from picrust2.split_domains import get_lowest_nsti

test_dir_path = path.join(path.dirname(path.abspath(__file__)), "test_data",
                          "split_domains")

study_seqs = path.join(test_dir_path, "seqs.fna")

exp_combined_marker_nsti = path.join(test_dir_path, "output", "combined_marker_predicted_and_nsti.tsv.gz")

@unittest.skipIf(os.environ.get("CI"),
                "Skip this one in Circle CI because the reference needs ~27 GB RAM. 
                Still running it on other instances of pytest.")
class split_domains_placement_test(unittest.TestCase):
    '''Check that the sequences in the test input file get assigned to the expected
    domain when they're placed in the bacteria and archaea reference trees.'''

    def test_place_seqs(self):
        '''Place the sequences into the bacterial and archaeal trees and check which
        has the lowest NSTI.'''

        with TemporaryDirectory() as temp_dir:

            # place sequences in the archaeal tree - doing this one first as the files are smaller so if it fails then it should be quicker
            out_tree_arc = path.join(temp_dir, "arc.tre")
            place_seqs_dir_arc = path.join(temp_dir, "place_seqs_arc")
            make_output_dir(place_seqs_dir_arc)

            place_seqs_pipeline(study_fasta=study_seqs,
                                ref_dir=default_ref_dir_arc,
                                placement_tool="epa-ng",
                                out_tree=out_tree_arc,
                                threads=1,
                                out_dir=place_seqs_dir_arc,
                                min_align=0.8,
                                chunk_size=500,
                                verbose=False)

            # now in the bacterial tree
            out_tree_bac = path.join(temp_dir, "bac.tre")
            place_seqs_dir_bac = path.join(temp_dir, "place_seqs_bac")
            make_output_dir(place_seqs_dir_bac)

            place_seqs_pipeline(study_fasta=study_seqs,
                                ref_dir=default_ref_dir_bac,
                                placement_tool="epa-ng",
                                out_tree=out_tree_bac,
                                threads=1,
                                out_dir=place_seqs_dir_bac,
                                min_align=0.8,
                                chunk_size=500,
                                verbose=False)

            # now get the NSTI predictions
            hsp_table_arc, ci_table_arc = castor_hsp_workflow(tree_path=out_tree_arc,
                                                              trait_table_path=default_tables_arc["16S"],
                                                              hsp_method="mp",
                                                              calc_nsti=True,
                                                              num_proc=1,
                                                              ran_seed=100)
            
            marker_nsti_arc = path.join(temp_dir, "arc_marker_predicted_and_nsti.tsv")
            hsp_table_arc.to_csv(marker_nsti_arc, sep="\t", index_label="sequence")

            hsp_table_bac, ci_table_bac = castor_hsp_workflow(tree_path=out_tree_bac,
                                                              trait_table_path=default_tables_bac["16S"],
                                                              hsp_method="mp",
                                                              calc_nsti=True,
                                                              num_proc=1,
                                                              ran_seed=100)
            
            marker_nsti_bac = path.join(temp_dir, "bac_marker_predicted_and_nsti.tsv")
            hsp_table_bac.to_csv(marker_nsti_bac, sep="\t", index_label="sequence")

            # and compare them
            obs, nsti_bac, nsti_arc = get_lowest_nsti(marker_nsti_bac, marker_nsti_arc, "bac", "arc")

            exp = pd.read_csv(exp_combined_marker_nsti, sep="\t", index_col="sequence")

            pd.testing.assert_frame_equal(obs.sort_index(), exp.sort_index(),
                                            check_dtype=False, atol=1e-4)


if __name__ == '__main__':
    unittest.main()
