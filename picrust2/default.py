#!/usr/bin/env python

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.0.0-b.1"

from picrust2.util import get_picrust_project_dir
from os import path

# Default support files packaged with PICRUSt2.
project_dir = get_picrust_project_dir()

default_fasta = path.join(project_dir, "default_files", "prokaryotic",
                          "img_centroid_16S_aligned.fna")

default_tree = path.join(project_dir, "default_files", "prokaryotic",
                         "img_centroid_16S_aligned.tree")


# Inititalize default trait table files for hsp.py.
prokaryotic_dir = path.join(project_dir, "default_files", "prokaryotic")

default_tables = {"16S": path.join(prokaryotic_dir,
                                   "16S_counts_mean_round_var.txt.gz"),

                  "COG": path.join(prokaryotic_dir,
                                   "cog_counts_mean_round_var.txt.gz"),

                  "EC": path.join(prokaryotic_dir,
                                  "ec_counts_mean_round_var.txt.gz"),

                  "KO": path.join(prokaryotic_dir,
                                  "ko_counts_mean_round_var.txt.gz"),

                  "PFAM": path.join(prokaryotic_dir,
                                    "pfam_counts_mean_round_var.txt.gz"),

                  "TIGRFAM": path.join(prokaryotic_dir,
                                       "tfam_counts_mean_round_var.txt.gz")}
