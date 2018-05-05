#!/usr/bin/env python

__license__ = "GPL"
__version__ = "2-alpha.9"

from picrust2.util import get_picrust_project_dir
from os import path

# Default precalculated files and support files packaged with PICRUSt2.
project_dir = get_picrust_project_dir()

default_fasta = path.join(project_dir, "precalculated", "prokaryotic",
                          "img_centroid_16S_aligned.fna")

default_tree = path.join(project_dir, "precalculated", "prokaryotic",
                         "img_centroid_16S_aligned.tree")


# Inititalize default trait table files for hsp.py.
prokaryotic_dir = path.join(project_dir, "precalculated", "prokaryotic")

default_tables = {"16S" : path.join(prokaryotic_dir,
                                    "16S_counts_mean_round_var.txt"),

                  "COG" : path.join(prokaryotic_dir,
                                    "cog_counts_mean_round_var.txt"),

                  "EC" : path.join(prokaryotic_dir,
                                   "ec_counts_mean_round_var.txt"),

                  "KO" : path.join(prokaryotic_dir,
                                   "ko_counts_mean_round_var.txt"),

                  "PFAM" : path.join(prokaryotic_dir,
                                     "pfam_counts_mean_round_var.txt"),

                  "TIGRFAM": path.join(prokaryotic_dir,
                                       "tfam_counts_mean_round_var.txt")}
