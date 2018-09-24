#!/usr/bin/env python

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.0.0-b.9"

from picrust2.util import get_picrust_project_dir
from os import path

# Default support files packaged with PICRUSt2.
project_dir = get_picrust_project_dir()

default_fasta = path.join(project_dir, "default_files", "prokaryotic",
                          "reference.fna")

default_tree = path.join(project_dir, "default_files", "prokaryotic",
                         "reference.tre")

default_regroup_map = path.join(project_dir, "default_files",
                                "pathway_mapfiles",
                                "ec_level4_to_metacyc_rxn.tsv")

default_pathway_map = path.join(project_dir, "default_files",
                                "pathway_mapfiles",
                                "metacyc_path2rxn_struc_filt_pro.txt")

# Inititalize default trait table files for hsp.py.
prokaryotic_dir = path.join(project_dir, "default_files", "prokaryotic")

default_tables = {"16S": path.join(prokaryotic_dir, "16S.txt.gz"),

                  "COG": path.join(prokaryotic_dir, "cog.txt.gz"),

                  "EC": path.join(prokaryotic_dir, "ec.txt.gz"),

                  "KO": path.join(prokaryotic_dir, "ko.txt.gz"),

                  "PFAM": path.join(prokaryotic_dir, "pfam.txt.gz"),

                  "TIGRFAM": path.join(prokaryotic_dir, "tigrfam.txt.gz")}

# Initialize default mapfiles to be used with add_descriptions.py
map_dir = path.join(project_dir, "default_files", "description_mapfiles")

default_map =    {"METACYC": path.join(map_dir, "metacyc_pathways_info_prokaryotes.txt.gz"),

                  "COG": path.join(map_dir, "cog_info.tsv.gz"),

                  "EC": path.join(map_dir, "ec_level4_info.tsv.gz"),

                  "KO": path.join(map_dir, "ko_info.tsv.gz"),

                  "PFAM": path.join(map_dir, "pfam_info.tsv.gz"),

                  "TIGRFAM": path.join(map_dir, "tigrfam_info.tsv.gz")}

