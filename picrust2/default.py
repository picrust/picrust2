#!/usr/bin/env python

from os import path

project_dir = path.dirname(path.abspath(__file__))

default_ref_dir_bac = path.join(project_dir, "default_files", "bacteria",
                            "bac_ref")

default_ref_dir_arc = path.join(project_dir, "default_files", "archaea",
                            "arc_ref")

default_fasta_bac = path.join(default_ref_dir_bac, "bac_ref.fna")

default_fasta_arc = path.join(default_ref_dir_arc, "arc_ref.fna")

default_tree_bac = path.join(default_ref_dir_bac, "bac_ref.tre")

default_tree_arc = path.join(default_ref_dir_arc, "arc_ref.tre")

default_hmm_bac = path.join(default_ref_dir_bac, "bac_ref.hmm")

default_hmm_arc = path.join(default_ref_dir_arc, "arc_ref.hmm")

default_model_bac = path.join(default_ref_dir_bac, "bac_ref.model")

default_model_arc = path.join(default_ref_dir_arc, "arc_ref.model")

default_raxml_info_bac = path.join(default_ref_dir_bac, "bac_ref.raxml_info")

default_raxml_info_arc = path.join(default_ref_dir_arc, "arc_ref.raxml_info")

default_regroup_map = path.join(project_dir, "default_files",
                                "pathway_mapfiles",
                                "ec_level4_to_metacyc_rxn_new.tsv")

default_pathway_map = path.join(project_dir, "default_files",
                                "pathway_mapfiles",
                                "metacyc_pathways_structured_filtered_v24_subreactions.txt")

#fungi_pathway_map = path.join(project_dir, "default_files", "pathway_mapfiles",
#                              "metacyc_path2rxn_struc_filt_fungi.txt")

# Inititalize default trait table files for hsp.py.
bacteria_dir = path.join(project_dir, "default_files", "bacteria")

default_tables_bac = {"16S": path.join(bacteria_dir, "16S.txt.gz"),

                  "EC": path.join(bacteria_dir, "ec.txt.gz"),

                  "KO": path.join(bacteria_dir, "ko.txt.gz"),

                  "GO": path.join(bacteria_dir, "go.txt.gz"),

                  "PFAM": path.join(bacteria_dir, "pfam.txt.gz"),

                  "BIGG": path.join(bacteria_dir, "bigg_reaction.txt.gz"),

                  "CAZY": path.join(bacteria_dir, "cazy.txt.gz"),

                  "GENE_NAMES": path.join(bacteria_dir, "preferred_name.txt.gz")}

archaea_dir = path.join(project_dir, "default_files", "archaea")

default_tables_arc = {"16S": path.join(archaea_dir, "16S.txt.gz"),

                  "EC": path.join(archaea_dir, "ec.txt.gz"),

                  "KO": path.join(archaea_dir, "ko.txt.gz"),

                  "GO": path.join(archaea_dir, "go.txt.gz"),

                  "PFAM": path.join(archaea_dir, "pfam.txt.gz"),

                  "BIGG": path.join(archaea_dir, "bigg_reaction.txt.gz"),

                  "CAZY": path.join(archaea_dir, "cazy.txt.gz"),

                  "GENE_NAMES": path.join(archaea_dir, "preferred_name.txt.gz")}


# Initialize default mapfiles to be used with add_descriptions.py
map_dir = path.join(project_dir, "default_files", "description_mapfiles")

default_map = {"METACYC": path.join(map_dir,
                                    "metacyc-pwy_name.txt.gz"),

               "EC": path.join(map_dir, "ec_name.txt.gz"),

               "KO": path.join(map_dir, "ko_name.txt.gz"),

               "GO": path.join(map_dir, "map_go_name.txt.gz")}
