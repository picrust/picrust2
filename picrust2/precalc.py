#!/usr/bin/env python

__license__ = "GPL"
__version__ = "2-alpha.8"

from picrust2.util import get_picrust_project_dir
from os import path

# Default precalculated files and support files packaged with PICRUSt2.
project_dir = get_picrust_project_dir()

default_fasta = path.join(project_dir, "precalculated", "prokaryotic",
                          "img_centroid_16S_aligned.fna")

default_tree = path.join(project_dir, "precalculated", "prokaryotic",
                         "img_centroid_16S_aligned.tree")
