#!/usr/bin/env python

from setuptools import setup
from glob import glob

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.0.4-b"
__maintainer__ = "Gavin Douglas"

long_description = ("Please visit the google group here if you have questions: "
                    "https://groups.google.com/forum/#!forum/picrust-users. "
                    "Citation: Phylogenetic Investigation of Communities by "
                    "Reconstruction of Unobserved States\n\n"
                    "Predictive functional profiling of microbial communities "
                    "using 16S rRNA marker gene sequences. Langille, M. "
                    "G.I.*; Zaneveld, J.*; Caporaso, J. G.; McDonald, D.; "
                    "Knights, D.; a Reyes, J.; Clemente, J. C.; Burkepile, D. "
                    "E.; Vega Thurber, R. L.; Knight, R.; Beiko, R. G.; and "
                    "Huttenhower, C. Nature Biotechnology, 1-10. 8 2013.")

setup(name='PICRUSt2',
      version=__version__,
      description=('PICRUSt: Phylogenetic Investigation of Communities by '
                   'Reconstruction of Unobserved States'),
      maintainer=__maintainer__,
      url='https://github.com/picrust/picrust2/wiki',
      packages=['picrust2'],
      scripts=glob('scripts/*py'),
      install_requires=['numpy',
			'h5py',
                        'joblib',
                        'biom-format'],
      package_data={'picrust2':
                    ['MinPath/MinPath12hmp.py',
                     'Rscripts/*R']},
      long_description=long_description)
