#!/usr/bin/env python

from setuptools import setup
from glob import glob

__copyright__ = "Copyright 2018-2019, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.1.4-b"
__maintainer__ = "Gavin Douglas"

long_description = ("Please visit the google group here if you have questions: "
                    "https://groups.google.com/forum/#!forum/picrust-users. "
                    "Pre-print: Douglas et al. 2019. PICRUSt2: An improved and "
                    "extensible approach for metagenome inference. bioRxiv. "
                    "doi: https://doi.org/10.1101/672295"
                    )

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
