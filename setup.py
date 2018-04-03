#!/usr/bin/env python

from setuptools import setup
from glob import glob


__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2-alpha.5"
__maintainer__ = "Gavin Douglas"
__email__ = "gavin.douglas@dal.ca"

long_description = ("PICRUSt: Phylogenetic Investigation of Communities by "
                    "Reconstruction of Unobserved States\n\n"
                    "Predictive functional profiling of microbial communities "
                    "using 16S rRNA marker gene sequences. Langille, M. "
                    "G.I.*; Zaneveld, J.*; Caporaso, J. G.; McDonald, D.; "
                    "Knights, D.; a Reyes, J.; Clemente, J. C.; Burkepile, D. "
                    "E.; Vega Thurber, R. L.; Knight, R.; Beiko, R. G.; and "
                    "Huttenhower, C. Nature Biotechnology, 1-10. 8 2013.")


setup(name='PICRUSt',
      version=__version__,
      description=('PICRUSt: Phylogenetic Investigation of Communities by '
                   'Reconstruction of Unobserved States'),
      author=__maintainer__,
      author_email=__email__,
      maintainer=__maintainer__,
      maintainer_email=__email__,
      url='https://github.com/picrust/picrust2/wiki',
      packages=['picrust'],
      scripts=glob('scripts/*py'),
      install_requires=['numpy >= 1.5.1',
                        'cogent == 1.5.3',
			'h5py == 2.7.1',
                        'biom-format >= 2.1.4, < 2.2.0',
                        'future == 0.16'],
      package_data={'precalculated':
		    ['prokaryotic/*',
                     'eukaryotic/*'],
                    'picrust':
                    ['support_files/R/*.R']},
      long_description=long_description)
