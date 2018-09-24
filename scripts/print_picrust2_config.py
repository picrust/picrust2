#!/usr/bin/env python

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.0.0-b.9"

from sys import platform, version as python_version, executable

try:
    from numpy import __version__ as numpy_lib_version
except ImportError:
    numpy_lib_version = "ERROR: Not installed - this is required! " +\
                        "(This will also cause the BIOM library to not be importable.)"

try:
    from biom import __version__ as biom_lib_version
except ImportError:
    biom_lib_version = "ERROR: Can't find the BIOM library code (or numpy) " +\
                       "- is it installed and in your $PYTHONPATH?"

try:
    from picrust2 import __version__ as picrust_lib_version
except ImportError:
    picrust_lib_version = "ERROR: Can't find the PICRUSt library code " +\
                          "- is it installed and in your $PYTHONPATH?"

import argparse

parser = argparse.ArgumentParser(

    description="A simple script that prints out the PICRUSt config " +
                "settings and does some sanity checks.",

    formatter_class=argparse.RawDescriptionHelpFormatter)


def get_script_version():
    return __version__


def print_picrust_config():
    system_info = [
     ("Platform", platform),
     ("Python/GCC version", python_version.replace('\n', ' ')),
     ("Python executable", executable)]

    max_len = max([len(e[0]) for e in system_info])
    print("\nSystem information")
    print("==================")
    for v in system_info:
        print("%*s:\t%s" % (max_len, v[0], v[1]))

    version_info = [
     ("NumPy version", numpy_lib_version),
     ("biom-format version", biom_lib_version),
     ("PICRUSt version", picrust_lib_version),
     ("PICRUSt script version", get_script_version()), ]

    max_len = max([len(e[0]) for e in version_info])

    print("\nDependency versions")
    print("===================")
    for v in version_info:
        print("%*s:\t%s" % (max_len, v[0], v[1]))
    print("")


def main():

    print_picrust_config()


if __name__ == '__main__':
    main()
