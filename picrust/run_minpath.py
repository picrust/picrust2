#!/usr/bin/env python
from __future__ import division
from collections import defaultdict

__author__ = "Gavin Douglas"
__copyright__ = "Copyright 2018, The PICRUSt Project"
__credits__ = ["Gavin Douglas", "Morgan Langille"]
__license__ = "GPL"
__version__ = "1.1.3"
__maintainer__ = "Gavin Douglas"
__email__ = "gavinmdouglas@gmail.com"
__status__ = "Development"


# Class that contains pathway abundances for each sample.
class pathway_counts():

    def __init__(self, filename):
        self.__counts = defaultdict(float)

    def set_pathway_abun(self, path_id, abun):
        self.__counts[path_id] = abun

    def return_pathway_abun(self, path_id):
        return self.__counts[path_id]


def harmonic_mean(in_num):
    """
    Return the harmonic mean for a list of numbers.
    """
    recip_sum = sum((1.0/num) for num in in_num)
    h_mean = len(in_num)/recip_sum

    return h_mean
