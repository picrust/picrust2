#!/usr/bin/env python

from __future__ import division
from collections import defaultdict
from picrust.util import system_call_check

__author__ = "Gavin Douglas"
__copyright__ = "Copyright 2018, The PICRUSt Project"
__credits__ = ["Gavin Douglas", "Morgan Langille"]
__license__ = "GPL"
__version__ = "2-alpha.4"

# Class that contains pathway abundances for each sample.
class pathway_counts():

	def __init__(self, filename):
		self.__counts = defaultdict(float)

	def set_pathway_abun(self, path_id, abun):
		self.__counts[path_id] = abun

	def return_pathway_abun(self, path_id):
		return self.__counts[path_id]


def minpath_wrapper(sample_id, biom_in, minpath_map, tmp_dir, functions,
	                print_opt):
	'''Read in sample_id, tmp_dir, list of functions to loop through.'''

	# Define MinPath input and outout filenames.
	minpath_in = str(tmp_dir + "/" + sample_id + "_minpath_in.txt")
	minpath_report = str(tmp_dir + "/" + sample_id + "_minpath_report.txt")
	minpath_details = str(tmp_dir + "/" + sample_id + "_minpath_details.txt")
	minpath_mps = str(tmp_dir + "/" + sample_id + "_minpath.mps")
	minpath_output = open(str(tmp_dir + "/" + sample_id + 
                          "_minpath_out.txt"), "w")

	id_minpath_fh = open(minpath_in, "w")

	# Counter to give each "read" in MinPath input a different id.
	func_num = 0

	for func_id in functions:
    	# Get count of each sequence in sample and write that sequence out
    	# along with count if non-zero abundance.
		func_count = int(biom_in.get_value_by_ids(obs_id=func_id,
                                                      samp_id=sample_id))
		# If 0 then skip.
		if func_count == 0:
			continue

		id_minpath_fh.write(func_id + "\t" + str(func_count) + "\n")

	id_minpath_fh.close()

	# Run MinPath on this sample.
	minpath_cmd = "MinPath12hmp.py -any " + minpath_in + " -map " +\
                  minpath_map + " -report " + minpath_report +\
                  " -details " + minpath_details + " -mps " + minpath_mps

	system_call_check(minpath_cmd, print_out=print_opt,
                      stdout=minpath_output)

	# Read through MinPath report and keep track of pathways identified
	# to be present.
	path_present = set()

	with open(minpath_report, "r") as minpath_report_in:
		for line in minpath_report_in:
			line_split = line.split()

			if int(line_split[7]) == 1:
				path_present.add(line_split[-1])

	# Now read in details file and take abundance of pathway to be
	# harmonic mean of gene families in pathway to be abundance of pathway.
	# Abundances of 0 will be added in for gene families not found.

	# Initialize dictionary that will contain pathway abundances.
	path_abun = defaultdict(float)

	# Initialize dictionary that will contain gene family abundance per
	# pathway.
	gf_abundances = {}

	# Boolean specifying that pathway in details file was called as
	# present by MinPath.
	present = False

	with open(minpath_details, "r") as minpath_details_in:
		for line in minpath_details_in:
			line_split = line.split()

			# If line starts with "path" then keep track of pathway name if
			# it was called as present in report file.
			if line_split[0] == "path":
				if line_split[-1] not in path_present:
					present = False
					continue
        
				present = True
				current_pathway = line_split[-1]

				# Initialize list containing gene family abundances.
				gf_abundances[current_pathway] = []

				# Add in abundances of 0 for missing genes.
				for i in range(int(line_split[3]) - int(line_split[5])):
					gf_abundances[current_pathway] += [0]

			# If line does not start with "path" then only proceed if current
			# pathway is present.
			elif present:
				gf_abundances[current_pathway] += [int(float(line_split[2]))]

	# Loop through all pathways present and get mean of 1/2 most abundant.
	for pathway in gf_abundances.keys():

		# Like HUMAnN2, sort enzyme reactions, take second half, and get 
		# their mean abundance.
		sorted_gf_abundances = sorted(gf_abundances[pathway])
		gf_abundances_subset = sorted_gf_abundances[int(len(sorted_gf_abundances)/ 2):]
		pathway_abun = sum(gf_abundances_subset)/len(gf_abundances_subset)

		path_abun[pathway] = pathway_abun

	return(path_abun)
