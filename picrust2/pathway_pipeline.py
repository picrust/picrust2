#!/usr/bin/env python

__copyright__ = "Copyright 2018-2019, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.2.1-b"

import sys
from collections import defaultdict
from joblib import Parallel, delayed
from os import path
import scipy.stats
import pandas as pd
import numpy as np
import copy
from picrust2.util import (system_call_check, check_files_exist,
                           make_output_dir, read_seqabun)
from picrust2.metagenome_pipeline import (strat_funcs_by_samples,
                                          metagenome_contributions,
                                          contrib_to_unstrat)


class PathwaysDatabase:
    '''Holds all of the reactions/pathways data from the file provided.
    This class was taken from HUMAnN2 v0.11.1 and was modified only slightly so
    that it could be used without additional classes and functions defined in
    HUMAnN2.'''

    def _is_optional_reaction(self, item, reaction_names=[]):
        '''Check if this reaction is optional.'''

        # If item matches reaction name perfectly then it can't be optional.
        # Similarly, some reactions can begin with "--" if a reaction begins
        # with that it likely would not be optional.
        if item in reaction_names or item[0:2] == "--":
            return False
        elif item[0] == "-":
            # Otherwise call it as optional if the item starts with "-".
            return True
        else:
            return False

    def _find_reaction_list_and_key_reactions(self, items,
                                              reaction_names=None):
        '''Find the reactions in the pathways items and also the key
        reactions.'''

        reaction_list = []
        key_reactions = []

        for item in items:

            # Ignore items that are not reactions as they are part of pathway
            # structure.
            if item not in ["(", ")", "", "+", ","]:

                # Check if this reaction is optional and remove the first
                # character if so.
                if self._is_optional_reaction(item, reaction_names):
                    item = item[1:]
                else:
                    # Record that this is a key reaction (i.e. not optional).
                    key_reactions.append(item)

                reaction_list.append(item)

        return reaction_list, key_reactions

    def _find_structure(self, items, reaction_names=None):
        '''Find the structure of the pathway from the string.'''

        # Initialize dictionary that will contain pathway structure (which can
        # include multiple levels).
        structure = [" "]
        levels = {0: structure}
        current_level = 0

        # Loop through items (which can include reactions and syntax of the
        # pathway structure)
        for item in items:
            if item:
                # Check if the item name indicates an optional reaction.
                # If so, remove the "-" at the beginning of the name.
                if self._is_optional_reaction(item, reaction_names):
                    item = item[1:]

                # Check if this is the start of a list.
                if item == "(":

                    # Initialize empty list at next level.
                    new_list = [" "]
                    levels[current_level + 1] = new_list

                    # Add the new list to the structure.
                    levels[current_level].append(new_list)

                    # Update the current level.
                    current_level += 1

                # Check if this is the end of a list.
                elif item == ")":
                    # If so, update the current level to close the list.
                    current_level -= 1
                # Check if this is a delimiter
                elif item in ["+", ","]:
                    # Update the delimiter at the beginning of the list.
                    levels[current_level][0] = item
                else:
                    levels[current_level].append(item)

        return structure

    def _set_pathways_structure(self, reactions, reaction_names=None):
        '''Determine the pathways structure from the input string.'''

        for pathway in reactions:
            # Check if the item is a list of items.
            if isinstance(reactions[pathway], list):
                reactions[pathway] = " ".join(reactions[pathway])

            # Split the reactions information by the structured pathways
            # delimiter.
            reactions[pathway] = reactions[pathway].split(" ")

            # Find and store the structure for the pathway.
            structure = self._find_structure(reactions[pathway],
                                             reaction_names)
            self.__pathways_structure[pathway] = structure

            # Find the list of reactions and the key reactions.
            reaction_list, key_reactions = self._find_reaction_list_and_key_reactions(reactions[pathway],
                                                                                      reaction_names)

            # Store the list of key reactions for the pathway.
            self.__key_reactions[pathway] = key_reactions

            # Update the reactions dictionary to contain the list of reactions
            # instead of the structure string.
            reactions[pathway] = reaction_list

        return reactions

    def _store_pathways(self, reactions):
        '''Create the dictionaries of reactions to pathways and pathways to
        reactions.'''

        for pathway in reactions:
            for reaction in reactions[pathway]:
                self.__pathways_to_reactions[pathway] = self.__pathways_to_reactions.get(pathway, []) + [reaction]
                self.__reactions_to_pathways[reaction] = self.__reactions_to_pathways.get(reaction, []) + [pathway]

    def __init__(self, database=None, reaction_names=[]):
        '''Load in the pathways data from the database file.'''

        self.__pathways_to_reactions = {}
        self.__reactions_to_pathways = {}
        self.__pathways_structure = {}
        self.__key_reactions = {}

        if database is not None:

            # Check that database file exists.
            check_files_exist([database])

            file_handle = open(database, "rt")

            line = file_handle.readline()

            # The database is expected to contain a single line per pathway.
            # This line begins with the pathway name and is followed by all
            # reactions. Alternatively it could also contain a since pathway
            # and reaction link per line if the pathways aren't structured.
            reactions = defaultdict(list)
            structured_pathway = False
            while line:
                data = line.strip().split("\t")

                if len(data) > 1:

                    # Remove pathway from this list.
                    pathway = data.pop(0)

                    # Add new key-value pair in reactions dict of pathway to
                    # all reactions.
                    reactions[pathway] += data

                    # Check to see if this pathway has structure.
                    if "(" in data[0]:
                        structured_pathway = True

                line = file_handle.readline()

            file_handle.close()

            # If this is a structured pathways set, then store the structure.
            if structured_pathway:
                reactions = self._set_pathways_structure(reactions,
                                                         reaction_names)

            self._store_pathways(reactions)

    def is_structured(self):
        '''Return True if this is a set of structured pathways.'''

        if self.__pathways_structure:
            return True
        else:
            return False

    def add_pathway_structure(self, pathway, structure,
                              reactions_database=None):
        '''Add the string structure for a pathway.'''

        reaction_names = None
        if reactions_database is not None:
            reaction_names = reactions_database.reaction_list()

        reactions = self._set_pathways_structure({pathway: structure},
                                                 reaction_names)
        self._store_pathways(reactions)

    def add_pathway(self, pathway, reactions):
        '''Add the unstructured pathway.'''

        self._store_pathways({pathway: reactions})

    def get_structure_for_pathway(self, pathway):
        '''Return the structure for a pathway.'''

        return copy.deepcopy(self.__pathways_structure.get(pathway, []))

    def get_key_reactions_for_pathway(self, pathway):
        '''Return the key reactions for a pathway.'''

        return copy.copy(self.__key_reactions.get(pathway, []))

    def find_reactions(self, pathway):
        '''Return the list of reactions associated with the pathway.'''

        return copy.copy(self.__pathways_to_reactions.get(pathway, []))

    def find_pathways(self, reaction):
        '''Return the list of pathways associated with the reaction.'''

        return copy.copy(self.__reactions_to_pathways.get(reaction, []))

    def reaction_list(self):
        '''Return the list of reactions included in the database.'''

        return list(self.__reactions_to_pathways.keys())

    def pathway_list(self):
        '''Return the list of pathways included in the database.'''

        return list(self.__pathways_to_reactions.keys())

    def get_database(self, min_reactions=1):
        '''Return the database as a flat file with a single pathway per
        line.'''

        data = []
        for pathway in sorted(self.__pathways_to_reactions):

            reactions = self.__pathways_to_reactions[pathway]

            if len(reactions) < min_reactions:
                continue
            else:
                data.append(pathway + " " + " ".join(reactions))

        return "\n".join(data)


def pathway_pipeline(inputfile,
                     mapfile,
                     out_dir,
                     proc=1,
                     run_minpath=True,
                     coverage=False,
                     no_regroup=False,
                     regroup_mapfile=None,
                     gap_fill_on=True,
                     per_sequence_contrib=False,
                     per_sequence_abun=None,
                     per_sequence_function=None,
                     wide_table=False,
                     print_cmds=False):
    '''Pipeline containing full pipeline for reading input files, making
    calls to functions to run MinPath and calculate pathway abundances and
    coverages. Will return: (1) unstratified pathway abundances, (2)
    unstratified pathway coverages, (3) stratified pathway abundances, (4)
    stratified pathway coverages, (5) pathway abundance predictions per
    sequence, (6) pathway coverage predictions per sequence, (7) unstratified
    pathway abundances based on per-sequence predictions. An object of class
    None will be returned for any non-applicable value.'''

    # If no regrouping flag set then set input regrouping mapfile to be None.
    if no_regroup:
        regroup_mapfile = None

    # Read in table of gene family abundances and determine if in unstratified,
    # stratified, or contribution format.
    in_metagenome, in_format = read_metagenome_input(inputfile)

    # Basic checks if --per_sequence_contrib set.
    if per_sequence_contrib:

        # Throw error if --per_sequence_contrib set, but --per_sequence_abun
        # and/or --per_sequence_function not set.
        if not per_sequence_abun or not per_sequence_function:
            sys.exit("Error: \"--per_sequence_contrib\" option set, but at "
                     "least one of \"per_sequence_abun\" or "
                     "\"--per_sequence_function\" were not set. These input "
                     "arguments need to be specified when "
                     "\"--per_sequence_contrib\" is used.")

        check_files_exist([per_sequence_abun, per_sequence_function])

    # Throw error file-format and wide table setting not compatible.
    if in_format == "strat" and not wide_table:
        sys.exit("Error: stratified table input (deprecated format), but "
                 "\"--wide_table\" option not set. You should input either an "
                 "unstratified or contributional table if you do not require "
                 "a wide-format table.")

    if in_format == "contrib" and wide_table and not per_sequence_contrib:
        sys.exit("Error: contributional table input, but \"--wide_table\" "
                 "option set. This option specifies that deprecated "
                 "wide-format stratified tables should be output, which "
                 "is only allowed when a wide-format stratified table is "
                 "input or the --per_sequence_contrib option is set.")

    # Remove 'description' column if it exists.
    if "description" in in_metagenome.columns:
        in_metagenome.drop("description", axis=1, inplace=True)

    # Get list of sample ids.
    if in_format == "contrib":
        samples = in_metagenome['sample'].unique()
    else:
        samples = [col for col in in_metagenome.columns
                   if col not in ["function", "sequence"]]

    # Initialize reactions to be empty unless regroup mapfile given.
    reactions = []

    # Regroup functions in input table to different ids if regroup mapfile is
    # provided.
    if regroup_mapfile:
        reactions = read_reaction_names(regroup_mapfile)

        in_metagenome = regroup_func_ids(in_metagenome, in_format,
                                         regroup_mapfile, proc)
        regrouped_outfile = path.join(out_dir, "regrouped_infile.tsv")
        in_metagenome.to_csv(path_or_buf=regrouped_outfile, sep="\t",
                             index=False)

    # Read in pathway structures.
    pathways_in = PathwaysDatabase(database=mapfile, reaction_names=reactions)

    # Write out mapfile with all structure removed.
    if run_minpath:
        minpath_mapfile = path.join(out_dir, "parsed_mapfile.tsv")
        with open(minpath_mapfile, "w") as out_map:
            out_map.write(pathways_in.get_database())
    else:
        minpath_mapfile = None

    # Subset input table of reactions to only those found in pathway database.
    in_metagenome = in_metagenome[in_metagenome.function.isin(pathways_in.reaction_list())]

    # Initialize output objects to be None (expect for unstratified abundance).
    path_cov_unstrat = None
    path_cov_strat = None
    path_abun_strat = None
    path_cov_by_seq = None
    path_abun_by_seq = None
    path_abun_unstrat_by_seq = None

    minpath_out_dir = path.join(out_dir, "minpath_running")
    make_output_dir(minpath_out_dir)

    if in_format == "contrib":
        # Get unstratified and stratified pathway levels.
        # Note that stratified tables will only be returned by this step (and
        # the "strat" option below) if per_sequence_contrib=False (extra step
        # required below).
        path_out_raw = Parallel(n_jobs=proc)(delayed(contrib_format_pathway_levels)(
                                                     sample_id,
                                                     in_metagenome.loc[in_metagenome['sample'] == sample_id],
                                                     minpath_mapfile,
                                                     minpath_out_dir,
                                                     pathways_in,
                                                     run_minpath,
                                                     coverage,
                                                     gap_fill_on,
                                                     per_sequence_contrib,
                                                     print_cmds)
                                                     for sample_id in samples)

    elif in_format == "strat":

        path_out_raw = Parallel(n_jobs=proc)(delayed(basic_strat_pathway_levels)(
                                                     sample_id,
                                                     in_metagenome[["function", "sequence", sample_id]],
                                                     minpath_mapfile,
                                                     minpath_out_dir,
                                                     pathways_in,
                                                     run_minpath,
                                                     coverage,
                                                     gap_fill_on,
                                                     per_sequence_contrib,
                                                     print_cmds)
                                                     for sample_id in samples)

    # Otherwise the data is in unstratified format, which is more straight-
    # forward to process.
    else:
        path_out_raw = Parallel(n_jobs=proc)(delayed(
                                               unstrat_pathway_levels)(
                                                   sample_id,
                                                   in_metagenome[["function", sample_id]],
                                                   minpath_mapfile,
                                                   minpath_out_dir,
                                                   pathways_in,
                                                   run_minpath,
                                                   coverage,
                                                   gap_fill_on,
                                                   print_cmds)
                                               for sample_id in samples)

    # Prep output unstratified DataFrames.
    path_raw_abun_unstrat = []
    path_raw_cov_unstrat = []

    for sample_output in path_out_raw:
        path_raw_abun_unstrat += [sample_output[0]]
        path_raw_cov_unstrat += [sample_output[1]]

    path_abun_unstrat = prep_pathway_df_out(path_raw_abun_unstrat)
    path_abun_unstrat.columns = samples
    path_abun_unstrat.sort_index(axis=0, inplace=True)

    if coverage:
        path_cov_unstrat = prep_pathway_df_out(path_raw_cov_unstrat,
                                               num_digits=10)
        path_cov_unstrat.columns = samples
        path_cov_unstrat.sort_index(axis=0, inplace=True)
    else:
        path_cov_unstrat = None

    # If --per_sequence_contrib not set then prep output stratified
    # table the same as the unstratified tables.
    if not per_sequence_contrib and in_format != "unstrat":
        path_raw_abun_strat = []

        for sample_output in path_out_raw:
            path_raw_abun_strat += [sample_output[2]]

        if in_format == "strat":
            path_abun_strat = prep_pathway_df_out(path_raw_abun_strat,
                                                  strat_index=True)
            path_abun_strat.columns = ["pathway", "sequence"] + samples
            path_abun_strat.sort_values(['pathway', 'sequence'], inplace=True)

        elif in_format == "contrib":
            path_abun_strat = pd.concat(path_raw_abun_strat)
            path_abun_strat.sort_values(['sample', 'function', 'taxon'],
                                        inplace=True)

    # Calculate pathway levels for each individual sequence (in parallel)
    # and then multiply this table by the abundance of each sequence
    # within each sample (using same approach as in metagenome pipeline).
    if per_sequence_contrib:

        per_seq_out_dir = path.join(out_dir, "minpath_running_per_seq")
        make_output_dir(per_seq_out_dir)

        path_abun_strat, \
        path_cov_strat, \
        path_abun_by_seq, \
        path_cov_by_seq = per_sequence_contrib_levels(sequence_abun=per_sequence_abun,
                                                      sequence_func=per_sequence_function,
                                                      minpath_map=minpath_mapfile,
                                                      per_seq_out_dir=per_seq_out_dir,
                                                      pathway_db=pathways_in,
                                                      run_minpath=run_minpath,
                                                      calc_coverage=coverage,
                                                      gap_fill_on=gap_fill_on,
                                                      nproc=proc,
                                                      regroup_map=regroup_mapfile,
                                                      wide_table=wide_table,
                                                      print_opt=print_cmds)

        if wide_table:
            path_abun_unstrat_by_seq = strat_to_unstrat_counts(strat_df=path_abun_strat,
                                                               func_col="pathway")
        else:
            path_abun_unstrat_by_seq = contrib_to_unstrat(contrib_table=path_abun_strat,
                                                          sample_order=list(path_abun_unstrat.columns.values))

    return(path_abun_unstrat, path_cov_unstrat, path_abun_strat,
           path_cov_strat, path_abun_by_seq, path_cov_by_seq,
           path_abun_unstrat_by_seq)


def prep_pathway_df_out(in_tab, strat_index=False, num_digits=4):
    '''Takes in lists of pathway abundances or coverages in series format and
    converts them into pandas dataframe to be output. If index labels are in
    stratified format then convert them to columns.'''

    # Convert these returned lists of series into pandas dataframes.
    in_tab_df = pd.concat(in_tab, axis=1, sort=True)

    # Replace all missing values (NaN)
    in_tab_df = in_tab_df.fillna(0)

    # Remove rows with all missing values.
    in_tab_df = in_tab_df.loc[~(in_tab_df == 0).all(axis=1)]

    if strat_index:
        orig_col = list(in_tab_df.columns)

        if not in_tab_df.empty:
            # Split stratified index into 2 new columns.
            in_tab_df['pathway'], in_tab_df['sequence'] = in_tab_df.index.str.split('\\|\\|\\|', 1).str

        else:
            in_tab_df['pathway'], in_tab_df['sequence'] = "", ""

        # Add these columns to be first.
        in_tab_df = in_tab_df[['pathway', 'sequence'] + orig_col]

    # Round to 4 decimals and return.
    return(in_tab_df.round(decimals=num_digits))


def read_metagenome_input(filename):
    '''Reads in gene abundance table which can be either unstratified,
    stratified, or contributional format (outputs of metagenome_pipeline.py).
    Will return the Pandas dataframe of this table and one of "unstrat",
    "strat" or "contrib" to indicate the format.'''

    # Read in input file as pandas dataframe.
    input_df = pd.read_csv(filename, sep="\t", dtype={'function': str})

    # Check that 'function' column is present, which is in all three input
    # formats.
    if 'function' not in input_df.columns:
        sys.exit("Error: required column \"function\" not found in input " +
                 "metagenome table")

    # Determine format based on presence of specific columns.
    if 'genome_function_count' in input_df.columns:
        table_type = "contrib"
        input_df['sample'] = input_df['sample'].astype('str')
        input_df['taxon'] = input_df['taxon'].astype('str')
    elif 'sequence' in input_df.columns:
        table_type = "strat"
        input_df['sequence'] = input_df['sequence'].astype('str')
    else:
        table_type = "unstrat"

    return(input_df, table_type)


def identify_minpath_present(report_file):
    '''Parse MinPath report output file and returns set containing all pathways
    ids that were called as present.'''

    path_present = set()

    with open(report_file, "r") as minpath_report_in:
        for line in minpath_report_in:
            line_split = line.split()

            if int(line_split[7]) == 1:
                path_present.add(line_split[-1])

    return(path_present)


def path_abun_weighted_by_seq(reaction_abun, func_ids, total_sum, path_abun,
                              pathway):
    '''Takes in a stratified dataframe, subsets functions to those of interest,
    and pivots by sequence column (takes sum over other columns). Also takes
    in total sum of gene families that went into calculating pathway abundance
    and the calculated pathway abundance. Will return the weighted pathway
    abundance contributed by each sequence as a Pandas Series.'''

    # Subset to genes in pathway.
    reaction_abun = reaction_abun.loc[reaction_abun['function'].isin(func_ids)]

    # Drop function column.
    reaction_abun = reaction_abun.drop(["function"], axis=1)

    # Get dataframe with sum of all genes per sequence.
    seq_path_abun = pd.pivot_table(data=reaction_abun, index="sequence",
                                   aggfunc=np.sum)

    # Get weighted pathway abundance (rounded).
    strat_path_abun = np.around((seq_path_abun / total_sum) * path_abun,
                                decimals=4)

    # Remove rows that are all 0.
    strat_path_abun = strat_path_abun.loc[strat_path_abun[strat_path_abun.columns[0]] > 0, :]

    # Rename index labels to be "pathway|||sequence".
    strat_path_abun.index = ["|||".join([pathway, str(seq)]) for seq in strat_path_abun.index]

    # Return as series.
    return(strat_path_abun.iloc[:, 0])


def contributional_path_abun(reaction_abun, func_ids, total_sum, path_abun,
                             pathway):
    '''Takes in a contributional table and subsets functions involved in
    pathway of interest. Then pivot by taxon and taxon abundance columns. Also
    takes in total sum of gene families that went into calculating pathway
    abundance and the calculated pathway abundance. Will return the weighted
    pathway abundance contributed by each sequence as a Pandas Series.'''

    # Subset to genes in pathway.
    reaction_abun = reaction_abun.loc[reaction_abun['function'].isin(func_ids)]

    # Drop function column.
    reaction_abun = reaction_abun.drop(["function"], axis=1)

    # Get dataframe with sum of all genes per sequence.
    contrib_path = pd.pivot_table(data=reaction_abun,
                                  index=["taxon", "taxon_abun",
                                          "taxon_rel_abun"],
                                  aggfunc=np.sum)

    contrib_path.reset_index(inplace=True)

    contrib_path.columns.name = None

    # Get weighted pathway abundance (rounded).
    contrib_path['taxon_function_abun'] = np.around((contrib_path['taxon_function_abun'] / total_sum) * path_abun,
                                                    decimals=4)

    contrib_path['taxon_rel_function_abun'] = np.around((contrib_path['taxon_rel_function_abun'] / total_sum) * path_abun,
                                                        decimals=4)

    # Remove rows that are all 0.
    contrib_path = contrib_path.loc[contrib_path['taxon_function_abun'] > 0, :]

    contrib_path["genome_function_count"] = contrib_path['taxon_function_abun'] / contrib_path["taxon_abun"]

    contrib_path["function"] = pathway

    return(contrib_path)


def minpath_wrapper(sample_id, unstrat_input, minpath_map, minpath_outdir,
                    print_opt=False, extra_str=""):
    '''Run MinPath based on gene abundances in a single sample. Will return
    a set of all pathways called as present.'''

    # Define MinPath input and output filenames.
    minpath_in = path.join(minpath_outdir, str(sample_id) + extra_str +
                           "_minpath_in.txt")

    minpath_report = path.join(minpath_outdir, str(sample_id) +
                               extra_str + "_minpath_report.txt")

    minpath_details = path.join(minpath_outdir, str(sample_id) +
                                extra_str + "_minpath_details.txt")

    minpath_mps = path.join(minpath_outdir, str(sample_id) + extra_str +
                            "_minpath.mps")

    id_minpath_fh = open(minpath_in, "w")

    # Loop over all reactions (which are the index labels in unstrat table
    # unless regrouped).
    for reaction_id in unstrat_input.index.values:
        # Get count of each sequence in sample and write that sequence out
        # along with count if non-zero abundance.
        reaction_count = unstrat_input.loc[reaction_id, sample_id]

        # If 0 then skip.
        if reaction_count == 0:
            continue

        id_minpath_fh.write(reaction_id + "\t" + str(reaction_count) + "\n")

    id_minpath_fh.close()

    # Run MinPath on this sample.
    path2minpath = path.join(path.dirname(path.abspath(__file__)), 'MinPath',
                             'MinPath12hmp.py')

    minpath_cmd = path2minpath + " -any " + minpath_in + " -map " +\
                  minpath_map + " -report " + minpath_report +\
                  " -details " + minpath_details + " -mps " + minpath_mps

    system_call_check(minpath_cmd, print_out=print_opt)

    # Read through MinPath report and keep track of pathways identified
    # to be present.
    path_present = identify_minpath_present(minpath_report)

    # Return list of which pathways are present.
    return(path_present)


def per_sequence_contrib_levels(sequence_abun, sequence_func,
                                minpath_map, per_seq_out_dir, pathway_db,
                                run_minpath, calc_coverage, gap_fill_on, nproc,
                                regroup_map, wide_table, print_opt=False):
    '''Reads in sequence abundance table (e.g. BIOM), predicted genomes per
    sequence, and mappings for linking functions to pathway abundances. Will
    first infer the pathway abundances for each predicted genome. Will then
    generate stratified output (either contribution format, i.e. long, or wide
    format). Will also output coverage tables if option specified. Returns
    the predicted pathways per genome and the stratified output.'''

    # Read sequence abundance table.
    study_seq_counts = read_seqabun(sequence_abun)

    if study_seq_counts.index.name != "normalized":
        print("Input abundance table is not normalized by marker gene copy "
              "number. Be sure to point to the normalized table instead if "
              "this was not intentional.")

    # Set index name to be sequence for rest of pipeline.
    study_seq_counts.index.name = "sequence"

    # Read in predicted function abundances by sequence.
    pred_function = pd.read_csv(sequence_func, sep="\t", index_col="sequence",
                                dtype={'sequence': str})

    # Drop metadata_NSTI column if it in this table.
    if "metadata_NSTI" in pred_function.columns:
        pred_function.drop(columns="metadata_NSTI", inplace=True)

    # Transpose table so that sequence ids are columns.
    pred_function = pred_function.transpose()

    pred_function.index.name = "function"

    # Restrict to sequences overlapping in sequence table as well (which may
    # have been removed if they were above the NSTI cut-off.
    overlapping_sequence_ids = [col for col in pred_function.columns if col in study_seq_counts.index.values]

    pred_function = pred_function[overlapping_sequence_ids]

    pred_function.reset_index(inplace=True)

    # Regroup function table to be reactions if specified.
    if regroup_map:
        pred_function = regroup_func_ids(pred_function, "unstrat", regroup_map,
                                         nproc)

    # Get pathway levels for each sequence (note that sample info is not used
    # here).
    per_seq_raw_out = Parallel(n_jobs=nproc)(delayed(unstrat_pathway_levels)(
                                               sequence,
                                               pred_function[["function", sequence]],
                                               minpath_map,
                                               per_seq_out_dir,
                                               pathway_db,
                                               run_minpath,
                                               calc_coverage,
                                               gap_fill_on,
                                               print_opt)
                                               for sequence in overlapping_sequence_ids)

    # Create dataframes from these outputted lists (series per sequence).
    # Prep output df. Then get stratified table with sample as columns
    # and multiply by sequence abundance per sample (for abundance).
    raw_abun = []
    raw_cov = []

    for seq_output in per_seq_raw_out:
        raw_abun += [seq_output[0]]
        raw_cov += [seq_output[1]]

    path_abun_by_seq = prep_pathway_df_out(raw_abun)
    path_abun_by_seq.columns = overlapping_sequence_ids

    path_abun_by_seq = path_abun_by_seq.transpose()

    # Subset and order study sequence table to be same as in stratified table.
    study_seq_counts = study_seq_counts.loc[overlapping_sequence_ids]

    if wide_table:
        strat_abun = strat_funcs_by_samples(func_abun=path_abun_by_seq,
                                            sample_abun=study_seq_counts,
                                            rare_seqs=[],
                                            return_unstrat=False)

        strat_abun.index.set_names("pathway", level=0, inplace=True)
        strat_abun.reset_index(drop=False, inplace=True)
        strat_abun.sort_values(['pathway', 'sequence'], inplace=True)
    else:
        strat_abun = metagenome_contributions(func_abun=path_abun_by_seq,
                                              sample_abun=study_seq_counts,
                                              rare_seqs=[])

    if calc_coverage:

        path_cov_by_seq = prep_pathway_df_out(raw_cov,
                                                  num_digits=10)
        path_cov_by_seq.columns = overlapping_sequence_ids

        # Convert study sequence abundances to be binary 1 and 0 for present
        # and absent rather than multiplying the coverages by abundances.
        study_seq_counts[study_seq_counts != 0] = 1

        path_cov_by_seq = path_cov_by_seq.transpose()

        if wide_table:

            strat_cov = strat_funcs_by_samples(func_abun=path_cov_by_seq,
                                               sample_abun=study_seq_counts,
                                               rare_seqs=[],
                                               return_unstrat=False)
            strat_cov.index.set_names("pathway", level=0, inplace=True)
            strat_cov.reset_index(drop=False, inplace=True)
            strat_cov.sort_values(['pathway', 'sequence'], inplace=True)
        else:
            strat_cov = metagenome_contributions(func_abun=path_cov_by_seq,
                                                 sample_abun=study_seq_counts,
                                                 rare_seqs=[],
                                                 skip_abun=True)
    else:
        path_cov_by_seq = None
        strat_cov = None

    return(strat_abun, strat_cov, path_abun_by_seq, path_cov_by_seq)


def contrib_format_pathway_levels(sample_id, contrib_input, minpath_map, out_dir,
                                  pathway_db, run_minpath, calc_coverage,
                                  gap_fill_on, per_sequence_contrib,
                                  print_opt=False):
    '''Read in sample_id, gene family abundance table, and out_dir, and run
    MinPath based on the gene family abundances. Returns both unstratified and
    contributional pathway abundances as dictionaries in a list when
    per_sequence_contrib=False. In this case will compute the simplistic
    "community-wide contributions" for contributional output. When
    per_sequence_contrib=True only the contributional values returned will be
    None (the stratified output will be calculated by a different step).
    Pathway coverages will also be returned when calc_coverage=True
    (unstratified only).'''

    unstrat_input = contrib_to_unstrat(contrib_input)

    # Define dictionary for keeping track of reaction abundances.
    reaction_abun = unstrat_input[sample_id].to_dict(defaultdict(int))

    if run_minpath:
        pathways_present = minpath_wrapper(sample_id, unstrat_input,
                                           minpath_map, out_dir,
                                           print_opt)
    else:
        pathways_present = set(pathway_db.pathway_list())

    # Initialize series and dataframe that will contain pathway abundances and
    # coverage scores.
    unstrat_abun = pd.Series()
    unstrat_cov = pd.Series()
    strat_abun = pd.DataFrame()

    # Return empty series if no pathways are present.
    if len(pathways_present) == 0:
        return([unstrat_abun, unstrat_cov, pd.Series(), pd.Series()])

    # Get median reaction/gene family abundance for sample, which is used for
    # calculating coverage.
    if calc_coverage:
        median_abun = calc_median_reaction_abun(reaction_abun,
                                                pathways_present, pathway_db)
    else:
        median_abun = None

    # Loop through all pathways present and get abundance and coverage.
    pathway_i = 0
    for pathway in pathways_present:

        # Get ALL reactions in pathway (which could include optional ones).
        reactions = pathway_db.find_reactions(pathway)

        # Get abundances of all of these reactions.
        path_reaction_abun = {reaction_id: reaction_abun[reaction_id] \
                              for reaction_id in reactions}

        # Get pathway abundance and coverage
        pathway_abun, pathway_cov = pathway_abun_and_coverage(pathway,
                                                              pathway_db,
                                                              path_reaction_abun,
                                                              median_abun,
                                                              calc_coverage,
                                                              gap_fill_on)

        if pathway_abun == 0:
            continue

        # Add these values to each respective pandas Series.
        unstrat_abun[pathway] = pathway_abun
        unstrat_cov[pathway] = pathway_cov

        if not per_sequence_contrib:
            # If --per_sequence_contrib not set then get stratified pathway
            # abundances simply by weighting community-wide pathway abundances
            # by the abundances of all the predicted abundances of reactions in
            # these pathways contributed by each sequence (i.e. predicted
            # genome).

            strat_path_abun = contributional_path_abun(contrib_input,
                                                       reactions,
                                                       sum(list(path_reaction_abun.values())),
                                                       unstrat_abun[pathway],
                                                       pathway)

            if pathway_i == 0:
                strat_abun = strat_path_abun
                pathway_i += 1
            else:
                strat_abun = pd.concat([strat_abun, strat_path_abun],
                                       sort=False)

            strat_abun["sample"] = sample_id

            strat_abun = strat_abun[['sample', 'function', 'taxon',
                                     'taxon_abun', 'taxon_rel_abun',
                                     'genome_function_count',
                                     'taxon_function_abun',
                                     'taxon_rel_function_abun']]

    # Return unstratified and stratified abundances and coverage scores.
    return([unstrat_abun, unstrat_cov, strat_abun])


def basic_strat_pathway_levels(sample_id, strat_input, minpath_map, out_dir,
                               pathway_db, run_minpath, calc_coverage,
                               gap_fill_on, per_sequence_contrib,
                               print_opt=False):
    '''Read in sample_id, gene family table, and out_dir, and run MinPath based
    on the gene family abundances. Returns both unstratified and stratified
    pathway abundances as dictionaries in a list when
    per_sequence_contrib=False. In this case will compute the simplistic
    "community-wide contributions" for stratified output. When
    per_sequence_contrib=True only the stratified values returned will be None
    (the stratified output will be calculated by a different step). Pathway
    coverages will also be returned when calc_coverage=True (unstratified
    only).'''

    # Get gene family abundances summed over all sequences for this sample.
    unstrat_input = strat_to_unstrat_counts(strat_input)

    # Define dictionary for keeping track of reaction abundances.
    reaction_abun = unstrat_input[sample_id].to_dict(defaultdict(int))

    if run_minpath:

        pathways_present = minpath_wrapper(sample_id, unstrat_input,
                                           minpath_map, out_dir,
                                           print_opt)
    else:
        pathways_present = set(pathway_db.pathway_list())

    # Initialize series and dataframe that will contain pathway abundances and
    # coverage scores.
    unstrat_abun = pd.Series()
    unstrat_cov = pd.Series()
    strat_abun = pd.Series()

    # Return empty series if no pathways are present.
    if len(pathways_present) == 0:
        return([unstrat_abun, unstrat_cov, pd.Series(), pd.Series()])

    # Get median reaction/gene family abundance for sample, which is used for
    # calculating coverage.
    if calc_coverage:
        median_abun = calc_median_reaction_abun(reaction_abun, pathways_present,
                                                pathway_db)
    else:
        median_abun = None

    # Loop through all pathways present and get abundance and coverage.
    for pathway in pathways_present:

        # Get ALL reactions in pathway (which could include optional ones).
        reactions = pathway_db.find_reactions(pathway)

        # Get abundances of all of these reactions.
        path_reaction_abun = {reaction_id: reaction_abun[reaction_id] for reaction_id in reactions}

        # Get pathway abundance and coverage
        pathway_abun, pathway_cov = pathway_abun_and_coverage(pathway,
                                                              pathway_db,
                                                              path_reaction_abun,
                                                              median_abun,
                                                              calc_coverage,
                                                              gap_fill_on)

        if pathway_abun == 0:
            continue

        # Add these values to each respective pandas Series.
        unstrat_abun[pathway] = pathway_abun
        unstrat_cov[pathway] = pathway_cov

        if not per_sequence_contrib:
            # If --per_sequence_contrib not set then get stratified pathway
            # abundances simply by weighting community-wide pathway abundances
            # by the abundances of all the predicted abundances of reactions in
            # these pathways contributed by each sequence (i.e. predicted
            # genome).

            strat_path_abun = path_abun_weighted_by_seq(strat_input,
                                                        reactions,
                                                        sum(list(path_reaction_abun.values())),
                                                        unstrat_abun[pathway],
                                                        pathway)

            strat_abun = pd.concat([strat_abun, strat_path_abun], sort=True)

    # Return unstratified and stratified abundances and coverage scores.
    return([unstrat_abun, unstrat_cov, strat_abun])


def unstrat_pathway_levels(sample_id, unstrat_input, minpath_map, out_dir,
                           pathway_db, run_minpath, calc_coverage,
                           gap_fill_on, print_opt=False, extra_str=""):
    '''Read in sample_id, gene family table, and out_dir. Returns unstratified
    pathway abundances as dictionaries in a list. Also returns the coverage of
    each unstratified pathway as the a different dictionary in a list (if
    calc_coverage=True)'''

    unstrat_input.set_index("function", inplace=True)

    # Define dictionary for keeping track of reaction abundances.
    reaction_abun = unstrat_input[sample_id].to_dict(defaultdict(int))

    if run_minpath:

        pathways_present = minpath_wrapper(sample_id, unstrat_input,
                                           minpath_map, out_dir,
                                           print_opt)
    else:
        pathways_present = set(pathway_db.pathway_list())

    # Initialize series that will contain pathway abundances and coverage.
    unstrat_abun = pd.Series([])
    unstrat_cov = pd.Series([])

    # Return empty series if no pathways are present.
    if len(pathways_present) == 0:
        return([unstrat_abun, unstrat_cov])

    # Get median reaction/gene family abundance for sample, which is used for
    # calculating coverage.
    if calc_coverage:
        median_abun = calc_median_reaction_abun(reaction_abun, pathways_present,
                                                pathway_db)
    else:
        median_abun = None

    # Loop through all pathways present and get abundance and coverage.
    for pathway in pathways_present:
        # Note that here the term "reactions" is used interchangeably with
        # gene families (reactions is the term used by HUMAnN2) - the gene
        # families above may have already been converted to reactions or not.

        # Get ALL reactions in pathway (which could include optional ones).
        reactions = pathway_db.find_reactions(pathway)

        # Get abundances of all of these reactions.
        path_reaction_abun = {reaction_id: reaction_abun[reaction_id] for reaction_id in reactions}

        # Get pathway abundance and coverage
        pathway_abun, pathway_cov = pathway_abun_and_coverage(pathway,
                                                              pathway_db,
                                                              path_reaction_abun,
                                                              median_abun,
                                                              calc_coverage,
                                                              gap_fill_on)

        # Add these values to each respective pandas Series.
        unstrat_abun[pathway] = pathway_abun
        unstrat_cov[pathway] = pathway_cov

    return([unstrat_abun, unstrat_cov])


def pathway_abun_and_coverage(pathway, pathway_db, reaction_abun, median_value,
                              calc_coverage, gap_fill_on):
    '''Determine pathway abundance and coverage for either structured or
    unstructured pathway. Calculating coverage is off by default.'''

    pathway_cov = None

    # Check if pathway database is structured. If so then get structure and
    # use it and the key reactions to get pathway abundance and coverage.
    if pathway_db.is_structured():
        structure = pathway_db.get_structure_for_pathway(pathway)
        key_reactions = pathway_db.get_key_reactions_for_pathway(pathway)

        # Run gap filling if set.
        if gap_fill_on:
            reaction_abun = gap_fill(key_reactions, reaction_abun)

        # Calculate pathway abundance.
        pathway_abun = compute_structured_pathway_abundance_or_coverage(structure,
                                                                        key_reactions,
                                                                        reaction_abun,
                                                                        False,
                                                                        median_value)
        if calc_coverage:
            # Calculate pathway coverage.
            pathway_cov = compute_structured_pathway_abundance_or_coverage(structure,
                                                                           key_reactions,
                                                                           reaction_abun,
                                                                           True,
                                                                           median_value)

    else:
        # Otherwise, sort enzyme reactions, take second half, and get their
        # mean abundance.

        # First get indices of sorted list.
        reaction_abun_only = list(reaction_abun.values())

        sorted_index = list(np.argsort(reaction_abun_only))
        sorted_reaction_abun = [reaction_abun_only[i] for i in sorted_index]

        # Take second half of gene family abundances.
        half_i = int(len(sorted_reaction_abun) / 2)
        subset_reaction_abun = sorted_reaction_abun[half_i:]

        # Take mean for unstratified pathway abundance.
        pathway_abun = sum(subset_reaction_abun) / len(subset_reaction_abun)

        if calc_coverage:
            # Get coverage of unstructured pathway by getting the proportion
            # of reactions with scores greater than the median.
            count_higher_than_median = 0
            for abun in reaction_abun_only:
                if abun > median_value:
                    count_higher_than_median += 1

            pathway_cov = count_higher_than_median / len(reaction_abun_only)

    return(pathway_abun, pathway_cov)


def regroup_func_ids(in_df, in_format, mapfile, proc):
    '''Reads in df and mapfile of function ids to new ids. Will
    return df with new functions ids under column "function".'''

    # Read throw id mapfile and get links between two id sets.
    func_map = defaultdict(list)

    with open(mapfile, 'r') as in_map:
        for line in in_map:
            line = line.rstrip()
            line_split = line.split()
            func_map[line_split[0]] += line_split[1].split(",")

    # Get set of all unique functions.
    functions = list(set(in_df['function']))

    chunk_size = int(len(functions) / proc) + 1

    function_chunks = [functions[i:i + chunk_size]
                       for i in range(0, len(functions), chunk_size)]

    # Loop over all functions in parallel and return pandas dataframe for each
    # function with regrouped ids (which will are combined together afterwards).
    raw_new_ids_dfs = Parallel(n_jobs=proc)(delayed(
                                    convert_func_ids)(func_subset,
                                                      in_df,
                                                      func_map)
                                    for func_subset in function_chunks)

    # Combine all returned DFs into a single DF.
    raw_new_ids_combined = pd.concat(raw_new_ids_dfs, sort=False)

    if in_format == "contrib":
        regrouped_table = pd.pivot_table(raw_new_ids_combined,
                                         index=['sample', 'function', 'taxon',
                                                'taxon_abun',
                                                'taxon_rel_abun'],
                                         aggfunc=np.sum)

    elif in_format == "strat":
        regrouped_table = pd.pivot_table(raw_new_ids_combined,
                                         index=["function", "sequence"],
                                         aggfunc=np.sum)
    elif in_format == "unstrat":
        regrouped_table = pd.pivot_table(raw_new_ids_combined,
                                         index="function",
                                         aggfunc=np.sum)

    regrouped_table.reset_index(inplace=True)
    regrouped_table.columns.name = None

    # If no rows remain then throw error.
    if regrouped_table.shape[0] == 0:
        sys.exit("\nError - no rows remain after regrouping input table. The "
                 "default pathway and regroup mapfiles are meant for EC "
                 "numbers. Note that KEGG pathways are not supported since "
                 "KEGG is a closed-source database, but you can input custom "
                 "pathway mapfiles if you have access. If you are using a "
                 "custom function database did you mean to set the "
                 "--no-regroup flag and/or change the default pathways "
                 "mapfile used?\n")

    return(regrouped_table)


def convert_func_ids(functions, in_df, func_map):
    '''Will return dataframe with all new ids replacing the value in
    "function" column for each input function.'''

    new_dfs = []

    # Loop over all functions.
    for func in functions:

        # Get subset of dataframe for this function and ids to regroup by.
        in_df_subset = in_df.loc[in_df['function'] == func]
        new_ids = func_map[func]

        # Inititalize empty dataframe to add new rows.
        new_df = pd.DataFrame(columns=in_df_subset.columns.values)

        # For each new id to regroup by replace all original function ids with
        # this new id and add to new df.
        for new_id in new_ids:
            new_subset = in_df_subset.copy()
            new_subset['function'] = new_id
            new_df = new_df.append(new_subset)

        new_dfs += [new_df]

    # Concatenate these new dfs together for all input functions.
    combined_new_df = pd.concat(new_dfs, sort=False)

    return(combined_new_df)


def read_reaction_names(reactions_database):
    '''Read in the reactions from a table that contains links between reactions
    and gene family ids. Will return a list of reactions (which are assumed to
    be the first field of the file after splitting by " ").'''

    # Check that the input file exists.
    check_files_exist([reactions_database])

    reactions = []

    with open(reactions_database, "rt") as infile:
        for line in infile:
            line_split = line.strip().split(" ")
            if len(line_split) > 1:
                reactions += line_split[0]

    return reactions


def compute_structured_pathway_abundance_or_coverage(structure,
                                                     key_reactions,
                                                     reaction_abun,
                                                     calc_coverage,
                                                     median_value):

    '''Returns the abundance or the coverage of a structured pathway based on
    the pathway structure, list of key reactions, and abundance of reactions in
    the pathway. Coverage will be calculated when calc_coverage=True, which
    uses the median abundance of reactions in the pathway.
    This calculation is based on the approach implemented in HUMAnN2 and below
    code was modified from the compute_structured_pathway_abundance_or_coverage
    function in HUMAnN2 v0.11.1.'''

    # Initialize lists to contain scores (abundances OR coverages) for required
    # and optional reactions separately.
    key_reaction_scores = []
    optional_reaction_scores = []

    # Process through the structure to compute the abundance.
    # Select the join instead of removing from the list to not alter the list
    # for the calling function.
    join = structure[0]
    for item in structure[1:]:
        if isinstance(item, list):
            # If the item is a list then recursively determine that part of the
            # pathway's abundance.
            key_reaction_scores.append(compute_structured_pathway_abundance_or_coverage(item,
                key_reactions, reaction_abun, calc_coverage, median_value))
        else:
            reaction_score = reaction_abun[item]

            # If this is a coverage calculation then get the value from the
            # chi-square cdf based on this abundance and median value as
            # parameters.
            if calc_coverage:
                reaction_score = scipy.stats.chi2.cdf(reaction_score,
                                                      median_value)

            # Check if this is a key reaction and append to the
            # appropriate list.
            if item in key_reactions:
                key_reaction_scores.append(reaction_score)
            else:
                optional_reaction_scores.append(reaction_score)

    # If this is an OR join then use the max of all of the reaction abundances.
    if join == ",":
        all_reaction_scores = key_reaction_scores + optional_reaction_scores
        abundance = 0

        if all_reaction_scores:
            abundance = max(all_reaction_scores)
    else:
        # If this is not an OR, then take the harmonic mean of the reactions.
        # First get the harmonic mean of the key reactions.

        abundance = harmonic_mean(key_reaction_scores)

        # Then, if there are any optional reaction abundances get the harmonic
        # mean of all reactions that have a great abundance than the harmonic
        # mean of the key reactions alone.
        if optional_reaction_scores:
            optional_reaction_scores_filt = [value for value in optional_reaction_scores if value > abundance]
            abundance = harmonic_mean(key_reaction_scores + optional_reaction_scores_filt)

    return abundance


def gap_fill(key_reactions, reaction_abun):
    '''Based on approach implemented in HUMAnN2: if at most only 1 key reaction
    has an abundance of 0 then set the reaction with lowest abundance to be
    that of the 2nd least abundant reaction.'''

    reaction_abun_gap_filled = reaction_abun.copy()

    missing_key_reaction = 0

    # Get all abundances for the key reactions.
    key_reactions_abun = []
    for reaction in key_reactions:
        abun = reaction_abun[reaction]

        if abun == 0:
            missing_key_reaction += 1

        key_reactions_abun.append(abun)

    # Return same abundances as before if more than 1 key reaction is zero.
    if missing_key_reaction > 1:
        return reaction_abun_gap_filled

    # Determine sorted order of key reaction abundances.
    key_abun_order = list(np.argsort(key_reactions_abun))

    # Reorder key reactions name and abundance lists.
    key_reactions = [key_reactions[i] for i in key_abun_order]
    key_reactions_abun = [key_reactions_abun[i] for i in key_abun_order]

    # Set key reaction with lowest abundance to have the same abundance as the
    # second least abundant key reaction.
    reaction_abun_gap_filled[key_reactions[0]] = key_reactions_abun[1]

    return reaction_abun_gap_filled


def harmonic_mean(values, decimal_points=16):
    '''Returns harmonic mean of input list of numbers. Will round to specified
    number of decimal points first. Modified from HUMAnN2 v0.11.1.'''

    hmean = 0

    # Return 0 if values is empty.
    if not values:
        return(hmean)

    values = list(np.around(values, decimals=decimal_points))

    # Return 0 if any values are 0 (since 1/0 is undefined).

    if min(values) > 0:
        reciprocal_sum = sum((1.0 / v) for v in values)
        hmean = len(values) / reciprocal_sum

    return hmean


def calc_median_reaction_abun(reaction_abun, pathways_present, pathway_db):
    '''Calculate the median reaction abundance across all pathways (reactions
    in multiple pathways are counted multiple times).'''

    all_reaction_abun = []

    # Loop over all pathways present and get the abundances for all reactions
    # within.
    for pathway in pathways_present:
        for reaction in pathway_db.find_reactions(pathway):
            abun = reaction_abun[reaction]
            if abun > 0:
                all_reaction_abun.append(abun)

    # Calculate and return median.
    return(np.median(all_reaction_abun))


def strat_to_unstrat_counts(strat_df, func_col="function"):
    '''Given a pandas dataframe with the columns "sequence", "function" (by
    default), and at least 1 sample column, will return the dataframe after
    removing sequence column and summing all functions per sample. Functions
    will be new index labels.'''

    strat_df = strat_df.drop(["sequence"], axis=1)

    return(pd.pivot_table(data=strat_df, index=func_col, aggfunc=np.sum))
