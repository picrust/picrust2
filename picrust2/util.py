#!/usr/bin/env python

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.0.0-b.5"

from os import makedirs
from os.path import abspath, dirname, isdir, join, exists
from subprocess import call
import pandas as pd
import numpy as np
import tempfile
import sys


def read_fasta(filename, cut_header=False):

    '''Read in FASTA file and return dictionary with each independent sequence
    id as a key and the corresponding sequence string as the value.
    '''

    # Intitialize empty dict.
    seq = {}

    # Intitialize undefined str variable to contain the most recently parsed
    # header name.
    name = None

    # Read in FASTA line-by-line.
    with open(filename, "r") as fasta:

        for line in fasta:

            # If header-line then split by whitespace, take the first element,
            # and define the sequence name as everything after the ">".
            if line[0] == ">":

                if cut_header:
                    name = line.split()[0][1:]
                else:
                    name = line[1:]

                name = name.rstrip("\r\n")

                # Intitialize empty sequence with this id.
                seq[name] = ""

            else:
                # Remove line terminator/newline characters.
                line = line.rstrip("\r\n")

                # Add sequence to dictionary.
                seq[name] += line

    return seq


def write_fasta(seq, outfile):
    out_fasta = open(outfile, "w")

    for s in seq:
        out_fasta.write(">" + s + "\n")
        out_fasta.write(seq[s] + "\n")

    out_fasta.close()


def read_phylip(filename, check_input=True):

    '''Reads in Phylip formatted multiple-sequence alignment and stores in
    dictionary with each key as a sequence name and each value as the sequence.
    If check_input=True then will check whether the number and length of each
    sequence is correctly specified (first line of file).
    '''

    # Intitialize empty dict.
    seq = {}

    line_count = 0

    # Read in phylip line-by-line.
    with open(filename, "r") as phylip:

        for line in phylip:

            line_split = line.split()

            # Header-line, store the # of seqs and length of seqs and skip.
            if line_count == 0:
                num_seqs = int(line_split[0])
                seq_length = int(line_split[1].rstrip("\r\n"))
                line_count += 1
                continue

            # Add sequence to dictionary.
            seq[line_split[0]] = line_split[1].rstrip("\r\n")

            if check_input:
                # Check length of sequence is correct.
                if len(seq[line_split[0]]) != seq_length:
                    raise SystemExit("Expected sequence " + line_split[0] +
                                     " to be length " + str(seq_length) +
                                     ", but was " +
                                     str(len(seq[line_split[0]])))

    if check_input:
        # Check number of sequences read in.
        if len(seq) != num_seqs:
            raise SystemExit("Expected " + str(num_seqs) + " sequences, but " +
                             "found " + str(len(seq)))

    return seq


def write_phylip(seq, outfile):
    out_phylip = open(outfile, "w")

    seq_count = 0

    for s in seq:

        # Write header if this is first line.
        if seq_count == 0:

            seq_length = len(seq[s])

            out_phylip.write(str(len(seq)) + " " + str(seq_length) + "\n")

            seq_count += 1

        if seq_length != len(seq[s]):

            # Throw error if sequence is unexpected length.
            raise SystemExit("Expected sequence " + s + " to be length " +
                             str(seq_length) + ", but was " +
                             str(len(seq[s])))

        out_phylip.write(s + " " + seq[s] + "\n")

    out_phylip.close()


def system_call_check(cmd, print_out=False, stdout=None, stderr=None):
    """Run system command and throw and error if return is not 0. Input command
    can be a list containing the command or a string."""

    # Print command out if option set.
    if print_out:
        if type(cmd) is list:
            print(" ".join(cmd))
        else:
            print(cmd)

    # Convert command to list if input as string.
    if type(cmd) is str:
        cmd = cmd.split()

    return_value = call(cmd, stdout=stdout, stderr=stderr)

    # Exit with error if command did not finish successfully.
    if return_value != 0:
        raise SystemExit("Error running this command:\n" + " ".join(cmd))
    else:
        return(return_value)


def get_picrust_project_dir():
    """ Returns the top-level PICRUST directory
    """
    # Get the full path of util.py
    current_file_path = abspath(__file__)
    # Get the directory containing util.py
    current_dir_path = dirname(current_file_path)
    # Return the directory containing the directory containing util.py
    return dirname(current_dir_path)


def make_output_dir(dirpath, strict=False):
    """Make an output directory if it doesn't exist

    Returns the path to the directory
    dirpath -- a string describing the path to the directory
    strict -- if True, raise an exception if dir already
    exists
    """
    dirpath = abspath(dirpath)

    # Check if directory already exists
    if isdir(dirpath):
        if strict:
            err_str = "Directory '%s' already exists" % dirpath
            raise IOError(err_str)

        return dirpath
    try:
        makedirs(dirpath)
    except IOError as e:
        err_str = "Could not create directory '%s'. Are permissions set " +\
                  "correctly? Got error: '%s'" %e
        raise IOError(err_str)

    return dirpath


def make_output_dir_for_file(filepath):
    """Create sub-directories for a new file if they don't already exist"""
    dirpath = dirname(filepath)
    if not isdir(dirpath) and not dirpath == '':
        makedirs(dirpath)


def generate_temp_filename(temp_dir=None, prefix="", suffix=""):
    '''Function to generate path to temporary filenames (does not create the)
    files). The input arguments can be used to customize the temporary
    filename. If no temporary directory is specified then the default temprary
    directory will be used.'''

    # If temp_dir not set then get default directory.
    if not temp_dir:
        temp_dir = tempfile._get_default_tempdir()

    return(join(temp_dir, prefix +
                next(tempfile._get_candidate_names()) + suffix))


def biom_to_pandas_df(biom_tab):
    '''Will convert from biom Table object to pandas dataframe.'''

    # Note this is based on James Morton's blog post:
    # http://mortonjt.blogspot.ca/2016/07/behind-scenes-with-biom-tables.html)

    return(pd.DataFrame(np.array(biom_tab.matrix_data.todense()),
                                 index=biom_tab.ids(axis='observation'),
                                 columns=biom_tab.ids(axis='sample')))


def three_df_index_overlap_sort(df1, df2, df3):
    '''Given 3 pandas dataframes, will first determine which index labels
    overlap across all dataframes and will subset the labels to this set and
    then will sort the dataframes to be in the same order'''

    label_overlap = df1.index.intersection(df2.index.intersection(df3.index))

    # If there are no overlapping labels then throw error.
    if len(label_overlap) == 0:
        raise ValueError("No sequence ids overlap between all three of the " +
                         "input files.")

    df1 = df1.reindex(label_overlap)
    df2 = df2.reindex(label_overlap)
    df3 = df3.reindex(label_overlap)

    return(df1, df2, df3)


def check_files_exist(filepaths):
    '''Takes in a list of filepaths and checks whether they exist. Will
    throw error describing which files do not exist if applicable.'''

    num_nonexist = 0

    missing_files = []

    for filepath in filepaths:

        if not exists(filepath):
            missing_files += [filepath]
            num_nonexist += 1

    if num_nonexist == 0:
        pass
    elif num_nonexist == 1:
        raise ValueError("This input file was not found: " + missing_files[0])
    elif num_nonexist > 1:
        raise ValueError("These input files were not found: " +
                         ", ".join(missing_files))


def add_descrip_col(inputfile : str, mapfile : str):
    '''Takes paths to input table and mapfile of function ids to descriptions.
    Will read both of these files in as pandas dataframes and will add
    descriptions as a separate column in a new pandas dataframe, which will be
    returned. Note that the first column of the input function abundance table
    is assumed to contain the functions ids.'''

    # Read in input tables.
    function_tab = pd.read_table(inputfile, sep="\t")
    map_tab = pd.read_table(mapfile, sep="\t", index_col=0, header=None,
                            names=["function", "description"])

    # Check to see if any of the mapfile row indices are in the function table
    # id column and throw an error if not.
    if not any(map_tab.index.isin(function_tab.iloc[:, 0])):
        sys.exit("Error: no function ids in input table are found " +
                 "in the provided mapfile")

    # Reindex mapfile to match order of functions in function table (and for
    # ids to be duplicated if the table is stratified.
    map_tab = map_tab.reindex(function_tab.iloc[:, 0], fill_value="not_found")

    function_tab.insert(1, "description", list(map_tab["description"]))

    # Add description column to function table and return.
    return(function_tab)
