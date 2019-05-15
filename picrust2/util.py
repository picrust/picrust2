#!/usr/bin/env python

__copyright__ = "Copyright 2018-2019, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.1.4-b"

from os import makedirs, chmod
from os.path import abspath, dirname, isdir, join, exists, splitext
from collections import defaultdict
from subprocess import call
import stat
import shutil as _shutil
import weakref as _weakref
import warnings as _warnings
import pandas as pd
import numpy as np
import biom
import tempfile
import gzip
import sys


def read_fasta(filename, cut_header=False):
    '''Read in FASTA file (gzipped or not) and return dictionary with each
    independent sequence id as a key and the corresponding sequence string as
    the value.'''

    # Intitialize empty dict.
    seq = {}

    # Intitialize undefined str variable to contain the most recently parsed
    # header name.
    name = None

    # Read in FASTA line-by-line.
    if filename[-3:] == ".gz":
        fasta_in = gzip.open(filename, "rt")
    else:
        fasta_in = open(filename, "r")

    for line in fasta_in:

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

    fasta_in.close()

    return seq


def write_fasta(seq, outfile):
    out_fasta = open(outfile, "w")

    # Look through sequence ids (sorted alphabetically so output file is
    # reproducible).
    for s in sorted(seq.keys()):
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
    '''Will write a dictionary containing id, sequence pairs in Phylip
    format. Originally written to run PaPaRa.'''

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


def read_stockholm(filename, clean_char=True):
    '''Reads in Stockholm formatted multiple sequence alignment and returns
    dictionary with ids as keys and full concatenated sequences as values. When
    clean_char=True this function will convert all characters to uppercase and
    convert "." characters to "-". This was originally written for converting
    hmmalign output files.'''

    # Intitialize defaultdict that will contain strings.
    seq = defaultdict(str)

    line_count = 0

    # Read in file line-by-line.
    with open(filename, "r") as stockholm:

        for line in stockholm:

            line = line.rstrip()

            # Header-line - check that it starts with "# STOCKHOLM".
            if line_count == 0 and "# STOCKHOLM" not in line:
                sys.exit("Error - stockholm format multiple-sequence "
                         "alignments should have \"# STOCKHOLM\" (and the "
                         "version number) on the first line")

            line_count += 1

            # Skip blank lines, lines that start with comment, and lines that
            # start with "//".
            if not line or line[0] == "#" or line[0:2] == "//":
                continue

            line_split = line.split()

            if clean_char:
                line_split[1] = line_split[1].upper()
                line_split[1] = line_split[1].replace(".", "-")

            # Add sequence to dictionary.
            seq[line_split[0]] += line_split[1]

    # Double-check that last line was "//"
    if line[0:2] != "//":
        sys.exit("Error - last line of stockholm file should have been "
                 "\"//\".")

    return seq


def system_call_check(cmd, print_out=False, print_stderr=False):
    '''Run system command and throw and error if return is not 0. Input command
    can be a list containing the command or a string.'''

    # Convert command to list if input as string.
    if type(cmd) is str:
        cmd = cmd.split()

    # Print command out if option set.
    if print_out:
        print(" ".join(cmd), file=sys.stderr)

    # Write stdout and stderr of command to temporary files.
    # Only output the content of these files if the job fails.
    with TemporaryDirectory() as temp_dir:

        stdout_file = join(temp_dir, "stdout.txt")
        stderr_file = join(temp_dir, "stderr.txt")

        with open(stdout_file, "wb") as stdout_fh, \
             open(stderr_file, "wb") as stderr_fh:

            return_value = call(cmd, stdout=stdout_fh,
                                stderr=stderr_fh)

        # Exit with error if command did not finish successfully.
        if return_value != 0:
            print("\nError running this command:\n" + " ".join(cmd),
                  file=sys.stderr)

            # Print out stdout and stderr.
            with open(stdout_file, 'r') as f:
                print("\nStandard output of failed command:", file=sys.stderr)
                print("\"" + f.read() + "\"", file=sys.stderr)

            with open(stderr_file, 'r') as f:
                print("\nStandard error of failed command:", file=sys.stderr)
                print("\"" + f.read() + "\"", file=sys.stderr)

            sys.exit(1)

        # Print log info to stderr is helpful for certain subprocesses.
        elif print_stderr:
            with open(stderr_file, 'r') as f:
                print(f.read(), file=sys.stderr)

    return(return_value)


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


def read_seqabun(infile):
    '''Will read in sequence abundance table in either TSV, BIOM, or mothur
    shared format.'''

    # First check extension of input file. If extension is "biom" then read in
    # as BIOM table and return. This is expected to be the most common input.
    in_name, in_ext = splitext(infile)
    if in_ext == "biom":
        return(biom.load_table(infile).to_dataframe(dense=True))

    # Next check if input file is a mothur shared file or not by read in first
    # row only.
    mothur_format = False
    try:
        in_test = pd.read_csv(filepath_or_buffer=infile, sep="\t", nrows=1)
        in_test_col = list(in_test.columns.values) 
        if len(in_test_col) >= 4 and (in_test_col[0] == "label" and \
                                      in_test_col[1] == "Group" and \
                                      in_test_col[2] == "numOtus"):
            mothur_format = True
    except Exception:
        pass

    # If identified to be mothur format then remove extra columns, set "Group"
    # to be index (i.e. row) names and then transpose.
    if mothur_format:
        input_seqabun = pd.read_csv(filepath_or_buffer=infile, sep="\t")
        input_seqabun.drop(labels=["label", "numOtus"], axis=1, inplace=True)
        input_seqabun.set_index(keys="Group", drop=True, inplace=True, 
                                verify_integrity=True)
        input_seqabun.index.name = None
        return(input_seqabun.transpose())
    else:
        return(biom.load_table(infile).to_dataframe(dense=True))


def three_df_index_overlap_sort(df1, df2, df3):
    '''Given 3 pandas dataframes, will first determine which index labels
    overlap across all dataframes and will subset the labels to this set and
    then will sort the dataframes to be in the same order'''

    label_overlap = df1.index.intersection(df2.index.intersection(df3.index)).sort_values()

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


def add_descrip_col(inputfile, mapfile, in_df=False):
    '''Takes paths to input table and mapfile of function ids to descriptions.
    Will read both of these files in as pandas dataframes and will add
    descriptions as a separate column in a new pandas dataframe, which will be
    returned. Note that the first column of the input function abundance table
    is assumed to contain the functions ids. An input dataframe rather than
    path to file can also be passed as the "inputfile", in which case in_df
    should be set to True.'''

    # Read in input tables.
    if in_df:
        function_tab = inputfile
    else:
        function_tab = pd.read_csv(inputfile, sep="\t")
    
    map_tab = pd.read_csv(mapfile, sep="\t", index_col=0, header=None,
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


def convert_humann2_to_picrust2(infiles, outfile, stratified):
    '''Reads in HUMAnN2 gene tables and will convert them to PICRUSt2
    format.'''

    humann2_samples = []

    # Loop over all sample infiles and add their data to this list.
    for infile in infiles:
        humann2_samples.append(pd.read_csv(infile, sep="\t", index_col=0))

    # Get the index name for each table and make sure they are identical.
    infile_index_names = []
    for sample_df in humann2_samples:
        infile_index_names.append(sample_df.index.name)

    infile_index_names_set = set(infile_index_names)

    if len(infile_index_names_set) > 1:
        sys.exit('Error input HUMAnN2 tables are not all for the same '
                 'datatype. The datatypes of the input files are: ' +
                 ', '.join(list(infile_index_names_set)))

    index_name = infile_index_names[0]

    # Concatenate all sample dfs into a single table.
    humann2_combined = pd.concat(humann2_samples, axis=1, sort=True)

    # Fill in zeros for missing data.
    humann2_combined = humann2_combined.fillna(0)

    # Determine name of first column of output table.
    if index_name == "# Pathway":
        first_col = "pathway"
    elif index_name == "# Gene Family":
        first_col = "function"
    else:
        sys.exit('Error: first column of input HUMAnN2 files is \"' +
                 index_name + '\". Either \"# Pathway\" or ' +
                 '\"# Gene Family\" was expected.')

    # Output table with current index if table isn't stratified.
    if not stratified:
        humann2_combined.to_csv(path_or_buf=outfile,  sep="\t",
                                index_label=first_col)

    # Otherwise if the table is stratified remove all rows that are in
    # unstratified format and split index into two columns.
    else:
        humann2_combined = humann2_combined[humann2_combined.index.str.contains('\\|')]

        original_col = list(humann2_combined.columns)

        humann2_combined[first_col], humann2_combined['sequence'] = humann2_combined.index.str.split('\\|', 1).str

        # Reorder columns.
        humann2_combined = humann2_combined.loc[:, [first_col, 'sequence'] +
                                                original_col]

        humann2_combined.to_csv(path_or_buf=outfile,  sep="\t", index=False)


def convert_picrust2_to_humann2(infiles, outfolder, stratified):
    '''Reads in a PICRUSt2 table(s) and splits each sample into different
    file as compatible with HUMAnN2.'''

    # If not stratified then read in single table.
    if not stratified:

        if len(infiles) > 1:
            sys.exit('Stopping - only expected one input file when converting '
                     'from PICRUSt2 unstratified table to HUMAnN2 format')

        in_tab = pd.read_csv(infiles[0], sep="\t", index_col=0)

        # Double-check that this table isn't stratified.
        if 'sequence' in in_tab.columns:
            sys.exit('Stopping - column named sequence was found in the input '
                     'table, but the unstratified conversion option was set.')

        index_name = in_tab.index.name

    # Otherwise read in stratified and unstratified tables.
    else:

        if len(infiles) != 2:

            sys.exit('Stopping - expected two input files (stratified and '
                     'unstratified PICRUSt2 tables) when converting to '
                     'HUMAnN2 stratified format')

        # Read in both input tables.
        in_tab1 = pd.read_csv(infiles[0], sep="\t")
        in_tab2 = pd.read_csv(infiles[1], sep="\t")

        # Make sure that only 1 input table is stratified.
        strat_table_count = 0

        if 'sequence' in in_tab1.columns:
            strat_table_count += 1

        if 'sequence' in in_tab2.columns:
            strat_table_count += 1

        if strat_table_count != 1:
            sys.exit('Stopping - exactly one stratified table should have been '
                     'input - found ' + str(strat_table_count) + '.')


        # Check that the first column labels are the same across both tables.
        if in_tab1.columns[0] != in_tab2.columns[0]:
            sys.exit('Stopping - label of first column does not match between '
                     'the two input tables: \"' + in_tab1.columns[0] +
                     '\" and \"' + in_tab2.columns[0] + '\".')

        index_name = in_tab1.columns[0]

        # Set index labels for each table appropriately and concatenate
        # together. First define convenience fucntion to set index labels
        # correctly.
        def set_picrust2_tab_index(in_table):

            if 'sequence' in in_table.columns:
                in_table.index = in_table.loc[:, in_table.columns[0]].str.cat(others=in_table.loc[:, in_table.columns[1]],
                                                                              sep='|')
                in_table.drop(in_table.columns[0:2], axis=1, inplace=True)
            else:
                in_table.index = in_table[in_table.columns[0]]
                in_table.drop(in_table.columns[0], axis=1, inplace=True)

            in_table.index.name = "combined"

            return(in_table)

        in_tab1 = set_picrust2_tab_index(in_tab1)
        in_tab2 = set_picrust2_tab_index(in_tab2)

        # Concatenate these tables.
        in_tab = pd.concat([in_tab1, in_tab2], sort=True)

    # Determine name of first column of output table.
    if index_name == 'pathway':
        first_col = '# Pathway'
    elif index_name == 'gene' or index_name == 'function':
        first_col = '# Gene Family'
    else:
        sys.exit('Error: first column of input PICRUSt2 tables is \"' +
                 index_name + '\". One of \"pathway\", \"gene\" or '
                 '\"function\" were expected.')

    # Remove description column if it is present.
    if 'description' in in_tab.columns:
        in_tab.drop('description', axis=1, inplace=True)

    # Sort index labels.
    in_tab.sort_index(axis=0, inplace=True)

    # Make output directory.
    make_output_dir(outfolder)

    # Write each sample to a different file in the output folder.
    for sample_id in in_tab.columns:
        outfile = join(outfolder, sample_id + "_humann2-format.tsv")
        
        # Subset to this sample only and remove all rows that are 0.
        tab_subset = in_tab[[sample_id]]
        tab_subset = tab_subset.loc[~(tab_subset == 0).all(axis=1)]

        tab_subset.to_csv(path_or_buf=outfile,  sep="\t",
                          index_label=first_col)


def convert_picrust2_to_humann2_merged(infiles, outfile):
    '''Reads in a PICRUSt2 table(s), combines them, and outputs in HUMAnN2
    tabular format.'''

    # Initialize empty object.
    new_tab = None

    # List for keeping track name of first column in each file.
    infile_index_names = []
        
    for infile in infiles:

        in_table = pd.read_csv(infile, sep="\t")

        infile_index_names.append(in_table.columns[0])

        # If sequence is a column then set the first 2 columns to be the index
        # labels. Otherwise just take the first column to be the index labels.
        if 'sequence' in in_table.columns:
            in_table.index = in_table.loc[:, in_table.columns[0]].str.cat(others=in_table.loc[:, in_table.columns[1]],
                                                                          sep='|')
            in_table.drop(in_table.columns[0:2], axis=1, inplace=True)
        else:
            in_table.index = in_table[in_table.columns[0]]
            in_table.drop(in_table.columns[0], axis=1, inplace=True)

        # Remove description column if it is present.
        if 'description' in in_table.columns:
            in_table.drop('description', axis=1, inplace=True)

        # Concantenate into new table if it exists already.
        if new_tab is not None:
            new_tab = pd.concat([new_tab, in_table], axis=1, sort=True)
        else:
            new_tab = in_table

    # Check for what first column name should be.
    infile_index_names_set = set(infile_index_names)

    if len(infile_index_names_set) > 1:
        sys.exit('Error input PICRUSt2 tables are not all for the same '
                 'datatype. The datatypes of the input files are: ' +
                 ', '.join(list(infile_index_names_set)))

    index_name = infile_index_names[0]

    # Fill in zeros for missing data.
    new_tab = new_tab.fillna(0)

    # Determine name of first column of output table.
    if index_name == "pathway":
        first_col = "# Pathway"
    elif index_name == "gene" or index_name == "function":
        first_col = "# Gene Family"
    else:
        sys.exit('Error: first column of input HUMAnN2 files is \"' +
                 index_name + '\". Either \"# Pathway\" or ' +
                 '\"# Gene Family\" was expected.')

    # Write out.
    new_tab.to_csv(path_or_buf=outfile,  sep="\t", index_label=first_col)


class TemporaryDirectory(object):
    '''Create and return a temporary directory.  This has the same
    behavior as mkdtemp but can be used as a context manager.  For
    example:
        with TemporaryDirectory() as tmpdir:
            ...
    Upon exiting the context, the directory and everything contained
    in it are removed.

    NOTE: This function was taken and modified from the tempfile package to
    first change permissions on folder to be deleted.'''

    def __init__(self, suffix=None, prefix=None, dir=None):
        self.name = tempfile.mkdtemp(suffix, prefix, dir)
        self._finalizer = _weakref.finalize(
            self, self._cleanup, self.name,
            warn_message="Implicitly cleaning up {!r}".format(self))

    @classmethod
    def _cleanup(cls, name, warn_message):
        _shutil.rmtree(name)
        _warnings.warn(warn_message, ResourceWarning)

    def __repr__(self):
        return "<{} {!r}>".format(self.__class__.__name__, self.name)

    def __enter__(self):
        return self.name

    def __exit__(self, exc, value, tb):
        self.cleanup()

    def cleanup(self):
        if self._finalizer.detach():

            # Line added by Gavin Douglas to change permissions to 777 before
            # deleting:
            call(["chmod", "-R", "777", self.name])

            _shutil.rmtree(self.name)
