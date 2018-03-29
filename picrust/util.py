#!/usr/bin/env python

from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2015, The PICRUSt Project"
__credits__ = ["Greg Caporaso", "Morgan Langille", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "2-alpha.2"

from biom import parse_table, Table
from biom.table import vlen_list_of_str_formatter
from biom.util import biom_open, HAVE_H5PY
from cogent.core.tree import PhyloNode, TreeError
from contextlib import contextmanager
from json import dumps
from numpy import array, asarray, atleast_1d
from os import fsync, makedirs, remove, rename
from os.path import abspath, dirname, isdir, join, split
import StringIO
from subprocess import PIPE, Popen, call


def make_sample_transformer(scaling_factors):
    def transform_sample(sample_value,sample_id,sample_metadata):
        scaling_factor = scaling_factors[sample_id]
        new_val = sample_value * scaling_factor
        return new_val
    return transform_sample


def scale_metagenomes(metagenome_table,scaling_factors):
    """ scale metagenomes from metagenome table and scaling factors
    """
    transform_sample_f = make_sample_transformer(scaling_factors)
    new_metagenome_table = metagenome_table.transform(transform_sample_f)
    return new_metagenome_table


def convert_precalc_to_biom(precalc_in, ids_to_load=None,transpose=True,md_prefix='metadata_'):
    """Loads PICRUSTs tab-delimited version of the precalc file and outputs a BIOM object"""

    #if given a string convert to a filehandle
    if type(precalc_in) ==str or type(precalc_in) == unicode:
        fh = StringIO.StringIO(precalc_in)
    else:
        fh=precalc_in

    #first line has to be header
    header_ids=fh.readline().strip().split('\t')

    col_meta_locs={}
    for idx,col_id in enumerate(header_ids):
        if col_id.startswith(md_prefix):
            col_meta_locs[col_id[len(md_prefix):]]=idx

    end_of_data=len(header_ids)-len(col_meta_locs)
    trait_ids = header_ids[1:end_of_data]

    col_meta=[]
    row_meta=[{} for i in trait_ids]

    if ids_to_load is not None and len(ids_to_load) > 0:
        ids_to_load=set(ids_to_load)
        load_all_ids=False
    else:
        load_all_ids=True

    matching=[]
    otu_ids=[]
    for line in fh:
        fields = line.strip().split('\t')
        row_id=fields[0]
        if(row_id.startswith(md_prefix)):
            #handle metadata

            #determine type of metadata (this may not be perfect)
            metadata_type=determine_metadata_type(line)
            for idx,trait_name in enumerate(trait_ids):
                row_meta[idx][row_id[len(md_prefix):]]=parse_metadata_field(fields[idx+1],metadata_type)

        elif load_all_ids or (row_id in set(ids_to_load)):
            otu_ids.append(row_id)
            matching.append(map(float,fields[1:end_of_data]))

            #add metadata
            col_meta_dict={}
            for meta_name in col_meta_locs:
                col_meta_dict[meta_name]=fields[col_meta_locs[meta_name]]
            col_meta.append(col_meta_dict)

            if not load_all_ids:
                ids_to_load.remove(row_id)

    if not otu_ids:
        raise ValueError,"No OTUs match identifiers in precalculated file. PICRUSt requires an OTU table reference/closed picked against GreenGenes.\nExample of the first 5 OTU ids from your table: {0}".format(', '.join(list(ids_to_load)[:5]))

    if ids_to_load:
       raise ValueError,"One or more OTU ids were not found in the precalculated file!\nAre you using the correct --gg_version?\nExample of (the {0}) unknown OTU ids: {1}".format(len(ids_to_load),', '.join(list(ids_to_load)[:5]))

    #note that we transpose the data before making biom obj
    matching = asarray(matching)
    if transpose:
        return Table(matching.T, trait_ids, otu_ids, row_meta, col_meta,
                     type='Gene table')
    else:
        return Table(matching, otu_ids, trait_ids, col_meta, row_meta,
                     type='Gene table')


def convert_biom_to_precalc(biom_table):
    """Converts a biom table into a PICRUSt precalculated tab-delimited file """
    col_ids = biom_table.ids(axis='observation')
    row_ids = biom_table.ids()

    lines = []
    header = ['#OTU_IDs'] + list(col_ids)

    col_metadata_names = []
    # peak at metadata for Samples (e.g. NSTI) so we can set the header
    if biom_table.metadata():
        col_metadata_names = biom_table.metadata()[0].keys()

    #add the metadata names to the header
    for col_metadata_name in col_metadata_names:
        header.append('metadata_' + col_metadata_name)

    lines.append(map(str, header))

    row_metadata_names = []
    # peak at metadata for observations (e.g. KEGG_Pathways)
    if biom_table.metadata(axis='observation'):
        row_metadata_names = biom_table.metadata(axis='observation')[0].keys()

    for metadata_name in row_metadata_names:
        metadata_line = ['metadata_' + metadata_name]

    # do the observation metadata now
        for col_id in col_ids:
            metadata = biom_table.metadata(axis='observation')[biom_table.index(col_id, axis='observation')]
            metadata_line.append(biom_meta_to_string(metadata[metadata_name]))
        lines.append(map(str, metadata_line))

    # transpose the actual count data
    transposed_table = biom_table._data.T
    for idx, count in enumerate(transposed_table.toarray()):
        line = [row_ids[idx]] + map(str, count)

        # add the metadata values to the end of the row now
        for meta_name in col_metadata_names:
            line.append(biom_table.metadata()[idx][meta_name])
        lines.append(line)

    return "\n".join("\t".join(map(str, x)) for x in lines)


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


def write_fasta(seq, outfile):
    out_fasta = open(outfile, "w")

    for s in seq:
        out_fasta.write(">" + s + "\n")
        out_fasta.write(seq[s] + "\n")

    out_fasta.close()


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


def determine_metadata_type(line):
    if ';' in line:
        if '|' in line:
            return 'list_of_lists'
        else:
            return 'list'
    else:
        return 'string'


def parse_metadata_field(metadata_str,metadata_format='string'):
    if metadata_format == 'string':
        return metadata_str
    elif metadata_format == 'list':
        return [e.strip() for e in metadata_str.split(';')]
    elif metadata_format == 'list_of_lists':
        return [[e.strip() for e in y.split(';')] for y in metadata_str.split('|')]


def biom_meta_to_string(metadata):
    """ Determine which format the metadata is and then convert to a string"""

    #Note that since ';' and '|' are used as seperators we must replace them if they exist
    if type(metadata) ==str or type(metadata)==unicode:
        return metadata.replace(';',':')
    elif type(metadata) == list:
        if type(metadata[0]) == list:
            return "|".join(";".join([y.replace(';',':').replace('|',':') for y in x]) for x in metadata)
        else:
            return ";".join(x.replace(';',':') for x in metadata)


def system_call(cmd, shell=True):
    """Call cmd and return (stdout, stderr, return_value).

    cmd can be either a string containing the command to be run, or a sequence
    of strings that are the tokens of the command.

    Please see Python's subprocess.Popen for a description of the shell
    parameter and how cmd is interpreted differently based on its value.

    This code was copied from QIIME's qiime_system_call() (util.py) function on June 3rd, 2013.
    """
    proc = Popen(cmd, shell=shell, universal_newlines=True, stdout=PIPE,
                 stderr=PIPE)
    # communicate pulls all stdout/stderr from the PIPEs to
    # avoid blocking -- don't remove this line!
    stdout, stderr = proc.communicate()
    return_value = proc.returncode
    return stdout, stderr, return_value


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


def file_contains_nulls(file):
    """Checks given file for null characters. These are sometimes created on SGE clusters when system IO is overloaded."""

    return '\x00' in open(file,'rb').read()


def parse_table_to_biom(table_lines, table_format="tab-delimited",\
    biom_format = 'otu table'):
    """Read the lines of an open trait table file, and output a .biom table object

    The trait table must be either a biom file, or a picrust tab-delimited file
    table_format -- must be either 'tab-delimited' or 'biom'

    """
    return parse_table(table_lines)


def get_picrust_project_dir():
    """ Returns the top-level PICRUST directory
    """
    # Get the full path of util.py
    current_file_path = abspath(__file__)
    # Get the directory containing util.py
    current_dir_path = dirname(current_file_path)
    # Return the directory containing the directory containing util.py
    return dirname(current_dir_path)


def transpose_trait_table_fields(data_fields,header,id_row_idx=0,\
    input_header_delimiter="\t",output_delimiter="\t"):
    """Transpose the fields of a trait table, returning new data_fields,header

    data_fields:  list of lists for data fields
    header:  a string describing the header_line
    id_row_idx:  index of row labels.  Almost always 0 but included for
    but included for completeness

    input_header_delimiter: delimiter for fields in the header string
    output_delimiter: use this delimiter to join header fields

    NOTE:  typically the header and data fields are generated
    by parse_trait_table in picrust.parse
    """
    header_fields = header.split(input_header_delimiter)

    # ensure no trailing newlines
    old_header_fields = [h.strip() for h in header_fields]
    new_header_fields = [old_header_fields[0]] + \
                        [df[id_row_idx].strip() for df in data_fields]

    non_label_data_fields = []
    for row in data_fields:
        non_label_fields = [e for i, e in enumerate(row) if i != id_row_idx]
        non_label_data_fields.append(non_label_fields)

    data_array = array(non_label_data_fields)
    new_data_array = data_array.T

    new_rows = []
    for i,row in enumerate(new_data_array):
        label = old_header_fields[i+1]
        # this is i+1 not i because i is the blank/meaningless
        # upper left corner entry.
        new_row = [label] + list(row)
        new_rows.append(new_row)
    new_header = output_delimiter.join(new_header_fields)

    return new_header + "\n", new_rows


def make_output_dir_for_file(filepath):
    """Create sub-directories for a new file if they don't already exist"""
    dirpath = dirname(filepath)
    if not isdir(dirpath) and not dirpath == '':
        makedirs(dirpath)


def write_biom_table(biom_table, biom_table_fp, compress=True,
                     write_hdf5=HAVE_H5PY, format_fs=None):
    """Writes a BIOM table to the specified filepath

    Parameters
    ----------
    biom_table : biom.Table
        The table object to write out
    biom_table_fp : str
        The path to the output file
    compress : bool, optional
        Defaults to ``True``. If True, built-in compression on the output HDF5
        file will be enabled. This option is only relevant if ``write_hdf5`` is
        ``True``.
    write_hdf5 : bool, optional
        Defaults to ``True`` if H5PY is installed and to ``False`` if H5PY is
        not installed. If ``True`` the output biom table will be written as an
        HDF5 binary file, otherwise it will be a JSON string.
    format_fs : dict, optional
        Formatting functions to be passed to `Table.to_hdf5`

    Notes
    -----
    This code was adapted from QIIME 1.9
    """
    generated_by = "PICRUSt " + __version__

    if write_hdf5:
        with biom_open(biom_table_fp, 'w') as biom_file:
            biom_table.to_hdf5(biom_file, generated_by, compress,
                               format_fs=format_fs)
    else:
        with open(biom_table_fp, 'w') as biom_file:
            biom_table.to_json(generated_by, biom_file)


def make_output_dir(dirpath, strict=False):
    """Make an output directory if it doesn't exist

    Returns the path to the directory
    dirpath -- a string describing the path to the directory
    strict -- if True, raise an exception if dir already
    exists
    """
    dirpath = abspath(dirpath)

    #Check if directory already exists
    if isdir(dirpath):
        if strict == True:
            err_str = "Directory '%s' already exists" % dirpath
            raise IOError(err_str)

        return dirpath
    try:
        makedirs(dirpath)
    except IOError,e:
        err_str = "Could not create directory '%s'. Are permissions set correctly? Got error: '%s'" %e
        raise IOError(err_str)

    return dirpath


class PicrustNode(PhyloNode):
    def multifurcating(self, num, eps=None, constructor=None):
        """Return a new tree with every node having num or few children

        num : the number of children a node can have max
        eps : default branch length to set if self or constructor is of
            PhyloNode type
        constructor : a TreeNode or subclass constructor. If None, uses self
        """
        if num < 2:
            raise TreeError, "Minimum number of children must be >= 2"

        if eps is None:
            eps = 0.0

        if constructor is None:
            constructor = self.__class__

        if hasattr(constructor, 'Length'):
            set_branchlength = True
        else:
            set_branchlength = False

        new_tree = self.copy()

        for n in new_tree.preorder(include_self=True):
            while len(n.Children) > num:
                new_node = constructor(Children=n.Children[-num:])

                if set_branchlength:
                    new_node.Length = eps

                n.append(new_node)

        return new_tree

    def bifurcating(self, eps=None, constructor=None):
        """Wrap multifurcating with a num of 2"""
        return self.multifurcating(2, eps, constructor)

    def nameUnnamedNodes(self):
        """sets the Data property of unnamed nodes to an arbitrary value

        Internal nodes are often unnamed and so this function assigns a
        value for referencing.
        Note*: This method is faster then pycogent nameUnnamedNodes()
        because it uses a dict instead of an array. Also, we traverse
        only over internal nodes (and not including tips)
        """

        #make a list of the names that are already in the tree
        names_in_use = {}
        for node in self.iterNontips(include_self=True):
            if node.Name:
                names_in_use[node.Name]=1

        #assign unique names to the Data property of nodes where Data = None
        name_index = 1
        for node in self.iterNontips(include_self=True):
            #if (not node.Name) or re.match('edge',node.Name):
            if not node.Name:
                new_name = 'node' + str(name_index)
                #choose a new name if name is already in tree
                while new_name in names_in_use:
                    name_index += 1
                    new_name = 'node' + str(name_index)
                node.Name = new_name
                names_in_use[node.Name]=1
                name_index += 1

    def getSubTree(self,names):
        """return a new subtree with just the tips in names

        assumes names is a set
        assumes all names in names are present as tips in tree
        """
        tcopy = self.deepcopy()

        while len(tcopy.tips()) != len(names):
            # for each tip, remove it if we do not want to keep it
            for n in tcopy.tips():
                if n.Name not in names:
                    n.Parent.removeNode(n)

            # reduce single-child nodes
            tcopy.prune()

        return tcopy

def list_of_list_of_str_formatter(grp, header, md, compression):
    """Serialize [[str]] into a BIOM hdf5 compatible form

    Parameters
    ----------
    grp : h5py.Group
        This is ignored. Provided for passthrough
    header : str
        The key in each dict to pull out
    md : list of dict
        The axis metadata
    compression : bool
        Whether to enable dataset compression. This is ignored, provided for
        passthrough

    Returns
    -------
    grp : h5py.Group
        The h5py.Group
    header : str
        The key in each dict to pull out
    md : list of dict
        The modified metadata that can be formatted in hdf5
    compression : bool
        Whether to enable dataset compression.

    Notes
    -----
    This method is intended to be a "passthrough" to BIOM's
    vlen_list_of_str_formatter method. It is a transform method.
    """
    new_md = [{header: atleast_1d(asarray(dumps(m[header])))} for m in md]
    return (grp, header, new_md, compression)


def picrust_formatter(*args):
    """Transform, and format"""
    return vlen_list_of_str_formatter(*list_of_list_of_str_formatter(*args))


@contextmanager
def atomic_write(file):
    """
    Yields an open temporary file and renames to ``file`` on exit

    This context manager aims to make it more convenient to write to a file
    only when the write operation has completed (e.g. downloading a file from
    the internet).

    The yielded file will be open in 'wb' mode.
    """
    dir_name, basename = split(file)
    tmp_path = join(dir_name, '~{}'.format(basename))

    try:
        with open(tmp_path, 'wb') as f:
            yield f
            f.flush()
            fsync(f.fileno())
    except:
        remove(tmp_path)
        raise
    else:
        rename(tmp_path, file)
