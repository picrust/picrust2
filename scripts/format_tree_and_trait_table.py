#!/usr/bin/env python
# File created on 15 Jul 2011
from __future__ import division

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
__credits__ = ["Jesse Zaneveld","Morgan Langille"]
__license__ = "GPL"
__version__ = "2-alpha.1"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"

from os.path import join,splitext
from cogent.parse.tree import DndParser
from cogent.util.option_parsing import parse_command_line_parameters,\
    make_option
from picrust.format_tree_and_trait_table import *
from picrust.util import make_output_dir, PicrustNode
from picrust.parse import parse_trait_table


# Set up commandline parameters
script_info = {}
script_info['brief_description'] = "Formatting script for filtering and reformatting trees and trait tables."
script_info['script_description'] =\
  """Reformats scripts and trait tables.  Optional fixes include:
        -- Add short (epsilon) branch lengths in place of 0 length branches
        -- Filter out taxa that don't match between tree and trait table
        -- Output tree in NEXUS format
        -- Ensure tree is bifurcating (remove polytomies using very short branches)
        -- Convert floating point trait values to integers
        -- Add a short branch length to the root branch (required by BayesTraits)
        -- Remove internal node names (required by BayesTraits)
        """

script_info['script_usage'] = [\
    ("Example 1","Reformat a tree and trait table with default options:","%prog -i traits.tab -t tree.nwk -o ./format_output/")]
script_info['output_description']= "Outputs a reformatted tree and trait table."
script_info['required_options'] = [\
          make_option('-t','--input_tree',type="existing_filepath",help='the input tree (Newick format)'),\
          make_option('-i','--input_trait_table',type="existing_filepath",help='the input trait table (QIIME OTU table format)')
                  ]

delimiter_choices = ['tab','space','comma']
script_info['optional_options'] = [\
          make_option('-m','--tree_to_trait_mapping',default=None,type="existing_filepath",help='a two-column, tab-delimited text file mapping identifiers in the tree(column 1) to identifiers in the trait table (column 2). If supplied, the identifiers in the trait table will be converted to match the identifiers in the tree. (This mapping does not need to be supplied if the tree and trait table already use a common set of identifiers.) [default: %default]'),\
          make_option('-o','--output_dir',default='./formatted/',type="new_filepath",help='the output directory [default: %default]'),\
          make_option('--input_table_delimiter',default='tab',type="choice",choices=delimiter_choices,\
            help='The character delimiting fields in the input trait table. Valid choices are:'+','.join(delimiter_choices)+' [default: %default]'),\
          make_option('--output_table_delimiter',default='tab',type="choice",choices=delimiter_choices,\
            help='The character delimiting fields in the output trait table. Valid choices are:'+','.join(delimiter_choices)+' [default: %default]'),\
          make_option('--suppress_bifurcating',default=False,action="store_true",help="If set, don't ensure that tree is fully bifurcating. [default: %default]"),\
          make_option('-n','--convert_to_nexus',default=False,action="store_true",help='Convert tree to NEXUS format, including a translate block mapping tip names to numbers. [default: %default]'),\
          make_option('-c','--convert_values_to_ints',default=False,action="store_true",help='Convert the values for each character state to integers. [default: %default]'),\
          make_option('--no_minimum_branch_length',default=False,action="store_true",help="If set, don't ensure all branches have at least a small but non-zero branchlength. [default: %default]"),\
          make_option('--supress_tree_filter',default=False,action="store_true",help="If set, don't filter out tree tips that aren't listed in the trait table [default: %default]"),\
          make_option('--supress_table_filter',default=False,action="store_true",help="If set, don't filter out trait table entries that aren't listed in the tree [default: %default]"),\
          make_option('-r','--add_branch_length_to_root',default=False,action="store_true", help='Add a short branch to the root node (this is required by some phylogeny programs).  The length of the branch is determined by the min_branch_length option  [default: %default]'),\
              make_option('-l','--limit_tree_to_otus_fp',type="existing_filepath",help='Will prune the reference tree to contain only those tips that are within the given OTU table')\

           ]
script_info['version'] = __version__



def main():

    # Parse input to get parameters
    option_parser, opts, args =\
        parse_command_line_parameters(**script_info)

    tree_file = opts.input_tree
    trait_table_fp = opts.input_trait_table
    verbose = opts.verbose

    #Set output base file names
    trait_table_base = 'trait_table.tab'
    pruned_tree_base = 'pruned_tree.newick'
    reference_tree_base = 'reference_tree.newick'

    output_dir = make_output_dir(opts.output_dir,strict=False)
    output_table_fp = join(output_dir,trait_table_base)
    output_tree_fp = join(output_dir,pruned_tree_base)
    output_reference_tree_fp = join(output_dir,reference_tree_base)

    #Handle parameters with more complex defaults
    delimiter_map = {"space":" ","tab":"\t","comma":","}
    input_delimiter = delimiter_map[opts.input_table_delimiter]
    output_delimiter = delimiter_map[opts.output_table_delimiter]

    if verbose:
        print "Running with options:"
        print "\t%s:%s" %("Tree file",tree_file)
        print "\t%s:%s" %("Trait table",trait_table_fp)
        print "\t%s:%s" %("Output tree",output_tree_fp)
        print "\t%s:%s" %("Output reference tree",output_reference_tree_fp)
        print "\t%s:%s" %("Output trait table",output_table_fp)
        print "\t%s:%s" %("Add branch length to root",opts.add_branch_length_to_root)
        print "\t%s:%s" %("Convert to NEXUS?",opts.convert_to_nexus)
        print "\t%s:%s" %("Input trait table delimiter",opts.input_table_delimiter)
        print "\t%s:%s" %("Output trait table delimiter",opts.output_table_delimiter)

    # Begin reformatting

    root_name = "root"

    if opts.no_minimum_branch_length:
        min_branch_length = None
    else:
        min_branch_length = 0.0001

    #Load inputs
    if verbose:
        print "Loading tree...."

    input_tree = DndParser(open(tree_file))

    if verbose:
        print "Loading trait table..."
    trait_table = open(trait_table_fp,"U")
    trait_table_lines = trait_table.readlines()
    if not trait_table_lines:
        raise IOError("No lines could be loaded from file %s. Please check the input file." %trait_table_fp)

    #Get id mappings from mapping file
    if opts.tree_to_trait_mapping:
        if verbose:
            print "Loading tree to trait table mapping file..."

        mapping_file = open(opts.tree_to_trait_mapping,"U")

        trait_to_tree_mapping =\
          make_id_mapping_dict(parse_id_mapping_file(mapping_file))

    else:
        if verbose:
            print "No tree to trait mapping file specified.  Assuming tree tip names and trait table names will match exactly."
        trait_to_tree_mapping = None

    # Call reformatting function using specified parameters
    # to get reference tree
    if opts.verbose:
        print """**BUILDING REFERENCE TREE (without respect to trait table)**"""

    new_reference_tree, not_useful_trait_table_lines =\
      reformat_tree_and_trait_table(\
      tree=input_tree,\
      trait_table_lines = [],\
      trait_to_tree_mapping = None,\
      input_trait_table_delimiter= None,\
      output_trait_table_delimiter= None,\
      filter_table_by_tree_tips=False,\
      convert_trait_floats_to_ints=False,\
      filter_tree_by_table_entries=False,\
      convert_to_bifurcating=True,\
      add_branch_length_to_root=False,\
      name_unnamed_nodes=True,\
      min_branch_length=min_branch_length,\
      verbose=opts.verbose)

    #Make a copy
    new_reference_tree_copy=new_reference_tree.deepcopy()

    if opts.verbose:
        print """**BUILDING PRUNED TREE AND TRAIT TABLE**"""
    # Call reformatting function using specified parameters
    new_tree, new_trait_table_lines = \
       reformat_tree_and_trait_table(tree=new_reference_tree_copy,\
       trait_table_lines = trait_table_lines,\
       trait_to_tree_mapping = trait_to_tree_mapping,\
       input_trait_table_delimiter= input_delimiter,\
       output_trait_table_delimiter=output_delimiter,\
       filter_table_by_tree_tips=True,\
       convert_trait_floats_to_ints=False,\
       filter_tree_by_table_entries=True,\
       convert_to_bifurcating=False,\
       add_branch_length_to_root=False,\
       name_unnamed_nodes=False,\
       min_branch_length=min_branch_length,\
       verbose=opts.verbose)



    #Alter reference tree to only contain tips in OTU table (and of course trait table)
    if opts.limit_tree_to_otus_fp:
        if opts.verbose:
            print "Pruning reference tree to contain only tips in OTU table (and trait table)...."
        otu_table = open(opts.limit_tree_to_otus_fp,"U")
        otu_table_lines = otu_table.readlines()
        header_line,otu_table_fields =parse_trait_table(otu_table_lines,delimiter = input_delimiter,has_header=False)
        header_line,trait_table_fields =\
         parse_trait_table(new_trait_table_lines,delimiter = input_delimiter)


        tips_to_keep = list(otu_table_fields) + list(trait_table_fields)
        tips_to_keep_in_tree = filter_table_by_presence_in_tree(new_reference_tree_copy,tips_to_keep)
        new_reference_tree = filter_tree_tips_by_presence_in_table(new_reference_tree_copy,\
          tips_to_keep_in_tree,verbose=opts.verbose)


    if opts.verbose:
        print "Almost finished. Writing trees and trait table to files..."
    #Write results to files

    # Open output files
    output_trait_table_file = open(output_table_fp,"w+")
    output_tree_file  = open(output_tree_fp,"w+")
    output_reference_tree_file  = open(output_reference_tree_fp,"w+")


    #Output trait table file

    if opts.verbose:
        print "Writing trait table to:", output_table_fp

    output_trait_table_file.write("\n".join(new_trait_table_lines))
    trait_table.close()
    output_trait_table_file.close()

    #Output tree file
    if opts.verbose:
        print "Writing pruned tree to:", output_tree_fp

    if opts.convert_to_nexus is True:
        lines = nexus_lines_from_tree(new_tree)
        output_tree_file.write("\n".join(map(str,lines)))
    else:
        output_tree_file.write(new_tree.getNewick(with_distances=True))

    output_tree_file.close()


    if opts.verbose:
        print "Writing reference tree to:", output_reference_tree_fp
    #Output reference tree file
    output_reference_tree_file.write(new_reference_tree.getNewick(with_distances=True))
    output_reference_tree_file.close()

def load_picrust_tree(tree_fp, verbose):
    """Safely load a tree for picrust"""
    if verbose:
        print "Loading tree..."
    #PicrustNode seems to run into very slow/memory intentsive perfromance...
    #tree = DndParser(open(opts.input_tree),constructor=PicrustNode)
    tree = DndParser(open(tree_fp),constructor=PicrustNode)
    label_conversion_fns = set_label_conversion_fns(verbose=verbose)

    tree = fix_tree_labels(tree,label_conversion_fns)
    return tree

def load_tab_delimited_trait_table(trait_table_fp,verbose=False):
    """Load a tab delimited trait table for picrust"""
    input_trait_table = open(trait_table_fp,"U")
    if verbose:
        print "Parsing trait table..."
    #Find which taxa are to be used in tests
    #(by default trait table taxa)
    trait_table_header,trait_table_fields = \
            parse_trait_table(input_trait_table)

    label_conversion_fns = set_label_conversion_fns(verbose=verbose)
    trait_table_fields = convert_trait_table_entries(trait_table_fields,\
      value_conversion_fns = [],\
      label_conversion_fns = label_conversion_fns)

    trait_table_fields = [t for t in trait_table_fields]

    if verbose:
        print "Number of trait table fields with single quotes:",\
          len([t for t in trait_table_fields if "'" in t[0]])

    return trait_table_header,trait_table_fields

if __name__ == "__main__":
    main()


