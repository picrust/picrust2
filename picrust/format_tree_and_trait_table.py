#!/usr/bin/env python
# File created on 15 Jul 2011
from __future__ import division

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2015, The PICRUSt Project"
__credits__ = ["Jesse Zaneveld", "Morgan Langille"]
__license__ = "GPL"
__version__ = "2-alpha.1"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"

from os.path import splitext
from string import maketrans
from sys import getrecursionlimit,setrecursionlimit
import re
from cogent.parse.tree import DndParser
from cogent.util.option_parsing import parse_command_line_parameters,\
    make_option

from picrust.parse import parse_trait_table,yield_trait_table_fields
from util import PicrustNode

def reformat_tree_and_trait_table(tree,trait_table_lines,trait_to_tree_mapping,\
    input_trait_table_delimiter="\t", output_trait_table_delimiter="\t",\
    filter_table_by_tree_tips=True, convert_trait_floats_to_ints=False,\
    filter_tree_by_table_entries=True,convert_to_bifurcating=False,\
    add_branch_length_to_root=False, name_unnamed_nodes=True,\
    remove_whitespace_from_labels = True,replace_ambiguous_states=True,\
    replace_problematic_label_characters = True,min_branch_length=0.0001,\
    verbose=True):
    """Return a full reformatted tree,pruned reformatted tree  and set of trait table lines

    tree - a PyCogent PhyloNode tree object

    trait_table_lines -- the lines of a trait table, where
      the rows are organisms and the columns are traits (e.g. gene counts).

    trait_id_to_tree_mapping -- a dict keyed by trait table ids, with
      values of tree ids.   If provided, trait table ids will be mapped to
      tree ids

    filter_table_by_tree_tips -- if True, remove trait table rows that don't map to ids on the
    tree

    convert_trait_floats_to_ints -- if True, convert floating point values in trait table cells to integers.

    filter_tree_by_table_entries -- if True, save only the subtree that encompasses organisms in the trait table.
    (equivalent to removing all tips in the tree that don't map to the trait table)

    convert_to_bifurcating -- if True, ensure that the tree is fully bifurcating by resolving polytomies with very short
    branches.

    add_branch_length_to_root -- if True, ensure that the root node has a minimum branch length

    name_unnamed_nodes -- if True, name unnamed nodes in the tree.   (Useful for ensuring internal nodes can be
    consistently identified in both the reference and pruned trees)

    remove_whitespace_from_labels -- if True, replace whitespace in organism labels with underscores

    replace_ambiguous_states -- if True, replace various strings representing ambiguous character states,
    as well as '-1' or -1 (used by IMG to represent a lack of data) with 0 values.

    replace_problematic_table_chars -- if True, replace ':' and ';' in the results with '_', and remove double quotes.
    (AncSR methods like ace can't handle these characters in organism labels)

    min_branch_length -- set the minimum branch length for all edges in the tree.

    This function combines the various reformatting functions in the
    library into a catch-all reformatter.

    TODO: This function is monolithic, so despite the individual
    parts being tested seperately, it probably needs to be broken
    down into several modular parts.  This would need to be done
    with care however, as the order of steps matters quite a bit.


    """



    input_tree = tree

    #Parse lines to fields once


    if trait_table_lines:
        if verbose:
            print "Parsing trait table...."
        header_line,trait_table_fields =\
          parse_trait_table(trait_table_lines,delimiter = input_trait_table_delimiter)
    else:
        if verbose:
            print "Found no trait table lines. Setting data and header to empty"
        trait_table_fields = []
        header_line = ''

    # Tree reformatting
    if convert_to_bifurcating:
        if verbose:
            print "Converting tree to bifurcating...."

        #maximum recursion depth on large trees
        #Try working around this issue with a large
        #recursion depth limit
        old_recursion_limit = getrecursionlimit()
        setrecursionlimit(50000)
        input_tree = input_tree.bifurcating() # Required by most ancSR programs
        setrecursionlimit(old_recursion_limit)

        #input_tree = ensure_root_is_bifurcating(input_tree)

        # The below nutty-looking re-filtering step is necessary
        # When ensuring the root is bifurcating, internal nodes can
        #get moved to the tips so without additional filtering we
        #get unannotated tip nodes

        #if filter_tree_by_table_entries:
        #    input_tree = filter_tree_tips_by_presence_in_table(input_tree,\
        #      trait_table_fields,delimiter=input_trait_table_delimiter)



    #Name unnamed nodes
    if name_unnamed_nodes:
        if verbose:
            print "Naming unnamed nodes in the reference tree...."
        input_tree=make_internal_nodes_unique(input_tree)
        #input_tree.nameUnnamedNodes()
        check_node_labels(input_tree,verbose=verbose)
        #Paranoid check for missing names:
        #if verbose:
        #    print "Checking that all nodes were named..."
        #for i,n in enumerate(input_tree.preorder()):
        #    if n.Name is None:
        #        raise ValueError('Node #%s (in tree.preorder()) was not named!'%str(i))


    #map trait table ids to tree ids
    if trait_to_tree_mapping:
        #if verbose:
        #    print "Validating that trait --> tree mappings match tree ids..."
        #    good,bad = validate_trait_table_to_tree_mappings(input_tree,\
        #      trait_to_tree_mapping.values(), verbose = True)
        #    print "Found %i valid ids." %(len(good))
        #    print "Found %i invalid ids." %(len(bad))
        #    #if bad:
        #    #    raise RuntimeError("The following putative tree ids in mapping file aren't actually in the input tree: %s" % bad)


        if verbose:
            print "Remapping trait table ids to match tree ids...."

        trait_table_fields =\
          remap_trait_table_organisms(trait_table_fields,trait_to_tree_mapping,\
          verbose = verbose)

    label_conversion_fns =\
      set_label_conversion_fns(remove_whitespace_from_labels=remove_whitespace_from_labels,\
        replace_problematic_label_characters=replace_problematic_label_characters)

    value_conversion_fns = set_value_conversion_fns(replace_ambiguous_states=replace_ambiguous_states,\
      convert_trait_floats_to_ints=convert_trait_floats_to_ints)


    #Apply both label and value converters to the trait table
    trait_table_fields = convert_trait_table_entries(\
      trait_table_fields,\
      value_conversion_fns = value_conversion_fns,\
      label_conversion_fns = label_conversion_fns)


    #We now need to apply any formatting functions to the tree nodes as well, to ensure
    #that names are consistent between the two.

    if label_conversion_fns:
        input_tree = fix_tree_labels(input_tree, label_conversion_fns)

    #Then filter the trait table to include only tree tips
    if filter_table_by_tree_tips:
        if verbose:
            print "Filtering trait table ids to include only those that match tree ids...."
        trait_table_fields = filter_table_by_presence_in_tree(input_tree,\
          trait_table_fields,delimiter=input_trait_table_delimiter)

        #if verbose:
        #    print "Verifying that new trait table ids match tree:"
        #    print "# of trait_table_lines: %i" %len(trait_table_lines)
        #    all_tip_ids = [tip.Name for tip in input_tree.iterTips()]
        #    print "example tree tip ids:",all_tip_ids[0:10]
    if filter_tree_by_table_entries:
        if verbose:
            print "filtering tree tips to match entries in trait table...."
        input_tree = filter_tree_tips_by_presence_in_table(input_tree,\
          trait_table_fields,delimiter=input_trait_table_delimiter,\
          verbose=verbose)

    if min_branch_length:
        if verbose:
            print "Setting a min branch length of %f throughout tree...." \
              % min_branch_length
        input_tree = set_min_branch_length(input_tree,min_length = min_branch_length)

    if add_branch_length_to_root:
        if vebose:
            print "Adding a min branch length of %f to the root node...." \
              % min_branch_length
        input_tree = add_branch_length_to_root(input_tree,root_name=input_tree.Name,\
          root_length=min_branch_length)
    if verbose:
        print "Performing a final round of tree pruning to remove internal nodes with only one child...."

    input_tree.prune()




    #Format resulting trait table lines
    result_trait_table_lines = [header_line]
    result_trait_table_lines.extend([output_trait_table_delimiter.join(f) for f in trait_table_fields])

    if verbose:
        print "Final reprocessing of trait table lines to remove trailing whitespace..."
    result_trait_table_lines =\
      [line.strip() for line in result_trait_table_lines if line.strip()]



    if verbose:
        print "Done reformatting tree and trait table"


    return input_tree, result_trait_table_lines

def check_node_labels(input_tree,verbose=False):
    """Check that all nodes are named!"""
    if verbose:
        print "Checking that all nodes were named..."
    for i,n in enumerate(input_tree.preorder()):
        print i,n.Name, n.NameLoaded
        if n.Name is None:
            err_text = 'WARNING: Node #%s (in tree.preorder()) was not named!.  Node properties: %s'%(str(i),str(dir(n)))
            print err_text

def set_label_conversion_fns(remove_whitespace_from_labels=True,\
  replace_problematic_label_characters=True,verbose=False):
    """Return a list of functions for formatting tree node or trait table labels"""
     #Set the functions that will be applied to trait table labels
    label_conversion_fns = []
    if remove_whitespace_from_labels:
        if verbose:
            print "Removing whitespace from trait table organism labels..."
        label_conversion_fns.append(remove_spaces)


    if replace_problematic_label_characters:
        #  Replace ambiguous characters with
        replacement_dict ={":":"_",";":"_"}
        if verbose:
            print "Replacing problematic labels in organism labels:"
            for k,v in replacement_dict.items():
                print k,'-->',v

        chars_to_delete = """'"'"""
        replace_problematic_chars_fn =\
           make_char_translation_fn(replacement_dict,chars_to_delete)
        label_conversion_fns.append(replace_problematic_chars_fn)
    return label_conversion_fns


def set_value_conversion_fns(replace_ambiguous_states=True,\
      convert_trait_floats_to_ints=False,verbose=False):
    """Return a list of value conversion functions for trait table values

     replace_ambiguous_states -- if True, replace values of -,
       -1,'-1','NULL' or None to 0

     convert_trait_floats_to_ints -- if True convert floats to ints

     verbose -- print verbose output describing the conversion fns

    """
    #Set the functions that will be applied to trait table values
    value_conversion_fns = []

    if replace_ambiguous_states:
        #  Replace ambiguous characters with 0's
        replacement_dict ={'-':0,'-1':0,-1:0,'NULL':0,None:0}
        if verbose:
            print "Replacing ambiguous characters:"
            for k,v in replacement_dict.items():
                print k,'-->',v

        replace_ambig_fn = make_translate_conversion_fn(replacement_dict)
        value_conversion_fns.append(replace_ambig_fn)


    if convert_trait_floats_to_ints:
        value_conversion_fns.append(lambda x: str(int(float(x))))

        if verbose:
            print "Converting floating point trait table values to integers...."

    return value_conversion_fns




def fix_tree_labels(tree,label_conversion_fns,verbose=False):
    """Fix tree labels by removing problematic characters"""
    if verbose:
        print "reformatting tree node names..."
    tree = format_tree_node_names(tree,label_conversion_fns)
    #print "Number of tree tips with single quotes:",len([t.Name for t in tree if "'" in t.Name])
    return tree

def make_internal_nodes_unique(tree,base_name='internal_node_%i'):
    """ Removes names that are not unique for internal nodes.
    First occurence of non-unique node is kept and subsequence ones are set to None"""
    #make a list of the names that are already in the tree
    names_in_use = set()
    for i,node in enumerate(tree.preorder(include_self=True)):
        if node.Name is not None:
            if node.Name in names_in_use:
                node.Name=None
            else:
                names_in_use.add(node.Name)

        if node.Name is None:
            while node.Name is None:
                #Find a unique name by adding integers
                proposed_name = base_name % i
                if proposed_name not in names_in_use:
                    node.Name = proposed_name
                    names_in_use.add(proposed_name)
                    break
                else:
                    i += 1
        #Set this so that the PhyloNode *actually* outputs the Name
        node.NameLoaded = True
    return tree


def format_tree_node_names(tree,label_formatting_fns=[]):
    """Return tree with node names formatted using specified fns

    tree -- a PyCogent PhyloNode tree object

    formatting_fns -- a list of formatting functions that are to
    be called on each node name in the tree, and which each return
    a new node name.

    """

    for n in tree.preorder():
        if n.Name is None:
            continue
        new_node_name = n.Name

        for formatting_fn in label_formatting_fns:
            new_node_name = formatting_fn(new_node_name)

        n.Name = new_node_name

    return tree




def nexus_lines_from_tree(tree):
    """Return NEXUS formatted lines from a PyCogent PhyloNode tree"""
    lines = ["#NEXUS"]
    lines.extend(make_nexus_trees_block(tree))
    return lines

def add_branch_length_to_root(tree, root_name ="root",root_length=0.0001):
    """Add branch length to the root of a tree if it's shorter than root_length
    tree -- A PyCogent PhyloNode object
    root_name -- the name of the root node
    root_length -- the desired minimum root length
    This is required by some programs such as BayesTraits"""

    root = tree.getNodeMatchingName(root_name)
    root.Length = max(root.Length,root_length)
    return tree


def set_min_branch_length(tree,min_length= 0.0001):
    """Return tree modified so that all branchlengths are >= min_length.

    tree -- a PyCogent PhyloNode object"""

    for node in tree.preorder():
        if not node.Parent:
            continue
        node.Length = max(node.Length,min_length)
    return tree


def make_nexus_trees_block(tree):
    """Generate a NEXUS format 'trees' block for a given tree

    WARNING:  Removes names from internal nodes, as these cause problems
    downstream
    """

    # First generate the mappings for the NEXUS translate command
    trees_block_template =\
      ["begin trees;",\
      "\ttranslate"]
    name_mappings = {}
    line = None
    for i,node in enumerate(tree.iterTips()):
        name_mappings[node.Name] = i
        if line:
            trees_block_template.append(line)

        line = "\t\t%i %s," %(i,node.Name)
    # The last line needs a semicolon rather than a comma
    line = "\t\t%i %s;" %(i,node.Name)
    trees_block_template.append(line)


    # Reformat tree newick such that names match NEXUS translation table
    for name_to_fix in name_mappings.keys():
        node_to_rename = tree.getNodeMatchingName(name_to_fix)
        node_to_rename.Name=name_mappings[name_to_fix]
    for nonTipNode in tree.iterNontips():
        nonTipNode.Name=''



    tree_newick = tree.getNewick(with_distances=True)
    #for name_to_fix in name_mappings.keys():
    #    tree_newick = tree_newick.replace(name_to_fix+":",str(name_mappings[name_to_fix])+":")
    #for nonTipNode in tree.iterNontips():
    #    tree_newick = tree_newick.replace(nonTipNode.Name+":","")
    #tree_newick = tree_newick.replace(root_name,"")


    tree_template  = "\t\ttree %s = %s" # tree name then newick string
    line = tree_template % ("PyCogent_tree",tree_newick)
    trees_block_template.append(line)

    trees_block_template.append("end;")
    return trees_block_template

def validate_trait_table_to_tree_mappings(tree,trait_table_ids,verbose=True):
    """Report whether tree ids are even in mapping file"""
    good = []
    bad = []
    nodes = [n.Name for n in tree.iterTips()]
    for tt_id in trait_table_ids:
        if tt_id in nodes:
            good.append(tt_id)
        else:
            bad.append(tt_id)
    if verbose:
        print "Of %i ids, %i were OK (mapped to tree)" %(len(trait_table_ids),len(good))
        print "Example good ids",good[0:min(len(good),10)]
        print "Example bad ids",bad[0:min(len(bad),10)]
        print "Example tip ids",nodes[0:min(len(nodes),10)]
    return good,bad

def filter_table_by_presence_in_tree(tree,trait_table_fields,name_field_index = 0,delimiter="\t"):
    """yield lines of a trait table lacking organisms missing from the tree"""

    tree_tips = [str(node.Name.strip()) for node in tree.preorder()]
    #print tree_tips
    result_fields = []
    for fields in trait_table_fields:
        curr_name = fields[name_field_index].strip()
        if curr_name not in tree_tips:
            #print curr_name,"could not be found in tree nodes"
            #print curr_name in tree_tips
            #try:
            #    print int(curr_name) in tree_tips
            #except:
            #    pass
            #print curr_name.strip() in tree_tips
            continue
        result_fields.append(fields)
    return result_fields



def make_translate_conversion_fn(translation_dict):
    """Return a new function that replaces values in input values with output_value
    translation_dict -- a dict that maps inputs that should be translated to
    their appropriate output
    """

    def translate_conversion_fn(trait_value_field):
        # Return translation, or the original value if no translation
        # is available
        try:
            trait_value_field = trait_value_field.strip()
        except AttributeError:
            trait_value_field = str(trait_value_field).strip()

        result = translation_dict.get(trait_value_field,trait_value_field)

        #print trait_value_field
        #print translation_dict.keys()

        if result in translation_dict.keys():
            raise RuntimeError("failed to translate value: %s" % result)

        return str(result)

    return translate_conversion_fn

def make_char_translation_fn(translation_dict,deletion_chars=''):
    """Return a new function that replaces values in input values with output_value
    translation_dict -- a dict that maps inputs that should be translated to
    their appropriate output
    """

    def translate_conversion_fn(trait_value_field):
        # Return translation, or the original value if no translation
        # is available
        trait_value_field = str(trait_value_field).strip()

        from_chars = ''
        to_chars = ''
        for k,v in translation_dict.items():
            from_chars += k
            to_chars += v

        translation_table = maketrans(from_chars,to_chars)
        #print trait_value_field
        #print translation_dict.keys()
        result = trait_value_field.translate(translation_table,deletion_chars)


        if result in translation_dict.keys():
            raise RuntimeError("failed to translate value: %s" % result)

        return str(result)

    return translate_conversion_fn


def remove_spaces(trait_label_field):
    """A conversion function that replaces spaces with underscores in a label
    """

    label = str(trait_label_field)
    fields = trait_label_field.lstrip().strip().split()
    return "_".join(fields)






def convert_trait_table_entries(trait_table_fields,\
        label_conversion_fns=[str],value_conversion_fns = [float]):
    """Convert trait values by running conversion_fns on labels and values

    trait_table_fields -- list of strings (from a trait table line)
      the first field is assumed to be an organism name, and so isn't
      formatted.

    label_conversion_fns -- a list of functions to be run on each
    organism name label (in the order they should be run).  Each
    function should need only a single entry as input, and output
    the resulting label

    value_conversion_fns -- another list of functions, but for
    trait values.  Again these will be run in order on each table
    value.


    """
    name_field_index = 0
    #print "Value conversion fns:",[f.__name__ for f in value_conversion_fns]
    #print "label_conversion_fns:",[f.__name__ for f in label_conversion_fns]
    for fields in trait_table_fields:
        new_fields = []
        for i,field in enumerate(fields):
            if i != name_field_index:
                converters_to_use = value_conversion_fns
            else:
                converters_to_use = label_conversion_fns

            #Run appropriate converters on this field
            new_val = field
            for curr_conv_fn in converters_to_use:
                new_val = str(curr_conv_fn(new_val))
            new_fields.append(new_val)

        yield new_fields



def ensure_root_is_bifurcating(tree,root_name='root',verbose=False):
    """Remove child node of root if it is a single child"""
    root_node = tree.getNodeMatchingName(root_name)
    if len(root_node.Children) == 1:
        if verbose:
            print "Rerooting to avoid monotomy at root"
        tree = tree.rootedAt(root_node.Children[0].Name)
        #tree.remove(root_node)
    tree.prune()

    return tree

def filter_tree_tips_by_presence_in_table(tree,trait_table_fields,name_field_index = 0,\
      delimiter="\t",verbose=True):
    """yield a tree lacking organisms missing from the trait table

    trait_table_fields -- a list of lists, containing the results of parsing the data
    lines of the trait table.  Each set of fields in the list should contain the organism name
    at index 0, and data values for the various traits at other positions

    """
    org_ids_in_trait_table = []
    new_tree = tree.deepcopy()

    for fields in trait_table_fields:
        curr_org = fields[name_field_index].strip()
        org_ids_in_trait_table.append(curr_org)


    # Build up a list of tips to prune
    tips_to_prune = []
    tips_not_to_prune = []
    n_tips_not_to_prune = 0
    for tip in tree.iterTips():
        if tip.Name.strip() not in org_ids_in_trait_table:
            tips_to_prune.append(tip.Name)
        else:
            n_tips_not_to_prune += 1
            tips_not_to_prune.append(tip.Name)

    if verbose and tips_to_prune:
        print "Found %i tips to prune." %(len(tips_to_prune))
        print "Example pruned tree tip names:",tips_to_prune[0:min(len(tips_to_prune),10)]
        print "Example valid org ids:",org_ids_in_trait_table[0:min(len(org_ids_in_trait_table),10)]
    if not n_tips_not_to_prune:
        raise RuntimeError(\
          "filter_tree_tips_by_presence_in_table:  operation would remove all tips.  Is this due to a formatting error in inputs?")
    if verbose:
        print "%i of %i tips will be removed (leaving %i)" %(len(tips_to_prune),\
          n_tips_not_to_prune + len(tips_to_prune), n_tips_not_to_prune)
        print "Example tips that will be removed (first 10):\n\n%s" % \
          tips_to_prune[0:min(len(tips_to_prune),10)]
    new_tree = get_sub_tree(tree,tips_not_to_prune)
    return new_tree


def get_sub_tree(tree,tips_not_to_prune):
    """Get sub tree, modifying recursion limit if necessary"""

    try:
        new_tree = tree.getSubTree(tips_not_to_prune)
    except RuntimeError:
        #NOTE:  getSubTree will hit
        #maximum recursion depth on large trees
        #Try working around this issue with a large
        #recursion depth limit
        old_recursion_limit = getrecursionlimit()
        setrecursionlimit(50000)
        new_tree = tree.getSubTree(tips_not_to_prune)
        setrecursionlimit(old_recursion_limit)
    return new_tree


def print_node_summary_table(input_tree):
    """Print a summary of the name,children,length, and parents of each node"""
    for node in input_tree.postorder():
        if node.Parent:
            parent_name = node.Parent.Name
        else:
            parent_name = None
        yield "\t".join(map(str,[node.Name,len(node.Children),node.Length,parent_name]))


def add_to_filename(filename,new_suffix,delimiter="_"):
    """Add to a filename, preserving the extension"""
    filename, ext = splitext(filename)
    new_filename = delimiter.join([filename,new_suffix])
    return "".join([new_filename,ext])


def make_id_mapping_dict(tree_to_trait_mappings):
    """Generates trait_to_tree mapping dictionary from a list of mapping tuples

    mappings -- in the format tree_id, trait_id

    """
    trait_to_tree_mapping_dict = {}

    for tree_id,trait_id in tree_to_trait_mappings:
        trait_to_tree_mapping_dict[trait_id] = tree_id

    return trait_to_tree_mapping_dict

def parse_id_mapping_file(file_lines,delimiter="\t"):
    """Parse two-column id mapping file, returning a generator of fields"""
    for line in file_lines:
        yield line.strip().split(delimiter)


def remap_trait_table_organisms(trait_table_fields,trait_to_tree_mapping_dict,verbose=False):
    """Yield trait table fields with organism ids substituted using the mapping dict

    An iterator containing lists for each trait.  The first field in each list
    should be the organism id, and the rest should be trait values.


    """

    remapped_fields = []
    bad_ids = []
    default_total = 0
    #if verbose:
    #    print trait_to_tree_mapping_dict
    #    print sorted(list(set(trait_to_tree_mapping_dict.keys())))
    for fields in trait_table_fields:

        try:
            fields[0] = trait_to_tree_mapping_dict[fields[0]]
        except KeyError:
            bad_ids.append(fields[0])
            continue

        remapped_fields.append(fields)

    if verbose and bad_ids:
        print "%i of %i trait table ids could not be mapped to tree" %(len(bad_ids),len(remapped_fields))
        print "Example trait table ids that could not be mapped to tree:" %(bad_ids[:min(len(bad_ids),10)])

    return remapped_fields

def load_picrust_tree(tree_fp, verbose=False):
    """Safely load a tree for picrust"""
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
