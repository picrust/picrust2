#!/usr/bin/env python

from __future__ import division

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2018, The PICRUSt Project"
__credits__ = ["Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "2-alpha.5"

from os.path import splitext
from cogent import LoadTree

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


def filter_table_by_presence_in_tree(tree,trait_table_lines,name_field_index = 0,delimiter="\t"):
    """yield lines of a trait table lacking organisms missing from the tree"""
    tree_tips = [node.Name.strip() for node in tree.preorder()]

    for fields in yield_trait_table_fields(trait_table_lines,delimiter):
        curr_name = fields[name_field_index].strip()
        if curr_name not in tree_tips:
            continue
        yield delimiter.join(fields)+"\n"

def convert_trait_values(trait_table_lines,name_field_index=0,delimiter="\t",conversion_fn = int):
    """Convert trait values by running conversion_fn on each"""
    for fields in yield_trait_table_fields(trait_table_lines,delimiter):
        new_fields = []
        for i,field in enumerate(fields):
            if i != name_field_index:
                new_fields.append(str(conversion_fn(float(field))))
            else:
                new_fields.append(field)
        yield delimiter.join(new_fields)+"\n"





def yield_trait_table_fields(trait_table_lines,delimiter="\t",skip_comment_lines=True,max_field_len=100):
    """Yield fields from trait table lines"""
    for line in trait_table_lines:
        #print "Parsing line:\n",line[0:min(100,len(line))],"..."
        if line.startswith("#") and skip_comment_lines:
            continue

        if delimiter not in line:
            delimiters_to_check = {"tab":"\t","space":"","comma":","}
            possible_delimiters = []
            for delim in delimiters_to_check.keys():
                if delimiters_to_check[delim] in line:
                    possible_delimiters.append(delim)
            error_line = "Delimiter '%s' not in line.  The following delimiters were found:  %s.  Is the correct delimiter one of these?"
            raise RuntimeError(error_line % (delimiter,",".join(possible_delimiters)))

        fields = line.split(delimiter)
        yield fields



def ensure_root_is_bifurcating(tree,root_name='root'):
    """Remove child node of root if it is a single child"""
    root_node = tree.getNodeMatchingName(root_name)
    if len(root_node.Children) == 1:
        print "rerooting to avoid monotomy at root"
        tree = tree.rootedAt(root_node.Children[0].Name)
        #tree.remove(root_node)
    tree.prune()

    return tree

def filter_tree_tips_by_presence_in_table(tree,trait_table_lines,name_field_index = 0,delimiter="\t"):
    """yield a tree lacking organisms missing from the trait table"""
    org_ids_in_trait_table = []
    new_tree = tree.deepcopy()

    for fields in yield_trait_table_fields(trait_table_lines, delimiter):
        curr_org = fields[name_field_index].strip()
        org_ids_in_trait_table.append(curr_org)


    # Build up a list of tips to prune
    tips_to_prune = []
    n_tips_not_to_prune = 0
    for tip in tree.iterTips():
        if tip.Name.strip() not in org_ids_in_trait_table:
            tips_to_prune.append(tip.Name)
            #print tip.Name
            #print org_ids_in_trait_table[0]
        else:
            n_tips_not_to_prune += 1
    if not n_tips_not_to_prune:
        raise RuntimeError(\
          "filter_tree_tips_by_presence_in_table:  operation would remove all tips.  Is this due to a formatting error in inputs?")

    # print "Tips to prune:\n\n%s" % tips_to_prune

    #TODO: This should be handled by the exclude_tip function (currently in make_test_trees.py)
    #(it has better error handling)
    for tip_name in tips_to_prune:
        tip = new_tree.getNodeMatchingName(tip_name)
        if tip.Parent is not None:
            removal_ok = tip.Parent.remove(tip)
        else:
            removal_ok = False
        new_tree.prune()

    return new_tree

def print_node_summary_table(input_tree):
    """Print a summary of the name,children,length, and parents of each node"""
    for node in input_tree.postorder():
        if node.Parent:
            parent_name = node.Parent.Name
        else:
            parent_name = None
        print "\t".join(map(str,[node.Name,len(node.Children),node.Length,parent_name]))


def add_to_filename(filename,new_suffix,delimiter="_"):
    """Add to a filename, preserving the extension"""
    filename, ext = splitext(filename)
    new_filename = delimiter.join([filename,new_suffix])
    return "".join([new_filename,ext])

def format_biom_table(biom_table):
    """ Given a biom-format Table object, returns that Table as a BIOM string"""
    generated_by_str = "PI-CRUST " + __version__
    return biom_table.getBiomFormatJsonString(generated_by_str)
