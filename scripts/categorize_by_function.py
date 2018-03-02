#!/usr/bin/env python
# File created on 09 Jan 2013
from __future__ import division

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
__credits__ = ["Daniel McDonald", "Morgan Langille", "Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "2-alpha.1"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

import json
import h5py

from cogent.util.option_parsing import parse_command_line_parameters, make_option
from biom import load_table
from biom.table import vlen_list_of_str_formatter
from picrust.util import write_biom_table

script_info = {}
script_info['brief_description'] = "Collapse table data to a specified level in a hierarchy."
script_info['script_description'] = "This script collapses hierarchical data to a specified level. For instance, often it is useful to examine KEGG results from a higher level within the pathway hierarchy. Many genes are sometimes involved in multiple pathways, and in these circumstances (also know as a one-to-many relationship), the gene is counted for each pathway. This has a side effect of increasing the total count of genes in the table."
script_info['script_usage'] = [\
("","Collapse predicted metagenome using KEGG Pathway metadata.","""%prog -i predicted_metagenomes.biom -c KEGG_Pathways -l 3 -o predicted_metagenomes.L3.biom"""),\
("","Change output to tab-delimited format (instead of BIOM).","""%prog -f -i predicted_metagenomes.biom -c KEGG_Pathways -l 3 -o predicted_metagenomes.L3.txt"""),\
("","Collapse COG Categories.","""%prog -i cog_predicted_metagenomes.biom -c COG_Category -l 2 -o cog_predicted_metagenomes.L2.biom"""),\
("","Collapse predicted metagenome using taxonomy metadata (not one-to-many).","""%prog -i observation_table.biom -c taxonomy -l 1 -o observation_table.L1.biom"""),\


]
script_info['output_description']= "Output table is contains gene counts at a higher level within a hierarchy."
script_info['required_options'] = [\
 make_option('-i','--input_fp',type="existing_filepath",help='the predicted metagenome table'),\
 make_option('-o','--output_fp',type='new_filepath', help='the resulting table'),
 make_option('-c','--metadata_category',type='string',help='the metadata category that describes the hierarchy (e.g. KEGG_Pathways, COG_Category, etc.). Note: RFAM predictions can not be collapsed because there are no categories to group them into.'),
 make_option('-l','--level',type='int',help='the level in the hierarchy to collapse to. A value of 0 is not allowed, a value of 1 is the highest level, and any higher value nears the leaves of the hierarchy. For instance, if the hierarchy contains 4 levels, specifying 3 would collapse at one level above being fully specified.')
]
script_info['optional_options'] = [
 make_option('--ignore',type='string',default=None, help="Ignore the comma separated list of names. For instance, specifying --ignore_unknown=unknown,unclassified will ignore those labels while collapsing. The default is to not ignore anything. [default: %default]"),
 make_option('-f','--format_tab_delimited',action="store_true",default=False,help='output the predicted metagenome table in tab-delimited format [default: %default]')]
script_info['version'] = __version__

def make_collapse_f(category, level, ignore):
    """produce a collapsing function for one-to-many relationships"""
    # adjust level such that, for instance, level 1 corresponds to index 0
    if level > 0:
        level -= 1

    def collapse(id_, md):
        if ignore is not None:
            ignore_labels = set(ignore.split(','))
        else:
            ignore_labels = None

        is_single_level = False

        for path in md[category]:

            # need to convert strings to lists (if needed) before checking
            #if they are in the ignore list
            if isinstance(path,basestring):
                # If we have a list of strings, we want the whole thing (only)
                path = md[category]
                is_single_level = True

            if ignore is not None and path[level].lower() in ignore_labels:
                continue

            yield (path[:(level+1)],path[level])

            #If we only have one list of strings, we're done - bail.
            if is_single_level:
                break
    return collapse


def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)

    if opts.level <= 0:
        option_parser.error("level must be greater than zero!")

    collapse_f = make_collapse_f(opts.metadata_category, opts.level,
                                 opts.ignore)
    table = load_table(opts.input_fp)

    if h5py.is_hdf5(opts.input_fp):
        # metadata are not deserializing correctly. Duct tape it.
        update_d = {}
        for i, md in zip(table.ids(axis='observation'),
                         table.metadata(axis='observation')):
            update_d[i] = {k: json.loads(v[0]) for k, v in md.items()}
        table.add_metadata(update_d, axis='observation')

    result = table.collapse(collapse_f, axis='observation', one_to_many=True,
                            norm=False,
                            one_to_many_md_key=opts.metadata_category)

    if(opts.format_tab_delimited):
        f = open(opts.output_fp, 'w')
        f.write(result.to_tsv(header_key=opts.metadata_category,
                              header_value=opts.metadata_category,
                              metadata_formatter=lambda s: '; '.join(s)))
        f.close()
    else:
        format_fs = {opts.metadata_category: vlen_list_of_str_formatter}
        write_biom_table(result, opts.output_fp, format_fs=format_fs)

if __name__ == "__main__":
    main()
