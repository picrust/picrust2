#!/usr/bin/env python
# File created on 22 Feb 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
__credits__ = ["Greg Caporaso","Morgan Langille"]
__license__ = "GPL"
__version__ = "2-alpha.1"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"


from cogent.util.option_parsing import parse_command_line_parameters, make_option
from biom import load_table, Table
from picrust.predict_metagenomes import transfer_observation_metadata
from os import path
from os.path import join
from picrust.util import get_picrust_project_dir, convert_precalc_to_biom,make_output_dir_for_file, write_biom_table
import gzip

script_info = {}
script_info['brief_description'] = "Normalize an OTU table by marker gene copy number"
script_info['script_description'] = ""
script_info['script_usage'] = [
("","Normalize the OTU abundances for a given OTU table picked against the newest version of Greengenes:","%prog -i closed_picked_otus.biom -o normalized_otus.biom"),
("","Change the version of Greengenes used for OTU picking:","%prog -g 18may2012 -i closed_picked_otus.biom -o normalized_otus.biom")
]
script_info['output_description']= "A normalized OTU table"
script_info['required_options'] = [
 make_option('-i','--input_otu_fp',type="existing_filepath",help='the input otu table filepath in biom format'),
 make_option('-o','--output_otu_fp',type="new_filepath",help='the output otu table filepath in biom format'),
]
gg_version_choices=['13_5','18may2012']
script_info['optional_options'] = [
    make_option('-g','--gg_version',default=gg_version_choices[0],type="choice",\
                    choices=gg_version_choices,\
                    help='Version of GreenGenes that was used for OTU picking. Valid choices are: '+\
                    ', '.join(gg_version_choices)+\
                    ' [default: %default]'),

    make_option('-c','--input_count_fp',default=None,type="existing_filepath",\
                    help='Precalculated input marker gene copy number predictions on per otu basis in biom format (can be gzipped).Note: using this option overrides --gg_version. [default: %default]'),
    make_option('--metadata_identifer',
             default='CopyNumber',
             help='identifier for copy number entry as observation metadata [default: %default]'),

    make_option('--load_precalc_file_in_biom',default=False,action="store_true",\
                    help='Instead of loading the precalculated file in tab-delimited format (with otu ids as row ids and traits as columns) load the data in biom format (with otu as SampleIds and traits as ObservationIds) [default: %default]'),

]
script_info['version'] = __version__


def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)

    otu_table = load_table(opts.input_otu_fp)

    ids_to_load = otu_table.ids(axis='observation')

    if(opts.input_count_fp is None):
        #precalc file has specific name (e.g. 16S_13_5_precalculated.tab.gz)
        precalc_file_name='_'.join(['16S',opts.gg_version,'precalculated.tab.gz'])
        input_count_table=join(get_picrust_project_dir(),'picrust','data',precalc_file_name)
    else:
        input_count_table=opts.input_count_fp

    if opts.verbose:
        print "Loading trait table: ", input_count_table

    ext=path.splitext(input_count_table)[1]

    if (ext == '.gz'):
        count_table_fh = gzip.open(input_count_table,'rb')
    else:
        count_table_fh = open(input_count_table,'U')

    if opts.load_precalc_file_in_biom:
        count_table = load_table(count_table_fh)
    else:
        count_table = convert_precalc_to_biom(count_table_fh, ids_to_load)

    #Need to only keep data relevant to our otu list
    ids=[]
    for x in otu_table.iter(axis='observation'):
        ids.append(str(x[1]))

    ob_id=count_table.ids(axis='observation')[0]

    filtered_otus=[]
    filtered_values=[]
    for x in ids:
        if count_table.exists(x, axis='sample'):
            filtered_otus.append(x)
            filtered_values.append(otu_table.data(x, axis='observation'))

    filtered_otu_table = Table(filtered_values, filtered_otus, otu_table.ids())

    copy_numbers_filtered={}
    for x in filtered_otus:
        value = count_table.get_value_by_ids(ob_id,x)
        try:
            #data can be floats so round them and make them integers
            value = int(round(float(value)))

        except ValueError:
            raise ValueError,\
                  "Invalid type passed as copy number for OTU ID %s. Must be int-able." % (value)
        if value < 1:
            raise ValueError, "Copy numbers must be greater than or equal to 1."

        copy_numbers_filtered[x]={opts.metadata_identifer:value}

    filtered_otu_table.add_metadata(copy_numbers_filtered, axis='observation')

    def metadata_norm(v, i, md):
        return v / float(md[opts.metadata_identifer])
    normalized_table = filtered_otu_table.transform(metadata_norm, axis='observation')

    #move Observation Metadata from original to filtered OTU table
    normalized_table = transfer_observation_metadata(otu_table, normalized_table, 'observation')

    make_output_dir_for_file(opts.output_otu_fp)
    write_biom_table(normalized_table, opts.output_otu_fp)


if __name__ == "__main__":
    main()
