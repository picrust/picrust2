#!/usr/bin/env python
# File created on 22 Feb 2012
from __future__ import division

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2011, The PICRUST project"
__credits__ = ["Greg Caporaso","Jesse Zaneveld","Morgan Langille"]
__license__ = "GPL"
__version__ = "2-alpha.1"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"


from cogent.util.option_parsing import parse_command_line_parameters, make_option
from biom import load_table
from picrust.predict_metagenomes import predict_metagenomes, calc_nsti
from picrust.metagenome_contributions import partition_metagenome_contributions
from picrust.util import make_output_dir_for_file, get_picrust_project_dir, convert_precalc_to_biom
from os import path
from os.path import join
import gzip
from sys import exit

script_info = {}
script_info['brief_description'] = "This script partitions metagenome functional contributions according to function, OTU, and sample, for a given OTU table."
script_info['script_description'] = ""
script_info['script_usage'] = [
("","Partition the predicted contribution to the  metagenomes from each organism in the given OTU table, limited to only K00001, K00002, and K00004.","%prog -i normalized_otus.biom -l K00001,K00002,K00004 -o ko_metagenome_contributions.tab"),
("","Partition the predicted contribution to the  metagenomes from each organism in the given OTU table, limited to only COG0001 and COG0002.","%prog -i normalized_otus.biom -l COG0001,COG0002 -t cog -o cog_metagenome_contributions.tab")
]
script_info['output_description']= "Output is a tab-delimited column indicating OTU contribution to each function."
script_info['required_options'] = [
 make_option('-i','--input_otu_table',type='existing_filepath',help='the input otu table in biom format'),
 make_option('-o','--output_fp',type="new_filepath",help='the output file for the metagenome contributions')
]
type_of_prediction_choices=['ko','cog','rfam']
gg_version_choices=['13_5','18may2012']
script_info['optional_options'] = [\
    make_option('-t','--type_of_prediction',default=type_of_prediction_choices[0],type="choice",\
                    choices=type_of_prediction_choices,\
                    help='Type of functional predictions. Valid choices are: '+\
                    ', '.join(type_of_prediction_choices)+\
                    ' [default: %default]'),
    make_option('-g','--gg_version',default=gg_version_choices[0],type="choice",\
                    choices=gg_version_choices,\
                    help='Version of GreenGenes that was used for OTU picking. Valid choices are: '+\
                    ', '.join(gg_version_choices)+\
                    ' [default: %default]'),

    make_option('-c','--input_count_table',default=None,type="existing_filepath",help='Precalculated function predictions on per otu basis in biom format (can be gzipped). Note: using this option overrides --type_of_prediction and --gg_version. [default: %default]'),
 make_option('--suppress_subset_loading',default=False,action="store_true",help='Normally, only counts for OTUs present in the sample are loaded.  If this flag is passed, the full biom table is loaded.  This makes no difference for the analysis, but may result in faster load times (at the cost of more memory usage)'),
    make_option('--load_precalc_file_in_biom',default=False,action="store_true",help='Instead of loading the precalculated file in tab-delimited format (with otu ids as row ids and traits as columns) load the data in biom format (with otu as SampleIds and traits as ObservationIds) [default: %default]'),
    make_option('-f','--limit_to_functional_categories',default=False,action="store",type='string',help='If provided only output prediction for functions that match the specified functional category. Multiple categories can be passed as a list separated by | [default: %default]'),
        make_option('-l','--limit_to_function',default=None,help='If provided, only output predictions for the specified function ids.  Multiple function ids can be passed using comma delimiters.')
]
script_info['version'] = __version__


def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)


    if opts.limit_to_function:
        limit_to_functions = opts.limit_to_function.split(',')
        if opts.verbose:
            print "Limiting output to only functions:",limit_to_functions
    else:
        limit_to_functions = []

    if opts.verbose:
        print "Loading otu table: ",opts.input_otu_table

    otu_table = load_table(opts.input_otu_table)
    ids_to_load = otu_table.ids(axis='observation')

    if(opts.input_count_table is None):
        #precalc file has specific name (e.g. ko_13_5_precalculated.tab.gz)
        precalc_file_name='_'.join([opts.type_of_prediction,opts.gg_version,'precalculated.tab.gz'])
        input_count_table=join(get_picrust_project_dir(),'picrust','data',precalc_file_name)
    else:
        input_count_table=opts.input_count_table

    if opts.verbose:
        print "Loading trait table: ", input_count_table

    ext=path.splitext(input_count_table)[1]

    if opts.verbose:
        print "Loading count table: ", input_count_table

    if (ext == '.gz'):
        genome_table_fh = gzip.open(input_count_table,'rb')
    else:
        genome_table_fh = open(input_count_table,'U')

    #In the genome/trait table genomes are the samples and
    #genes are the observations


    if opts.load_precalc_file_in_biom:
        if not opts.suppress_subset_loading:
            #Now we want to use the OTU table information
            #to load only rows in the count table corresponding
            #to relevant OTUs

            if opts.verbose:
                print "Loading traits for %i organisms from the trait table" %len(ids_to_load)

            genome_table = load_subset_from_biom_str(genome_table_fh.read(),ids_to_load,axis='samples')
        else:
            if opts.verbose:
                print "Loading *full* count table because --suppress_subset_loading was passed. This may result in high memory usage"
            genome_table = load_table(genome_table_fh)
    else:
        genome_table = convert_precalc_to_biom(genome_table_fh,ids_to_load)
    ok_functional_categories = None

    metadata_type = None
    if opts.limit_to_functional_categories:
        ok_functional_categories = opts.limit_to_functional_categories.split("|")
        if opts.verbose:
            print "Limiting to functional categories: %s" %(str(ok_functional_categories))

        # Either KEGG_Pathways or COG_Category needs
        # to be assigned to metadata_key to limit to
        # functional categories (not needed for 
        # individual functions) 

        if opts.type_of_prediction == "ko":
            metadata_type = "KEGG_Pathways"
        elif opts.type_of_prediction == "cog":
            metadata_type = "COG_Category"
        elif opts.type_of_prediction == "rfam":
            exit("Stopping program: when type of prediction is set to rfam you can only limit to individual functions (-l) rather than to functional categories (-f)")
              
    partitioned_metagenomes = partition_metagenome_contributions(otu_table,genome_table,limit_to_functions=limit_to_functions,\
      limit_to_functional_categories = ok_functional_categories ,  metadata_key = metadata_type )

    output_text = "\n".join(["\t".join(map(str,i)) for i in partitioned_metagenomes])
    if opts.verbose:
        print "Writing results to output file: ",opts.output_fp

    make_output_dir_for_file(opts.output_fp)
    open(opts.output_fp,'w').write(output_text)

if __name__ == "__main__":
    main()
