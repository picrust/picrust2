#!/usr/bin/env python
# File created on 10 April 2012
from __future__ import division

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
__credits__ = ["Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "1.1.3"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"

from collections import defaultdict
from os import listdir
from os.path import join
from cogent.util.option_parsing import parse_command_line_parameters,\
  make_option
from picrust.evaluate_test_datasets import unzip,evaluate_test_dataset,\
 update_pooled_data, run_accuracy_calculations_on_biom_table,run_accuracy_calculations_on_pooled_data,\
 format_scatter_data, format_correlation_data, run_and_format_roc_analysis

from biom import load_table
from picrust.util import make_output_dir


script_info = {}
script_info['brief_description'] = "Evaluate the accuracy of character predictions, given directories of expected vs. observed test results"
script_info['script_description'] =\
    """The script finds all paired expected and observed values in a set of directories and generates the following output: 1) data for a scatterplot of observed vs. expected values for each character (typically gene family count) within each organism (so one file per organism). 2) A summary of accuracy across all organisms.
    character """
script_info['script_usage'] = [("","Evaluate the accuracy of all predictions in a folder, and output summary statistics.","%prog -i obs_otu_table.biom -e exp_otu_table.txt -o./evaluation_results/")]
script_info['output_description']= "Outputs will be obs,exp data points for the comparison"
script_info['required_options'] = [
 make_option('-i','--trait_table_dir',type="existing_dirpath",help='the input trait table directory (files in biom format)'),\
 make_option('-e','--exp_trait_table_dir',type="existing_dirpath",help='the input expected trait table directory (files in biom format)'),\
 make_option('-o','--output_dir',type="new_dirpath",help='the output directory'),
]
script_info['optional_options'] = [
        make_option('-f','--field_order',\
                default='file_type,prediction_method,weighting_method,holdout_method,distance,organism',help='pass comma-separated categories, in the order they appear in file names.   Categories are "file_type","prediction_method","weighting_method","holdout_method" (randomization vs. holdout),"distance",and "organism".  Example:  "-f file_type,test_method,asr_method specifies that files will be in the form: predict_traits--distance_exclusion--wagner.  Any unspecified values are set to "not_specified".  [default: %default]'),\
        make_option('-p','--pool_by',\
          default='distance',help='pass comma-separated categories to pool results by those metadata categories. Valid categories are: holdout_method, prediction_method,weighting_method,distance and organism. For example, pass "distance" to output results pooled by holdout distance in addition to holdout method and prediction method  [default: %default]')
]
script_info['version'] = __version__



def evaluate_test_dataset_dir(obs_dir_fp,exp_dir_fp,file_name_delimiter="--",\
        file_name_field_order=\
        {'file_type':0,"prediction_method":1,"weighting_method":2,"holdout_method":3,\
          "distance":4,"organism":5},strict=False, verbose=True,pool_by=['distance'],\
          roc_success_criteria=['binary','exact']):
    """Return control evaluation results from the given directories

    obs_dir_fp -- directory containing PICRUST-predicted genomes.   These MUST start with
    'predict_traits', and must contain the values specified in file_name_field_order,\
    separated by the delimiter given in file_name_delimiter.  For example:

    predict_traits--exclude_tips_by_distance--0.87--'NC_000913|646311926'

    exp_dir_fp -- as obs_dir_fp above, but expectation file names (usually sequenced genomes
    with known gene content) must start with exp_biom_traits

    file_name_delimiter -- the delimiter that separates metadata stored in the filename

    NOTE: technically this isn't the best way of doing things.  We may want at some point
    to revisit this setup and store metadata about each comparison in a separate file.  But
    storing in the filename is convenient for our initial analysis.

    file_name_field_order -- the order of the required metadata fields in the filename.
    Required fields are file_type,method,distance,and organism

    pool_by -- if passed, concatenate traits from each trial that is identical in this category.  e.g. pool_by 'distance' will pool traits across individual test genomes with the same holdout distance.

    roc_success_criteria -- a list of methods for measuring 'success' of a prediction.  Separate  ROC curves will be created for each.
    Description:

    The method assumes that for each file type in the observed directory, a paired file
    is also found in the exp_dir with similar method, distance, and organism, but a varied
    file type (test_tree, test_trait_table)


    Process:
    1. Search test directory for all gene predictions in the correct format
    2. For each, find the corresponding expected trait table in the expectation file
    3. Organize all of the correlations and scatter points
    4. Return the following outputs:
        --- Table of obs,exp values with organism, method, distance metadata
        --- Summary of correlations (table by organism and distance)
        --- Summary of AUC (table by organism,method,distance)

    """
    trials = defaultdict(list)
    correlation_lines = []
    scatter_lines = []
    #We'll want a quick unzip fn for converting points to trials
    #TODO: separate out into a 'get_paired_data_from_dirs' function

    pooled_observations = {}
    pooled_expectations = {}
    input_files=sorted(listdir(obs_dir_fp))
    for file_number,f in enumerate(input_files):
        if verbose:
            print "\nExamining file {0} of {1}: {2}".format(file_number+1,len(input_files),f)
        if 'accuracy_metrics' in f:
            print "%s is an Accuracy file...skipping" %str(f)
            continue
        #filename_components_list = f.split(file_name_delimiter)
        filename_components = {}
        for i,field in enumerate(f.split(file_name_delimiter)):
            filename_components[i]=field
        #if verbose:
        #    print "Filename components:",filename_components
        try:
            file_type,holdout_method,weighting_method,\
            prediction_method,distance,organism = \
              filename_components.get(file_name_field_order.get('file_type','not_specified'),'not_specified'),\
              filename_components.get(file_name_field_order.get('holdout_method','not_specified'),'not_specified'),\
              filename_components.get(file_name_field_order.get('weighting_method','not_specified'),'not_specified'),\
              filename_components.get(file_name_field_order.get('prediction_method','not_specified'),'not_specified'),\
              filename_components.get(file_name_field_order.get('distance','not_specified'),'not_specified'),\
              filename_components.get(file_name_field_order.get('organism','not_specified'),'not_specified')

        except IndexError, e:
            print "Could not parse filename %s using delimiter: %s.  Skipping..." %(f,file_name_delimiter)
            continue

        #Get predicted traits
        if file_type == 'predict_traits':
            if verbose:
                #print "Found a prediction file"
                print "\tLoading .biom format observation table:",f

            try:
              obs_table =\
                load_table(join(obs_dir_fp,f))
              obs_table = obs_table.transpose()
            except ValueError:
                print 'Failed, skipping...'
                continue
        else:
            continue

        # Get paired observation file
        exp_filename = file_name_delimiter.join(['exp_biom_traits',holdout_method,distance,organism])
        exp_filepath = join(exp_dir_fp,exp_filename)
        if verbose:
            print "\tLooking for the expected trait file matching %s here: %s" %(f,exp_filepath)

        try:
            exp_table = load_table(exp_filepath)
        except IOError, e:
            if strict:
                raise IOError(e)
            else:
                if verbose:
                    print "Missing expectation file....skipping!"
                continue
        base_tag =  '%s\t%s\t' %(holdout_method,prediction_method)
        tags = [base_tag+'all_results']
        combined_tag = base_tag +\
                "\t".join([str(field)+"_"+str(filename_components[file_name_field_order[field]]) for field in pool_by])
        tags.append(combined_tag)

        #TODO: abstract out pooling into its own function
        non_pooled_fields = [filename_components.get(file_name_field_order[k],None) for k in file_name_field_order.keys() if k not in pool_by]
        pooled_observations,pooled_expectations =\
                update_pooled_data(obs_table,exp_table,tags,pooled_observations,\
          pooled_expectations,str(file_number),verbose=verbose)

    return run_accuracy_calculations_on_pooled_data(pooled_observations,\
      pooled_expectations,roc_success_criteria=roc_success_criteria,verbose=verbose)

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    pool_by = opts.pool_by.split(',')


    #create output directory
    make_output_dir(opts.output_dir)

    #Construct a dict from user specified field order
    file_name_field_order = {}
    for i,field in enumerate(opts.field_order.split(',')):
        file_name_field_order[field]=i
        if opts.verbose:
            print "Assuming file names are in this order:",file_name_field_order

    for k in pool_by:
        #Check that we're only pooling by values that exist
        if k not in file_name_field_order.keys():
            err_text=\
              "Bad value for option '--pool_by'.  Can't pool by '%s'.   Valid categories are: %s" %(k,\
              ",".join(file_name_field_order.keys()))
            raise ValueError(err_text)

    if opts.verbose:
        print "Pooling results by:",pool_by

    roc_success_criteria = ['binary','exact','int_exact']

    scatter_lines,correlation_lines,roc_result_lines,roc_auc_lines =\
      evaluate_test_dataset_dir(opts.trait_table_dir,\
      opts.exp_trait_table_dir,file_name_delimiter="--",\
      file_name_field_order=file_name_field_order,pool_by=pool_by,\
      roc_success_criteria=roc_success_criteria,verbose=opts.verbose)

    #Output scatter data

    output_fp = join(opts.output_dir,'evaluation_scatter_data.tab')
    if opts.verbose:
        print "Writing scatter plot data to:",output_fp
    file_lines = scatter_lines

    f = open(output_fp,"w+")
    f.writelines(file_lines)
    f.close()

    #Output correlation data

    output_fp = join(opts.output_dir,'evaluation_correlation_data.tab')

    if opts.verbose:
        print "Writing correlation data to:",output_fp

    file_lines = correlation_lines

    f = open(output_fp,"w+")
    f.writelines(file_lines)
    f.close()

    #Output raw ROC plot data
    if opts.verbose:
        print "Writing ROC data..."
    for c in roc_result_lines.keys():
        output_fp = join(opts.output_dir,'evaluation_roc_data_%s.tab' %c)
        if opts.verbose:
            print "Outputting ROC data for success criterion %s to: %s" %(c,output_fp)
        file_lines = roc_result_lines[c]

        f = open(output_fp,"w+")
        f.writelines(file_lines)
        f.close()

    #Output summary ROC AUC data
    if opts.verbose:
        print "Writing ROC AUC data..."

    for c in roc_auc_lines.keys():
        output_fp = join(opts.output_dir,'evaluation_roc_auc_data_%s.tab' %c)
        file_lines = roc_auc_lines[c]

        if opts.verbose:
            print "Outputting ROC AUC data for success criterion %s to: %s" %(c,output_fp)
        f = open(output_fp,"w+")
        f.writelines(file_lines)
        f.close()

if __name__ == "__main__":
    main()
