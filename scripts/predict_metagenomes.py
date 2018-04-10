#!/usr/bin/env python

from __future__ import division

__license__ = "GPL"
__version__ = "2-alpha.6"


import argparse
from biom import load_table
from picrust2.predict_metagenomes import predict_metagenomes, predict_metagenome_variances,\
  calc_nsti,load_subset_from_biom_str
from picrust2.util import make_output_dir_for_file, write_biom_table
from os import path
from os.path import split,join,splitext
from picrust2.util import get_picrust_project_dir, scale_metagenomes, \
    picrust_formatter
import gzip
import re

parser = argparse.ArgumentParser(

    description="This script produces metagenome functional predictions for " +
                "a sequence abundance table.",

    formatter_class=argparse.RawDescriptionHelpFormatter)


parser.add_argument('-i','--input', metavar='PATH', required=True, type=str,
                    help='The input sequence abundance table in BIOM format')

parser.add_argument('-o','--output', metavar='PATH', required=True, type=str,
                    help='The output table of gene family abundances in ' +
                         'BIOM format')

parser.add_argument('-c', '--input_count', metavar='PATH', required=True,
                    type=str, 
                    help='Precalculated function predictions on per ' +
                         'sequence basis in BIOME format (can be gzipped)')
    
parser.add_argument('-a', '--accuracy_metrics', default=None, metavar='PATH',
                    type=str,
                    help='If provided, calculate weighted NSTI values for ' +
                         'the predicted metagenomes. Note this requires NSTI' +
                         'values to be present in the precalculated file.')

parser.add_argument('--no_round', default=False, action="store_true", 
                    help='Disable rounding number of predicted functions to ' +
                         'the the nearest whole number. This option is ' +
                         'important if you are inputting abundances as ' +
                         'proportions')

parser.add_argument('--normalize_by_function', default=False,
                    action="store_true",
                    help='Normalizes the predicted functional abundances by ' +
                         'dividing each abundance by the sum of functional ' +
                         'abundances in the sample. Total sum of abundances ' +
                         'for each sample will equal 1.')

parser.add_argument('--normalize_by_otu', default=False, action="store_true",
                    help='Normalizes the predicted functional abundances by ' +
                         'dividing each abundance by the sum of sequences '
                         'in the sample. Note: total sum of abundances for ' +
                         'each sample will NOT equal 1.')

parser.add_argument('--suppress_subset_loading', default=False,
                    action="store_true",
                    help='Normally, only counts for sequences present in ' +
                         'the sample are loaded. If this flag is passed, ' +
                         'the full BIOM table is loaded. This makes no ' +
                         'difference for the analysis, but may result ' +
                         'in faster load times (at the cost of more memory ' +
                         'usage)')

parser.add_argument('--load_precalc_file_in_biom', default=False,
                    action="store_true",
                    help='Instead of loading the precalculated file in ' +
                         'tab-delimited format (with sequence ids as row ' +
                         'ids and traits as columns) load the data ' +
                         'in BIOM format (with sequences as SampleIds ' +
                         'and traits as ObservationIds)')

parser.add_argument('-f', '--format_tab_delimited', action="store_true",
                    default=False,
                    help='Output the predicted metagenome table in ' +
                         'tab-delimited format')

parser.add_argument('--verbose', default=False, action="store_true",
                    help='Details will be printed to screen while program ' +
                         'runs')

def determine_data_table_fp(precalc_data_dir,type_of_prediction,gg_version,\
      user_specified_table=None,precalc_file_suffix='precalculated.tab.gz',verbose=False):
    """Determine data table to load, allowing custom user tables or a choice of precalculated files

    precalc_data_dir -- the directory where precalculated tables of gene counts and variances are stored
    type_of_prediction -- a string describing the type of precalculated prediction file.
    gg_version -- the version of greengenes the precalculated prediction was generated against

    This function assumes that precalculated files are named based on the type of prediction,
    (KO, COG, PFAM, etc) and the greengenes version, and then end with a set suffix, which might
    vary between count tables and variance tables.
    """

    if(user_specified_table is None):
        #We assume the precalc file has a specific name (e.g. ko_13_5_precalculated.tab.gz)
        precalc_file_name='_'.join([type_of_prediction,gg_version,\
          precalc_file_suffix])

        input_count_table=join(precalc_data_dir,precalc_file_name)
    else:
        input_count_table=user_specified_table

    if verbose:
        print("Selected data table for loading: " + input_count_table)
    return input_count_table


def main():

    args = parser.parse_args()

    if args.verbose:
        print("Loading OTU table: " + args.input)

    otu_table = load_table(args.input)

    ids_to_load = otu_table.ids(axis='observation').tolist()

    # Determine whether user wants predictions round to nearest whole
    # number or not.
    if args.no_round:
        round_flag = False
    else:
        round_flag = True

    if args.verbose:
        print("Done loading OTU table containing %i samples and %i OTUs." \
          %(len(otu_table.ids()),len(otu_table.ids(axis='observation'))))

    if args.verbose:
        print("Loading gene count data from file: %s" %args.input_count)

    # Load a table of gene counts by OTUs.
    genome_table=load_table(args.input_count).transpose()

    if args.verbose:
        print("Loaded %i genes across %i OTUs from gene count table" \
          %(len(genome_table.ids(axis='observation')),len(genome_table.ids())))

    # if args.with_confidence:
    #     if args.input_variance_table:
    #         variance_table_fp = args.input_variance_table
    #     else:
    #         variance_table_fp = determine_data_table_fp(precalc_data_dir,\
    #           args.type_of_prediction,args.gg_version,\
    #           precalc_file_suffix='precalculated_variances.tab.gz',\
    #           user_specified_table=args.input_count_table)

        if args.verbose:
            print("Loading variance information from table: %s" \
            %variance_table_fp)

        variance_table= load_data_table(variance_table_fp,\
          load_data_table_in_biom=args.load_precalc_file_in_biom,\
          suppress_subset_loading=args.suppress_subset_loading,\
          ids_to_load=ids_to_load,transpose=True)

        if args.verbose:
            print("Loaded %i genes across %i OTUs from variance table" \
              %(len(variance_table.ids(axis='observation')),len(variance_table.ids())))
        #Raise an error if the genome table and variance table differ
        #in the genomes they contain.
        #better to find out now than have something obscure happen latter on
        if args.verbose:
            print("Checking that genome table and variance table are consistent")
        try:
            assert set(variance_table.ids(axis='observation')) == set(genome_table.ids(axis='observation'))
        except AssertionError as e:
            for var_id in variance_table.ids(axis='observation'):
                if var_id not in genome_table.ids(axis='observation'):
                    print("Variance table ObsId %s not in genome_table ObsIds" %var_id)
            raise AssertionError("Variance table and genome table contain different gene ids")
        try:
            assert set(variance_table.ids()) == set(genome_table.ids())
        except AssertionError as e:
            for var_id in variance_table.ids():
                if var_id not in genome_table.ids():
                    print("Variance table SampleId %s not in genome_table SampleIds" %var_id)
            raise AssertionError("Variance table and genome table contain different OTU ids")

        #sort the ObservationIds and SampleIds to be in the same order
        variance_table=variance_table.sort_order(genome_table.ids(axis='observation'), axis='observation')
        variance_table=variance_table.sort_order(genome_table.ids(), axis='sample')

    make_output_dir_for_file(args.output)

    if args.accuracy_metrics:
        # Calculate accuracy metrics
        weighted_nsti = calc_nsti(otu_table,genome_table,weighted=True)
        samples= weighted_nsti[0]
        nstis = list(weighted_nsti[1])
        samples_and_nstis = zip(samples,nstis)
        if args.verbose:
            print ("Writing NSTI information to file:" +  args.accuracy_metrics)
        accuracy_output_fh = open(args.accuracy_metrics,'w')
        accuracy_output_fh.write("#Sample\tMetric\tValue\n")
        for sample,nsti in samples_and_nstis:
            line = "%s\tWeighted NSTI\t%s\n" %(sample,str(nsti))
            accuracy_output_fh.write(line)

    # if args.with_confidence:
    #     #If we are calculating variance, we get the prediction as part
    #     #of the process

    #     if args.verbose:
    #         print "Predicting the metagenome, metagenome variance and confidence intervals for the metagenome..."

    #     predicted_metagenomes,predicted_metagenome_variances,\
    #     predicted_metagenomes_lower_CI_95,predicted_metagenomes_upper_CI_95=\
    #       predict_metagenome_variances(otu_table,genome_table,variance_table,whole_round=round_flag)
    # else:
    #If we don't need confidence intervals, we can do a faster pure numpy prediction

    if args.verbose:
        print("Predicting the metagenome...")
    predicted_metagenomes = predict_metagenomes(otu_table,genome_table,whole_round=round_flag)

    if args.normalize_by_otu:
        #normalize (e.g. divide) the abundances by the sum of the OTUs per sample
        if args.verbose:
            print("Normalizing functional abundances by sum of OTUs per sample")
        inverse_otu_sums = [1/x for x in otu_table.sum(axis='sample')]
        scaling_factors = dict(zip(otu_table.ids(),inverse_otu_sums))
        predicted_metagenomes = scale_metagenomes(predicted_metagenomes,scaling_factors)

    if args.normalize_by_function:
        #normalize (e.g. divide) the abundances by the sum of the functions per sample
        #Sum of functional abundances per sample will equal 1 (e.g. relative abundance).
        if args.verbose:
            print("Normalizing functional abundances by sum of functions per sample")
        predicted_metagenomes = predicted_metagenomes.norm(axis='sample', inplace=False)


    write_metagenome_to_file(predicted_metagenomes,args.output,\
        args.format_tab_delimited,"metagenome prediction",verbose=args.verbose)

    # if args.with_confidence:
    #     output_path,output_filename = split(args.output_metagenome_table)
    #     base_output_filename,ext = splitext(output_filename)
    #     variance_output_fp =\
    #       join(output_path,"%s_variances%s" %(base_output_filename,ext))
    #     upper_CI_95_output_fp =\
    #       join(output_path,"%s_upper_CI_95%s" %(base_output_filename,ext))
    #     lower_CI_95_output_fp =\
    #       join(output_path,"%s_lower_CI_95%s" %(base_output_filename,ext))

    #     write_metagenome_to_file(predicted_metagenome_variances,\
    #       variance_output_fp,args.format_tab_delimited,\
    #       "metagenome prediction variance",verbose=args.verbose)

    #     write_metagenome_to_file(predicted_metagenomes_upper_CI_95,\
    #       upper_CI_95_output_fp,args.format_tab_delimited,\
    #       "metagenome prediction upper 95% confidence interval",\
    #       verbose=args.verbose)

    #     write_metagenome_to_file(predicted_metagenomes_lower_CI_95,\
    #       lower_CI_95_output_fp,args.format_tab_delimited,\
    #       "metagenome prediction lower 95% confidence interval",\
    #       verbose=args.verbose)


def write_metagenome_to_file(predicted_metagenome,output_fp,\
    tab_delimited=False,verbose_filetype_message="metagenome prediction",\
    verbose=False):
    """Write a BIOM Table object to a file, creating the directory if needed
    predicted_metagenome -- a BIOM table object
    output_fp -- the filepath to write the output
    tab_delimited -- if False, write in BIOm format, otherwise write as a tab-delimited file
    verbose -- if True output verbose info to StdOut
    """

    if verbose:
        print("Writing %s results to output file: %s"\
          %(verbose_filetype_message,output_fp))

    make_output_dir_for_file(output_fp)
    if tab_delimited:
        #peek at first observation to decide on what observeration metadata
        #to output in tab-delimited format
        (obs_val,obs_id,obs_metadata)=\
          predicted_metagenome.iter(axis='observation').next()

        #see if there is a metadata field that contains the "Description"
        #(e.g. KEGG_Description or COG_Description)
        h = re.compile('.*Description')
        metadata_names=filter(h.search,obs_metadata.keys())
        if metadata_names:
            #use the "Description" field we found
            metadata_name=metadata_names[0]
        elif(obs_metadata.keys()):
            #if no "Description" metadata then just output the first
            #observation metadata
            metadata_name=(obs_metadata.keys())[0]
        else:
            #if no observation metadata then don't output any
            metadata_name=None

        open(output_fp,'w').write(predicted_metagenome.to_tsv(\
          header_key=metadata_name,header_value=metadata_name,metadata_formatter=biom_meta_to_string))
    else:
        #output in BIOM format
        format_fs = {'KEGG_Description': picrust_formatter,
                     'COG_Description': picrust_formatter,
                     'KEGG_Pathways': picrust_formatter,
                     'COG_Category': picrust_formatter
                     }
        write_biom_table(predicted_metagenome, output_fp, format_fs=format_fs)


def biom_meta_to_string(metadata, replace_str=':'):
    """ Determine which format the metadata is (e.g. str, list, or list of lists) and then convert to a string"""

    #Note that since ';' and '|' are used as seperators we must replace them if they exist
    if isinstance(metadata, (str, unicode)):
        return metadata.replace(';',replace_str)
    elif isinstance(metadata, list):
        if isinstance(metadata[0], list):
            return "|".join(";".join([y.replace(';',replace_str).replace('|',replace_str) for y in x]) for x in metadata)
        else:
            return ";".join(x.replace(';',replace_str) for x in metadata)


if __name__ == "__main__":
    main()
