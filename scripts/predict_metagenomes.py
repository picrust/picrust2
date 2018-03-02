#!/usr/bin/env python
# File created on 22 Feb 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
__credits__ = ["Greg Caporaso","Jesse Zaneveld","Morgan Langille"]
__license__ = "GPL"
__version__ = "2-alpha.1"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"


from cogent.util.option_parsing import parse_command_line_parameters, make_option
from biom import load_table
from picrust.predict_metagenomes import predict_metagenomes,predict_metagenome_variances,\
  calc_nsti,load_subset_from_biom_str
from picrust.util import make_output_dir_for_file,write_biom_table, convert_precalc_to_biom
from os import path
from os.path import split,join,splitext
from picrust.util import get_picrust_project_dir, scale_metagenomes, \
    picrust_formatter
import gzip
import re

script_info = {}
script_info['brief_description'] = "This script produces the actual metagenome functional predictions for a given OTU table."
script_info['script_description'] = ""
script_info['script_usage'] = [
                               ("","Predict KO abundances for a given OTU table picked against the newest version of GreenGenes.", "%prog -i normalized_otus.biom -o predicted_metagenomes.biom"),
                               ("","Change output format to plain tab-delimited:","%prog -f -i normalized_otus.biom -o predicted_metagenomes.txt"),
                               ("","Predict KO abundances for a given OTU table, but do not round predictions to nearest whole number (esp. needed for proportional abundances):","%prog --no_round -i normalized_otus.biom -o predicted_metagenomes.txt"),
                               ("","Predict COG abundances for a given OTU table.","%prog -i normalized_otus.biom -t cog -o cog_predicted_metagenomes.biom"),\
                               ("","Predict RFAM abundances for a given OTU table.","%prog -i normalized_otus.biom -t rfam -o rfam_predicted_metagenomes.biom"),\
                               ("","Output confidence intervals for each prediction.","%prog -i normalized_otus.biom -o predicted_metagenomes.biom --with_confidence"),\
                               ("","Predict metagenomes using a custom trait table in tab-delimited format.","%prog -i otu_table_for_custom_trait_table.biom -c custom_trait_table.tab -o output_metagenome_from_custom_trait_table.biom"),\
                               ("","Predict metagenomes,variances,and 95% confidence intervals for each gene category using a custom trait table in tab-delimited format.","%prog -i otu_table_for_custom_trait_table.biom --input_variance_table custom_trait_table_variances.tab -c custom_trait_table.tab -o output_metagenome_from_custom_trait_table.biom --with_confidence"),\
                               ("","Change the version of GG used to pick OTUs","%prog -i normalized_otus.biom -g 18may2012 -o predicted_metagenomes.biom")
                                ]
script_info['output_description']= "Output is a table of function counts (e.g. KEGG KOs) by sample ids."
script_info['required_options'] = [
 make_option('-i','--input_otu_table',type='existing_filepath',help='the input otu table in biom format'),
 make_option('-o','--output_metagenome_table',type="new_filepath",help='the output file for the predicted metagenome')
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
    make_option('-a','--accuracy_metrics',default=None,type="new_filepath",help='If provided, calculate accuracy metrics for the predicted metagenome.  NOTE: requires that per-genome accuracy metrics were calculated using predict_traits.py during genome prediction (e.g. there are "NSTI" values in the genome .biom file metadata)'),
    make_option('--no_round',default=False,action="store_true",help='Disable rounding number of predicted functions to the the nearest whole number. This option is important if you are inputting abundances as proportions [default: %default]'),
    make_option('--normalize_by_function',default=False,action="store_true",help='Normalizes the predicted functional abundances by dividing each abundance by the sum of functional abundances in the sample. Total sum of abundances for each sample will equal 1.'),
    make_option('--normalize_by_otu',default=False,action="store_true",help='Normalizes the predicted functional abundances by dividing each abundance by the sum of OTUs in the sample. Note: total sum of abundances for each sample will NOT equal 1.'),
    make_option('--suppress_subset_loading',default=False,action="store_true",help='Normally, only counts for OTUs present in the sample are loaded.  If this flag is passed, the full biom table is loaded.  This makes no difference for the analysis, but may result in faster load times (at the cost of more memory usage)'),
    make_option('--load_precalc_file_in_biom',default=False,action="store_true",help='Instead of loading the precalculated file in tab-delimited format (with otu ids as row ids and traits as columns) load the data in biom format (with otu as SampleIds and traits as ObservationIds) [default: %default]'),
    make_option('--input_variance_table',default=None,type="existing_filepath",help='Precalculated table of variances corresponding to the precalculated table of function predictions.  As with the count table, these are on a per otu basis and in BIOM format (can be gzipped). Note: using this option overrides --type_of_prediction and --gg_version. [default: %default]'),
    make_option('--with_confidence',default=False,action="store_true",help='Calculate 95% confidence intervals for metagenome predictions.  By default, this uses the confidence intervals for the precalculated table of genes for greengenes OTUs.  If you pass a custom count table with -c and select this option, you must also specify a corresponding table of confidence intervals for the gene content prediction using --input_variance_table. (these are generated by running predict_traits.py with the --with_confidence option). If this flag is set, three addtional output files will be generated, named the same as the metagenome prediction output, but with .variance .upper_CI or .lower_CI appended immediately before the file extension [default: %default]'),
    make_option('-f','--format_tab_delimited',action="store_true",default=False,help='output the predicted metagenome table in tab-delimited format [default: %default]')]
script_info['version'] = __version__


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
        print "Selected data table for loading: ", input_count_table
    return input_count_table


def load_data_table(data_table_fp,\
  load_data_table_in_biom=False,suppress_subset_loading=False,ids_to_load=None,\
  transpose=False,verbose=False):
    """Load a data table, detecting gziiped files and subset loading
    data_table_fp -- path to the input data table

    load_data_table_in_biom -- if True, load the data table as a BIOM table rather
    than as tab-delimited

    suppress_subset_loading -- if True, load the entire table, rather than just
    ids_of_interest

    ids_to_load -- a list of OTU ids for which data should be loaded

    gzipped files are detected based on the '.gz' suffix.
    """
    if not path.exists(data_table_fp):
        raise IOError("File "+data_table_fp+" doesn't exist! Did you forget to download it?")

    ext=path.splitext(data_table_fp)[1]
    if (ext == '.gz'):
        genome_table_fh = gzip.open(data_table_fp,'rb')
    else:
        genome_table_fh = open(data_table_fp,'U')

    if load_data_table_in_biom:
        if not suppress_subset_loading:
            #Now we want to use the OTU table information
            #to load only rows in the count table corresponding
            #to relevant OTUs

            if verbose:
                print "Loading traits for %i organisms from the trait table" %len(ids_to_load)

            genome_table = load_subset_from_biom_str(genome_table_fh.read(),ids_to_load,axis='samples')
        else:
            if verbose:
                print "Loading *full* count table because --suppress_subset_loading was passed. This may result in high memory usage"
            genome_table = load_table(data_table_fp)
    else:
        genome_table = convert_precalc_to_biom(genome_table_fh,ids_to_load,transpose=transpose)

    if verbose:
        print "Done loading trait table containing %i functions for %i organisms." %(len(genome_table.ids(axis='observation')),len(genome_table.ids()))

    return genome_table


def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)

    if opts.verbose:
        print "Loading OTU table: ",opts.input_otu_table

    otu_table = load_table(opts.input_otu_table)
    ids_to_load = otu_table.ids(axis='observation').tolist()

    # Determine whether user wants predictions round to nearest whole
    # number or not.
    if opts.no_round:
        round_flag = False
    else:
        round_flag = True

    if opts.verbose:
        print "Done loading OTU table containing %i samples and %i OTUs." \
          %(len(otu_table.ids()),len(otu_table.ids(axis='observation')))

    #Hardcoded loaction of the precalculated datasets for PICRUSt,
    #relative to the project directory
    precalc_data_dir=join(get_picrust_project_dir(),'picrust','data')

    # Load a table of gene counts by OTUs.
    #This can be either user-specified or precalculated
    genome_table_fp = determine_data_table_fp(precalc_data_dir,\
      opts.type_of_prediction,opts.gg_version,\
      user_specified_table=opts.input_count_table,verbose=opts.verbose)

    if opts.verbose:
        print "Loading gene count data from file: %s" %genome_table_fp

    genome_table= load_data_table(genome_table_fp,\
      load_data_table_in_biom=opts.load_precalc_file_in_biom,\
      suppress_subset_loading=opts.suppress_subset_loading,\
      ids_to_load=ids_to_load,verbose=opts.verbose,transpose=True)

    if opts.verbose:
        print "Loaded %i genes across %i OTUs from gene count table" \
          %(len(genome_table.ids(axis='observation')),len(genome_table.ids()))

    if opts.with_confidence:
        if opts.input_variance_table:
            variance_table_fp = opts.input_variance_table
        else:
            variance_table_fp = determine_data_table_fp(precalc_data_dir,\
              opts.type_of_prediction,opts.gg_version,\
              precalc_file_suffix='precalculated_variances.tab.gz',\
              user_specified_table=opts.input_count_table)

        if opts.verbose:
            print "Loading variance information from table: %s" \
            %variance_table_fp

        variance_table= load_data_table(variance_table_fp,\
          load_data_table_in_biom=opts.load_precalc_file_in_biom,\
          suppress_subset_loading=opts.suppress_subset_loading,\
          ids_to_load=ids_to_load,transpose=True)

        if opts.verbose:
            print "Loaded %i genes across %i OTUs from variance table" \
              %(len(variance_table.ids(axis='observation')),len(variance_table.ids()))
        #Raise an error if the genome table and variance table differ
        #in the genomes they contain.
        #better to find out now than have something obscure happen latter on
        if opts.verbose:
            print "Checking that genome table and variance table are consistent"
        try:
            assert set(variance_table.ids(axis='observation')) == set(genome_table.ids(axis='observation'))
        except AssertionError,e:
            for var_id in variance_table.ids(axis='observation'):
                if var_id not in genome_table.ids(axis='observation'):
                    print "Variance table ObsId %s not in genome_table ObsIds" %var_id
            raise AssertionError("Variance table and genome table contain different gene ids")
        try:
            assert set(variance_table.ids()) == set(genome_table.ids())
        except AssertionError,e:
            for var_id in variance_table.ids():
                if var_id not in genome_table.ids():
                    print "Variance table SampleId %s not in genome_table SampleIds" %var_id
            raise AssertionError("Variance table and genome table contain different OTU ids")

        #sort the ObservationIds and SampleIds to be in the same order
        variance_table=variance_table.sort_order(genome_table.ids(axis='observation'), axis='observation')
        variance_table=variance_table.sort_order(genome_table.ids(), axis='sample')

    make_output_dir_for_file(opts.output_metagenome_table)

    if opts.accuracy_metrics:
        # Calculate accuracy metrics
        weighted_nsti = calc_nsti(otu_table,genome_table,weighted=True)
        samples= weighted_nsti[0]
        nstis = list(weighted_nsti[1])
        samples_and_nstis = zip(samples,nstis)
        if opts.verbose:
            print "Writing NSTI information to file:", opts.accuracy_metrics
        accuracy_output_fh = open(opts.accuracy_metrics,'w')
        accuracy_output_fh.write("#Sample\tMetric\tValue\n")
        for sample,nsti in samples_and_nstis:
            line = "%s\tWeighted NSTI\t%s\n" %(sample,str(nsti))
            accuracy_output_fh.write(line)

    if opts.with_confidence:
        #If we are calculating variance, we get the prediction as part
        #of the process

        if opts.verbose:
            print "Predicting the metagenome, metagenome variance and confidence intervals for the metagenome..."

        predicted_metagenomes,predicted_metagenome_variances,\
        predicted_metagenomes_lower_CI_95,predicted_metagenomes_upper_CI_95=\
          predict_metagenome_variances(otu_table,genome_table,variance_table,whole_round=round_flag)
    else:
        #If we don't need confidence intervals, we can do a faster pure numpy prediction

        if opts.verbose:
            print "Predicting the metagenome..."
        predicted_metagenomes = predict_metagenomes(otu_table,genome_table,whole_round=round_flag)

    if opts.normalize_by_otu:
        #normalize (e.g. divide) the abundances by the sum of the OTUs per sample
        if opts.verbose:
            print "Normalizing functional abundances by sum of OTUs per sample"
        inverse_otu_sums = [1/x for x in otu_table.sum(axis='sample')]
        scaling_factors = dict(zip(otu_table.ids(),inverse_otu_sums))
        predicted_metagenomes = scale_metagenomes(predicted_metagenomes,scaling_factors)

    if opts.normalize_by_function:
        #normalize (e.g. divide) the abundances by the sum of the functions per sample
        #Sum of functional abundances per sample will equal 1 (e.g. relative abundance).
        if opts.verbose:
            print "Normalizing functional abundances by sum of functions per sample"
        predicted_metagenomes = predicted_metagenomes.norm(axis='sample', inplace=False)


    write_metagenome_to_file(predicted_metagenomes,opts.output_metagenome_table,\
        opts.format_tab_delimited,"metagenome prediction",verbose=opts.verbose)

    if opts.with_confidence:
        output_path,output_filename = split(opts.output_metagenome_table)
        base_output_filename,ext = splitext(output_filename)
        variance_output_fp =\
          join(output_path,"%s_variances%s" %(base_output_filename,ext))
        upper_CI_95_output_fp =\
          join(output_path,"%s_upper_CI_95%s" %(base_output_filename,ext))
        lower_CI_95_output_fp =\
          join(output_path,"%s_lower_CI_95%s" %(base_output_filename,ext))

        write_metagenome_to_file(predicted_metagenome_variances,\
          variance_output_fp,opts.format_tab_delimited,\
          "metagenome prediction variance",verbose=opts.verbose)

        write_metagenome_to_file(predicted_metagenomes_upper_CI_95,\
          upper_CI_95_output_fp,opts.format_tab_delimited,\
          "metagenome prediction upper 95% confidence interval",\
          verbose=opts.verbose)

        write_metagenome_to_file(predicted_metagenomes_lower_CI_95,\
          lower_CI_95_output_fp,opts.format_tab_delimited,\
          "metagenome prediction lower 95% confidence interval",\
          verbose=opts.verbose)


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
        print "Writing %s results to output file: %s"\
          %(verbose_filetype_message,output_fp)

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
