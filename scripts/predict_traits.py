#!/usr/bin/env python
# File created on 15 Jul 2011
from __future__ import division

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
__credits__ = ["Jesse RR Zaneveld", "Morgan Langille"]
__license__ = "GPL"
__version__ = "2-alpha.1"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"


from math import e
from os.path import splitext
from cogent.util.option_parsing import parse_command_line_parameters, make_option
from picrust.parse import extract_ids_from_table,\
  parse_asr_confidence_output
from picrust.predict_traits import assign_traits_to_tree,\
  predict_traits_from_ancestors, update_trait_dict_from_file,\
  make_neg_exponential_weight_fn, biom_table_from_predictions,\
  predict_random_neighbor,predict_nearest_neighbor,\
  calc_nearest_sequenced_taxon_index,weighted_average_tip_prediction, \
  get_brownian_motion_param_from_confidence_intervals
from picrust.util import make_output_dir_for_file, write_biom_table, convert_biom_to_precalc
from picrust.format_tree_and_trait_table import load_picrust_tree, set_label_conversion_fns

script_info = {}
script_info['brief_description'] = "Given a tree and a set of known character states (observed traits and reconstructions), output predictions for unobserved character states"
script_info['script_description'] =\
"""
This script produces predictions of unobserved traits given a phylogenetic tree and a table that summarizes which traits are present in which ancestral organisms.
In the most common usage, this script is used to predict which gene families are present in each OTU (Operational Taxonomic Unit; roughly equivalent to a bacterial 'species'), given a tree and a set of ancestral state reconstructions.

The output of the script is a trait prediction file, which summarizes the predicted traits of each organism of interest (by default, this is all of the organisms that are tips in the phylogenetic tree).

The prediction method works as follows:

    1.  For each terminal (tip) node where a prediction is to be performed, the algorithm through the reconstructed ancestral states, and finds the last node in the ancestry of our organism of interest for which a prediction is available

    2.  The trait for the organism is then predicted based on a branch-length weighted average of the ancestral node and it's close relatives. (This is necessary because technical limitations involving the handling of ambiguous characters in many Ancestral State Reconstruction programs prevent the parent node of the organism from being directly reconstructed in most cases.)

    The exact weight function to use can be specified from the commandline (see options below).

    In general, this approach causes the prediction to be a weighted average of the closest reconstructed ancestor, and the either reconstructed or directly observed trait value of the organism of interest's sibling node(s).
"""

#Define valid choices for 'choice parameters
METHOD_CHOICES = ['asr_and_weighting','nearest_neighbor','asr_only','weighting_only',\
  'random_neighbor']
WEIGHTING_CHOICES = ['exponential','linear','equal']
CONFIDENCE_FORMAT_CHOICES = ['sigma','confidence_interval']

#Add script information
script_info['script_usage'] = [\
("","Required options with NSTI:","%prog -a -i trait_table.tab -t reference_tree.newick -r asr_counts.tab -o predict_traits.tab"),\
("","Limit predictions to particular tips in OTU table:","%prog -a -i trait_table.tab -t reference_tree.newick -r asr_counts.tab -o predict_traits_limited.tab -l otu_table.tab"),
("","Reconstruct confidence","%prog -a -i trait_table.tab -t reference_tree.newick -r asr_counts.tab -c asr_ci.tab -o predict_traits.tab")
]
#Define commandline interface
script_info['output_description']= "Output is a table (tab-delimited or .biom) of predicted character states"
script_info['required_options'] = [\
make_option('-i','--observed_trait_table',type="existing_filepath",\
  help='the input trait table describing directly observed traits (e.g. sequenced genomes) in tab-delimited format'),\
make_option('-t','--tree',type="existing_filepath",\
  help='the full reference tree, in Newick format')
]
script_info['optional_options'] = [\
 make_option('-o','--output_trait_table',type="new_filepath",\
   default='predicted_states.tsv',help='the output filepath for trait predictions [default: %default]'),\
 make_option('-a','--calculate_accuracy_metrics',default=False,action="store_true",\
   help='if specified, calculate accuracy metrics (i.e. how accurate does PICRUSt expect its predictions to be?) and add to output file [default: %default]'),\
 make_option('--output_accuracy_metrics_only',type="new_filepath",\
   default=None,help='if specified, calculate accuracy metrics (e.g. NSTI), output them to this filepath, and do not do anything else. [default: %default]'),\

 make_option('-m','--prediction_method',default='asr_and_weighting',choices=METHOD_CHOICES,help='Specify prediction method to use.  The recommended prediction method is set as default, so other options are primarily useful for control experiments and methods validation, not typical use.  Valid choices are:'+",".join(METHOD_CHOICES)+'.  "asr_and_weighting"(recommended): use ancestral state reconstructions plus local weighting with known tip nodes.  "nearest_neighbor": predict the closest tip on the tree with trait information.  "random_annotated_neighbor": predict a random tip on the tree with trait information. "asr_only": predict the traits of the last reconstructed ancestor, without weighting. "weighting_only": weight all genomes by distance, to the organism of interest using the specified weighting function and predict the weighted average.   [default: %default]'),\

 make_option('-w','--weighting_method',default='exponential',choices=WEIGHTING_CHOICES,help='Specify prediction the weighting function to use.  This only applies to prediction methods that incorporate local weighting ("asr_and_weighting" or "weighting_only")  The recommended weighting  method is set as default, so other options are primarily useful for control experiments and methods validation, not typical use.  Valid choices are:'+",".join(WEIGHTING_CHOICES)+'.  "exponential"(recommended): weight genomes as a negative exponent of distance.  That is 2^-d, where d is the tip-to-tip distance from the genome to the tip.  "linear": weight tips as a linear function of weight, normalized to the maximum possible distance (max_d -d)/d. "equal_weights": set all weights to a constant (ignoring branch length).   [default: %default]'),
 make_option('-l','--limit_predictions_by_otu_table',type="existing_filepath",help='Specify a valid path to a legacy QIIME OTU table to perform predictions only for tips that are listed in the OTU table (regardless of abundance)'),\
 make_option('-g','--limit_predictions_to_organisms',help='Limit predictions to specific, comma-separated organims ids. (Generally only useful for lists of < 10 organism ids, for example when performing leave-one-out cross-validation).'),\
 make_option('-r','--reconstructed_trait_table',\
   type="existing_filepath",default=None,\
   help='the input trait table describing reconstructed traits (from ancestral_state_reconstruction.py) in tab-delimited format [default: %default]'),\

 make_option('--confidence_format',\
   choices=CONFIDENCE_FORMAT_CHOICES,default='sigma',\
   help='the format for the confidence intervals from ancestral state reconstruction. Only needed if passing a reconstruction confidence file with -c or --reconstruction_confidence.  These are typically sigma values for maximum likelihood ASR  methods, but 95% confidence intervals for phylogenetic independent contrasts (e.g. from the ape R packages ace function with pic as the reconstruction method).  Valid choices are:'+",".join(CONFIDENCE_FORMAT_CHOICES)+'. [default: %default]'),\
 make_option('-c','--reconstruction_confidence',\
   type="existing_filepath",default=None,\
   help='the input trait table describing confidence intervals for reconstructed traits (from ancestral_state_reconstruction.py) in tab-delimited format [default: %default]'),
   make_option('--output_precalc_file_in_biom',default=False,action="store_true",help='Instead of outputting the precalculated file in tab-delimited format (with otu ids as row ids and traits as columns) output the data in biom format (with otu as SampleIds and traits as ObservationIds) [default: %default]'),
   make_option('--no_round',default=False,action="store_true",help='Flag to set if you do not want predictions to be rounded to the nearest integer [default: %default]')
]

script_info['version'] = __version__

#Helper formatting functions

def write_results_to_file(f_out,headers,predictions,sep="\t"):
    """Write a set of predictions to a file

    headers -- a list of header columns
    predictions -- a dict of predictions, keyed by organism
    with arrays as values
    """
    f= open(f_out,"w")
    lines = [sep.join(headers)+"\n"]

    for pred in sorted(predictions.keys()):
        new_fields = [pred]
        value = predictions[pred]
        #print value
        if value is None or len(value) == 0:
            continue
        trait_value_fields = list(value)
        new_fields.extend(trait_value_fields)
        #print new_fields
        new_line = sep.join(map(str,new_fields))
        lines.append("%s\n" % new_line)
    f.writelines(lines)
    f.close()

#Main script

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)

    #if we specify we want NSTI only then we have to calculate it first
    if opts.output_accuracy_metrics_only:
        opts.calculate_accuracy_metrics=True

    if opts.verbose:
        print "Loading tree from file:", opts.tree

    if opts.no_round:
        round_opt = False 
    else:
        round_opt = True

    # Load Tree
    tree = load_picrust_tree(opts.tree, opts.verbose)

    table_headers=[]
    traits={}
    #load the asr trait table using the previous list of functions to order the arrays
    if opts.reconstructed_trait_table:
        table_headers,traits =\
                update_trait_dict_from_file(opts.reconstructed_trait_table)

        #Only load confidence intervals on the reconstruction
        #If we actually have ASR values in the analysis
        if opts.reconstruction_confidence:
            if opts.verbose:
                print "Loading ASR confidence data from file:",\
                opts.reconstruction_confidence
                print "Assuming confidence data is of type:",opts.confidence_format

            asr_confidence_output = open(opts.reconstruction_confidence)
            asr_min_vals,asr_max_vals, params,column_mapping =\
              parse_asr_confidence_output(asr_confidence_output,format=opts.confidence_format)
            if 'sigma' in params:
                brownian_motion_parameter = params['sigma'][0]
            else:
                brownian_motion_parameter = None

            if opts.verbose:
                print "Done. Loaded %i confidence interval values." %(len(asr_max_vals))
                print "Brownian motion parameter:",brownian_motion_parameter
        else:
            brownian_motion_parameter = None

    #load the trait table into a dict with organism names as keys and arrays as functions
    table_headers,genome_traits =\
            update_trait_dict_from_file(opts.observed_trait_table,table_headers)


    #Combine the trait tables overwriting the asr ones if they exist in the genome trait table.
    traits.update(genome_traits)

    # Specify the attribute where we'll store the reconstructions
    trait_label = "Reconstruction"

    if opts.verbose:
        print "Assigning traits to tree..."

    # Decorate tree using the traits
    tree = assign_traits_to_tree(traits,tree, trait_label=trait_label)


    if opts.reconstruction_confidence:
        if opts.verbose:
            print "Assigning trait confidence intervals to tree..."
        tree = assign_traits_to_tree(asr_min_vals,tree,\
            trait_label="lower_bound")

        tree = assign_traits_to_tree(asr_max_vals,tree,\
            trait_label="upper_bound")

        if brownian_motion_parameter is None:

             if opts.verbose:
                 print "No Brownian motion parameters loaded. Inferring these from 95% confidence intervals..."
             brownian_motion_parameter = get_brownian_motion_param_from_confidence_intervals(tree,\
                      upper_bound_trait_label="upper_bound",\
                      lower_bound_trait_label="lower_bound",\
                      trait_label=trait_label,\
                      confidence=0.95)
             if opts.verbose:
                 print "Inferred the following rate parameters:",brownian_motion_parameter
    if opts.verbose:
        print "Collecting list of nodes to predict..."

    #Start by predict all tip nodes.
    nodes_to_predict = [tip.Name for tip in tree.tips()]

    if opts.verbose:
        print "Found %i nodes to predict." % len(nodes_to_predict)

    if opts.limit_predictions_to_organisms:
        organism_id_str = opts.limit_predictions_to_organisms
        ok_organism_ids = organism_id_str.split(',')
        ok_organism_ids = [n.strip() for n in ok_organism_ids]
        for f in set_label_conversion_fns(True,True):
            ok_organism_ids = [f(i) for i in ok_organism_ids]

        if opts.verbose:
            print "Limiting predictions to user-specified ids:",\
              ",".join(ok_organism_ids)


        if not ok_organism_ids:
            raise RuntimeError(\
              "Found no valid ids in input: %s. Were comma-separated ids specified on the command line?"\
              % opts.limit_predictions_to_organisms)

        nodes_to_predict =\
          [n for n in nodes_to_predict if n in ok_organism_ids]

        if not nodes_to_predict:
            raise RuntimeError(\
              "Filtering by user-specified ids resulted in an empty set of nodes to predict.   Are the ids on the commmand-line and tree ids in the same format?  Example tree tip name: %s, example OTU id name: %s" %([tip.Name for tip in tree.tips()][0],ok_organism_ids[0]))

        if opts.verbose:
            print "After filtering organisms to predict by the ids specified on the commandline, %i nodes remain to be predicted" %(len(nodes_to_predict))

    if opts.limit_predictions_by_otu_table:
        if opts.verbose:
            print "Limiting predictions to ids in user-specified OTU table:",\
              opts.limit_predictions_by_otu_table
        otu_table = open(opts.limit_predictions_by_otu_table,"U")
        #Parse OTU table for ids

        otu_ids =\
          extract_ids_from_table(otu_table.readlines(),delimiter="\t")

        if not otu_ids:
            raise RuntimeError(\
              "Found no valid ids in input OTU table: %s.  Is the path correct?"\
              % opts.limit_predictions_by_otu_table)

        nodes_to_predict =\
          [n for n in nodes_to_predict if n in otu_ids]

        if not nodes_to_predict:
            raise RuntimeError(\
              "Filtering by OTU table resulted in an empty set of nodes to predict.   Are the OTU ids and tree ids in the same format?  Example tree tip name: %s, example OTU id name: %s" %([tip.Name for tip in tree.tips()][0],otu_ids[0]))

        if opts.verbose:
            print "After filtering by OTU table, %i nodes remain to be predicted" %(len(nodes_to_predict))

    # Calculate accuracy of PICRUST for the given tree, sequenced genomes
    # and set of ndoes to predict
    accuracy_metrics = ['NSTI']
    accuracy_metric_results = None
    if opts.calculate_accuracy_metrics:
        if opts.verbose:
            print "Calculating accuracy metrics: %s" %([",".join(accuracy_metrics)])
        accuracy_metric_results = {}
        if 'NSTI' in accuracy_metrics:

            nsti_result,min_distances =\
                calc_nearest_sequenced_taxon_index(tree,\
                limit_to_tips = nodes_to_predict,\
                trait_label = trait_label, verbose=opts.verbose)

            #accuracy_metric_results['NSTI'] = nsti_result
            for organism in min_distances.keys():
                accuracy_metric_results[organism] = {'NSTI': min_distances[organism]}

            if opts.verbose:
                print "NSTI:", nsti_result

        if opts.output_accuracy_metrics_only:
            #Write accuracy metrics to file
            if opts.verbose:
                print "Writing accuracy metrics to file:",opts.output_accuracy_metrics

            f = open(opts.output_accuracy_metrics_only,'w+')
            f.write("metric\torganism\tvalue\n")
            lines =[]
            for organism in accuracy_metric_results.keys():
                for metric in accuracy_metric_results[organism].keys():
                    lines.append('\t'.join([metric,organism,\
                      str(accuracy_metric_results[organism][metric])])+'\n')
            f.writelines(sorted(lines))
            f.close()
            exit()


    if opts.verbose:
        print "Generating predictions using method:",opts.prediction_method

    if opts.weighting_method == 'exponential':
        #For now, use exponential weighting
        weight_fn = make_neg_exponential_weight_fn(e)

    variances=None #Overwritten by methods that calc variance
    confidence_intervals=None #Overwritten by methods that calc variance

    if opts.prediction_method == 'asr_and_weighting':
        # Perform predictions using reconstructed ancestral states

        if opts.reconstruction_confidence:
            predictions,variances,confidence_intervals =\
              predict_traits_from_ancestors(tree,nodes_to_predict,\
              trait_label=trait_label,\
              lower_bound_trait_label="lower_bound",\
              upper_bound_trait_label="upper_bound",\
              calc_confidence_intervals = True,\
              brownian_motion_parameter=brownian_motion_parameter,\
              weight_fn=weight_fn,verbose=opts.verbose,
              round_predictions=round_opt)

        else:
             predictions =\
              predict_traits_from_ancestors(tree,nodes_to_predict,\
              trait_label=trait_label,\
              weight_fn =weight_fn,verbose=opts.verbose,
              round_predictions=round_opt)

    elif opts.prediction_method == 'weighting_only':
        #Ignore ancestral information
        predictions =\
          weighted_average_tip_prediction(tree,nodes_to_predict,\
          trait_label=trait_label,\
          weight_fn =weight_fn,verbose=opts.verbose)



    elif opts.prediction_method == 'nearest_neighbor':

        predictions = predict_nearest_neighbor(tree,nodes_to_predict,\
          trait_label=trait_label,tips_only = True)

    elif opts.prediction_method == 'random_neighbor':

        predictions = predict_random_neighbor(tree,\
          nodes_to_predict,trait_label=trait_label)

    if opts.verbose:
        print "Done making predictions."

    make_output_dir_for_file(opts.output_trait_table)

    out_fh=open(opts.output_trait_table,'w')
    #Generate the table of biom predictions
    if opts.verbose:
        print "Converting results to .biom format for output..."

    biom_predictions=biom_table_from_predictions(predictions,table_headers,\
                                                         observation_metadata=None,\
                                                         sample_metadata=accuracy_metric_results,convert_to_int=False)
    if opts.verbose:
        print "Writing prediction results to file: ",opts.output_trait_table

    if opts.output_precalc_file_in_biom:

        #write biom table to file
        write_biom_table(biom_predictions, opts.output_trait_table)

    else:
        #convert to precalc (tab-delimited) format

        out_fh = open(opts.output_trait_table, 'w')
        out_fh.write(convert_biom_to_precalc(biom_predictions))
        out_fh.close()

    #Write out variance information to file
    if variances:

        if opts.verbose:
            print "Converting variances to BIOM format"

        if opts.output_precalc_file_in_biom:
            suffix='.biom'
        else:
            suffix='.tab'

        biom_prediction_variances=biom_table_from_predictions({k:v['variance'] for k,v in variances.iteritems()},table_headers,\
        observation_metadata=None,\
        sample_metadata=None,convert_to_int=False)
        outfile_base,extension = splitext(opts.output_trait_table)
        variance_outfile = outfile_base+"_variances"+suffix
        make_output_dir_for_file(variance_outfile)

        if opts.verbose:
            print "Writing variance information to file:",variance_outfile

        if opts.output_precalc_file_in_biom:
            write_biom_table(biom_prediction_variances, variance_outfile)
        else:
            open(variance_outfile,'w').write(\
                convert_biom_to_precalc(biom_prediction_variances))


    if confidence_intervals:

        if opts.verbose:
            print "Converting upper confidence interval values to BIOM format"

        biom_prediction_upper_CI=biom_table_from_predictions({k:v['upper_CI'] for k,v in confidence_intervals.iteritems()},table_headers,\
          observation_metadata=None,\
          sample_metadata=None,convert_to_int=False)

        outfile_base,extension = splitext(opts.output_trait_table)
        upper_CI_outfile = outfile_base+"_upper_CI"+suffix
        make_output_dir_for_file(upper_CI_outfile)

        if opts.verbose:
            print "Writing upper confidence limit information to file:",upper_CI_outfile

        if opts.output_precalc_file_in_biom:
            write_biom_table(biom_prediction_upper_CI, upper_CI_outfile)
        else:
            open(upper_CI_outfile,'w').write(\
                convert_biom_to_precalc(biom_prediction_upper_CI))

        biom_prediction_lower_CI=biom_table_from_predictions({k:v['lower_CI'] for k,v in confidence_intervals.iteritems()},table_headers,\
          observation_metadata=None,\
          sample_metadata=None,convert_to_int=False)

        outfile_base,extension = splitext(opts.output_trait_table)
        lower_CI_outfile = outfile_base+"_lower_CI"+suffix
        make_output_dir_for_file(lower_CI_outfile)

        if opts.verbose:
            print "Writing lower confidence limit information to file",lower_CI_outfile

        if opts.output_precalc_file_in_biom:
            write_biom_table(biom_prediction_lower_CI, lower_CI_outfile)
        else:
            open(lower_CI_outfile,'w').write(\
                convert_biom_to_precalc(biom_prediction_lower_CI))


if __name__ == "__main__":
    main()
