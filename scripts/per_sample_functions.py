#!/usr/bin/env python

from __future__ import division

__license__ = "GPL"
__version__ = "2-alpha.7"

import argparse
import biom
import pandas as pd
import numpy as np
from rpy2.robjects import pandas2ri
from joblib import Parallel, delayed
from picrust2.per_sample_functions import(read_in_rds,
                                          norm_by_marker_copies,
                                          process_func_count_prob)

parser = argparse.ArgumentParser(

    description="Outputs expected counts and confidence intervals per " +
                "function for each sample.",

    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-i', '--input', metavar='PATH', required=True, type=str,
                    help='Input table in BIOM format of study sequences ' +
                         'abundances')

parser.add_argument('-r', '--rds', metavar='PATH', required=True, type=str,
                    help='Input RDS file of predicted function states per ' +
                         'study sequence')

parser.add_argument('-m', '--marker', metavar='PATH', required=True, type=str,
                    help='File containing predicted number of marker gene ' +
                         'copies per study sequence, which is used to ' +
                         'normalize table of sequence counts')

parser.add_argument('--exp_outfile', default='exp_predict_out.tsv', 
                    metavar='PATH', type=str,
                    help='Output file for predicted function abundances ' +
                         'per sample')

parser.add_argument('--ci_outfile', default='ci_predict_out.tsv', 
                    metavar='PATH', type=str,
                    help='Output file for confidence intervals of predicted ' +
                         'function abundances per sample')

parser.add_argument('-t', '--threads', default=1, type=int,
                    help='Number of threads')

parser.add_argument('--output_norm', default=False, action='store_true',
                    help='If specified, write out sequence abundance table ' +
                         'normalized (divided) by predicted abundance of ' +
                         'marker genes')

parser.add_argument('--norm_outfile', default='norm_seq_counts.tsv', 
                    metavar='PATH', type=str,
                    help='Optional output file of normalized sequence ' +
                         'abundance counts (only if \"--output_norm\" flag ' +
                         'set')


def main():

    args = parser.parse_args()

    pandas2ri.activate()

    predict_func_probs = read_in_rds(args.rds)

    func_names = predict_func_probs.names

    # Read in study sequences in HDF5 BIOM format.
    input_biom = biom.load_table(args.input)

    # Convert biom table to pandas dataframe.
    # (Based on James Morton's blog post here:
    # http://mortonjt.blogspot.ca/2016/07/behind-scenes-with-biom-tables.html)
    study_seq_counts = pd.DataFrame(np.array(
                                        input_biom.matrix_data.todense()),
                             index=input_biom.ids(axis='observation'),
                             columns=input_biom.ids(axis='sample'))

    exp_marker_copy = pd.read_table(filepath_or_buffer=args.marker,
                                    sep="\t",
                                    index_col="tips")

    study_seq_counts = norm_by_marker_copies(study_seq_counts,
                                             exp_marker_copy,
                                             args.output_norm,
                                             args.norm_outfile)

    if args.threads > 1:
        metagenome_out = Parallel(n_jobs=args.threads)(delayed(
                                  process_func_count_prob)(
                                  predict_func_probs[i],
                                  func_names[i],
                                  study_seq_counts)
                                  for i in range(len(predict_func_probs)))

    else:
        metagenome_out = []
        tmp = [0, 237]
        for i in tmp: #range(len(predict_func_probs)):
            metagenome_out += [process_func_count_prob(predict_func_probs[i],
                                                       func_names[i],
                                                       study_seq_counts)]

    # Convert returned list of dictionaries to dataframe.
    metagenome_out_df = pd.DataFrame.from_dict(metagenome_out)

    # Set index names to be functions.
    metagenome_out_df.set_index("function", inplace=True)

    # Split dataframe into expectations and CIs separately.
    metagenome_out_df_exp = metagenome_out_df.filter(items=study_seq_counts.columns.values)
    metagenome_out_df_ci = metagenome_out_df.filter(regex=("_ci\d+$"))


    # Write output files.
    metagenome_out_df_exp.to_csv(path_or_buf=args.exp_outfile,
                                 index_label="function",
                                 sep="\t")

    metagenome_out_df_ci.to_csv(path_or_buf=args.ci_outfile,
                            index_label="function",
                            sep="\t")

if __name__ == "__main__":
    main()
