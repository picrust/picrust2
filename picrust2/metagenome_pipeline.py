#!/usr/bin/env python

__license__ = "GPL"
__version__ = "2.5.1"

import sys
import pandas as pd
from joblib import Parallel, delayed
import numpy as np
from os import path
from picrust2.util import (read_seqabun, make_output_dir, check_files_exist,
                           three_df_index_overlap_sort)


def run_metagenome_pipeline(input_seqabun,
                            function,
                            max_nsti,
                            marker=None,
                            min_reads=1,
                            min_samples=1,
                            strat_out=False,
                            wide_table=False,
                            skip_norm=False,
                            out_dir='metagenome_out'):
    '''Main function to run full metagenome pipeline. Meant to run modular
    functions largely listed below. Will return predicted metagenomes
    straitifed and unstratified by contributing genomes (i.e. taxa).'''

    if not marker and not skip_norm:
        sys.exit("Table of predicted marker gene copy numbers is required "
                 "unless --skip_norm is specified.")
    elif marker and skip_norm:
        sys.exit("Table of predicted marker gene copy numbers should not be "
                 "specified when --skip_norm option is set.")

    make_output_dir(out_dir)

    # Initialize empty pandas dataframe to contain NSTI values.
    nsti_val = pd.DataFrame()

    study_seq_counts = read_seqabun(input_seqabun)

    pred_function = pd.read_csv(function, sep="\t", dtype={'sequence': str})
    pred_function.set_index('sequence', drop=True, inplace=True)

    # If NSTI column present then remove all rows with value above specified
    # max value. Also, remove NSTI column (in both dataframes).
    if 'metadata_NSTI' in pred_function.columns:
        pred_function, nsti_val = drop_tips_by_nsti(tab=pred_function,
                                                    nsti_col='metadata_NSTI',
                                                    max_nsti=max_nsti)
    if not skip_norm:
        check_files_exist([marker])
        pred_marker = pd.read_csv(marker, sep="\t", dtype={'sequence': str})
        pred_marker.set_index('sequence', drop=True, inplace=True)

        if 'metadata_NSTI' in pred_marker.columns:
            pred_marker, nsti_val = drop_tips_by_nsti(tab=pred_marker,
                                                      nsti_col='metadata_NSTI',
                                                      max_nsti=max_nsti)

        # Re-order predicted abundance tables to be in same order as study seqs.
        # Also, drop any sequence ids that don't overlap across all dataframes.
        study_seq_counts, pred_function, pred_marker = three_df_index_overlap_sort(study_seq_counts,
                                                                                   pred_function,
                                                                                   pred_marker)
        norm_output = path.join(out_dir, "seqtab_norm.tsv.gz")

        # Normalize input study sequence abundances by predicted abundance of
        # marker genes and output normalized table if specified.
        study_seq_counts = norm_by_marker_copies(input_seq_counts=study_seq_counts,
                                                 input_marker_num=pred_marker,
                                                 norm_filename=norm_output)
    else:
        # Get intersecting rows between input files and sort.
        label_overlap = pred_function.index.intersection(study_seq_counts.index).sort_values()

        if len(label_overlap) == 0:
            sys.exit("No sequence ids overlap between both input files.")

        pred_function = pred_function.reindex(label_overlap)
        study_seq_counts = study_seq_counts.reindex(label_overlap)

    # If NSTI column input then output weighted NSTI values.
    if not nsti_val.empty:
        weighted_nsti_out = path.join(out_dir, "weighted_nsti.tsv.gz")
        calc_weighted_nsti(seq_counts=study_seq_counts,
                           nsti_input=nsti_val,
                           outfile=weighted_nsti_out)

    # Determine which sequences should be in the "RARE" category if stratified
    # table is specified.
    if strat_out:
        rare_seqs = []

        if min_reads != 1 or min_samples != 1:
            rare_seqs = id_rare_seqs(in_counts=study_seq_counts,
                                     min_reads=min_reads,
                                     min_samples=min_samples)

    # Generate and return final tables.
    if not strat_out:
        return(None, unstrat_funcs_only_by_samples(pred_function,
                                                   study_seq_counts))

    elif strat_out and not wide_table:
        return(metagenome_contributions(pred_function, study_seq_counts,
                                        rare_seqs),
               unstrat_funcs_only_by_samples(pred_function, study_seq_counts))

    elif strat_out and wide_table:
        return(strat_funcs_by_samples(pred_function, study_seq_counts,
                                      rare_seqs))


def strat_funcs_by_samples(func_abun, sample_abun, rare_seqs=[],
                           return_unstrat=True):
    '''Take in function table and study sequence abundance table. Returns
    stratified table of function abundances per sequence (rows) per samples
    (columns). Will also return unstratified format (function by sample).
    Will collapse rare sequences into a single sequence category if
    given a non-empty list of ids.'''

    # Create combined stratified dataframe.
    strat_func = pd.concat((sample_abun.multiply(func_abun[func], axis=0) for func in func_abun.columns),
                           ignore_index=True, axis=0, sort=True)

    # Set multi-index labels to be sequence and function ids.
    strat_func.index = pd.MultiIndex.from_product((func_abun.columns,
                                                   func_abun.index))
    strat_func.index.names = ['function', 'sequence']

    if len(rare_seqs) > 0:
        # Sum rare rows and make a new row named "RARE" in each df.
        # Remove all the original rare rows.

        # First identify rows to remove, slice them out, and remove from original.
        rare_seqs_rows = strat_func.index.get_level_values('sequence').isin(rare_seqs)
        raw_seqs_slice = strat_func.iloc[rare_seqs_rows]
        strat_func = strat_func.iloc[~rare_seqs_rows]

        # Next remove sequence index and sum based on function id. Re-add "RARE"
        # to be the sequence id.
        raw_seqs_slice.index = raw_seqs_slice.index.droplevel('sequence')
        raw_seqs_slice = raw_seqs_slice.sum(level=['function'])
        raw_seqs_slice['sequence'] = 'RARE'
        raw_seqs_slice.set_index('sequence', append=True, drop=True,
                                 inplace=True)

        # Concat the RARE seqs to the full df.
        strat_func = pd.concat([strat_func, raw_seqs_slice], axis=0, sort=True)

    # Remove rows that are all 0 and keep sample ids in same order.
    strat_func = strat_func.loc[~(strat_func == 0).all(axis=1)]

    strat_func = strat_func[list(sample_abun.columns)]

    # Return dataframe and also unstratified dataframe if specified.
    if return_unstrat:
        return(strat_func, strat_func.groupby(level='function', axis=0).sum())
    else:
        return(strat_func)


def unstrat_funcs_only_by_samples(func_abun, sample_abun):
    '''Take in function table and study sequence abundance table. Returns
    unstratified table of function abundances by samples.'''

    # Sample ids are taken from sequence abundance table.
    sample_ids = sample_abun.columns.values

    # List that will contain a series for each sample's unstratified
    # function abundnaces.
    sample_funcs = []

    for sample in sample_ids:
        sample_funcs.append(func_abun.mul(sample_abun[sample], axis=0).sum(axis=0))

    # Build dataframe from these series.
    unstrat_func = pd.concat(sample_funcs, axis=1, sort=True)

    unstrat_func = unstrat_func.loc[~(unstrat_func == 0).all(axis=1)]

    unstrat_func.columns = sample_ids

    unstrat_func.index.name = 'function'

    return(unstrat_func)


def drop_tips_by_nsti(tab, nsti_col, max_nsti):
    '''Takes in (1) function table (columns are functions, rows are predicted
    genomes for ASVs, (2) which column of this function table corresponds to
    NSTI values, and (3) the maximum NSTI value for ASVs to be retained. Will
    return table with ASVs above the max NSTI cut-off excluded and a dataframe
    which just contains the NSTI values for ASVs that passed the cut-off. Will
    report how many ASVs were removed and will throw error if all ASVs are.'''

    orig_num_rows = tab.shape[0]

    tab = tab[tab[nsti_col] <= max_nsti]

    filt_num_rows = tab.shape[0]

    if orig_num_rows == filt_num_rows:
        print("All ASVs were below the max NSTI cut-off of " + str(max_nsti) +
              " and so all were retained for downstream analyses.",
              file=sys.stderr)

    elif filt_num_rows == 0:
        sys.exit("Stopping - all ASVs filtered from table when max NSTI "
                 "cut-off of " + str(max_nsti) + " used.")

    else:
        num_removed = orig_num_rows - filt_num_rows
        print(str(num_removed) + " of " + str(orig_num_rows) + " ASVs were "
              "above the max NSTI cut-off of " + str(max_nsti) + " and were "
              "removed from the downstream analyses.", file=sys.stderr)

    # Keep track of NSTI column as separate dataframe and remove this column
    # from the main dataframe.
    nsti_val = tab[[nsti_col]]

    return(tab.drop(nsti_col, axis=1, inplace=False), nsti_val)


def calc_weighted_nsti(seq_counts, nsti_input, outfile=None, return_df=False):
    '''Will calculate weighted NSTI values given sequence count table and NSTI
    value for each sequence. Will output these weighted values to a file if
    output file is specified. Will only return a df if specified.'''

    nsti_mult = seq_counts.mul(nsti_input.metadata_NSTI, axis=0)

    # Get column sums divided by total abundance per sample.
    weighted_nsti = pd.DataFrame(nsti_mult.sum(axis=0)/seq_counts.sum(axis=0))

    weighted_nsti.columns = ["weighted_NSTI"]

    weighted_nsti.fillna(0, inplace=True)

    # Write to outfile if specified.
    if outfile:
        weighted_nsti.to_csv(path_or_buf=outfile, sep="\t", header=True,
                             index_label="sample", compression="infer")
    if return_df:
        return(weighted_nsti)


def norm_by_marker_copies(input_seq_counts,
                          input_marker_num,
                          norm_filename=None,
                          round_decimal=2):

    '''Divides sequence counts (which correspond to amplicon sequence
    variants) by the predicted marker gene copies for each sequence. Will write
    out the normalized table if option specified.'''

    input_seq_counts = input_seq_counts.div(input_marker_num.loc[
                                                input_seq_counts.index.values,
                                                input_marker_num.columns.values[0]],
                                            axis="index")

    input_seq_counts = input_seq_counts.round(decimals=round_decimal)

    # Output normalized table if specified.
    if norm_filename:
        input_seq_counts.to_csv(path_or_buf=norm_filename,
                                index_label="normalized",
                                sep="\t",
                                compression="infer")

    return(input_seq_counts)


def id_rare_seqs(in_counts, min_reads, min_samples):
    '''Determine which rows of a sequence countfile are below either the
    cut-offs of min read counts or min samples present.'''

    # Check if "RARE" is the name of a sequence in this table.
    if "RARE" in in_counts.index:
        sys.exit("Stopping: the sequence called \"RARE\" in the sequence " +
                 "abundance table should be re-named.")

    low_freq_seq = set(in_counts[in_counts.sum(axis=1) < min_reads].index)
    few_samples_seq = set(in_counts[(in_counts != 0).astype(int).sum(axis=1) < min_samples].index)

    return(list(low_freq_seq.union(few_samples_seq)))


def metagenome_contributions(func_abun, sample_abun, rare_seqs=[],
                             skip_abun=False):
    '''Take in function table and study sequence abundance table. Returns
    long-form table of how each sequence contributes functions in each
    sample. Note that the old format of columns (such as calling the sequences
    "OTUs" is retained here for backwards compatability. A subset of input
    sequences will be collapsed to a single category called "RARE" if a
    non-empty list is input for the rare_seqs option. The skip_abun option
    can be set when the abundances columns are not needed.'''

    # Make copy of sample abundance that is in terms of relative abundances.
    sample_relabun = sample_abun.div(sample_abun.sum(axis=0), axis=1) * 100

    # Counter used to identify the first sample.
    s_i = 0

    for sample in sample_abun.columns:

        single_abun = sample_abun[sample]
        single_abun = single_abun.iloc[single_abun.to_numpy().nonzero()]

        single_relabun = sample_relabun[sample]
        single_relabun = single_relabun.iloc[single_relabun.to_numpy().nonzero()]

        intersecting_taxa = single_abun.index.intersection(func_abun.index)
        func_abun_subset = func_abun.loc[intersecting_taxa]
        single_abun = single_abun.loc[intersecting_taxa]
        single_relabun = single_relabun.loc[intersecting_taxa]

        # Melt function table to be long format.
        func_abun_subset['taxon'] = func_abun_subset.index.to_list()

        func_abun_subset_melt = pd.melt(func_abun_subset,
                                        id_vars='taxon',
                                        value_name='genome_function_count',
                                        var_name='function')

        # Remove rows where gene count is 0.
        func_abun_subset_melt = func_abun_subset_melt[func_abun_subset_melt['genome_function_count'] != 0]

        if not skip_abun:

            func_abun_subset_melt['taxon_abun'] = single_abun.loc[func_abun_subset_melt['taxon'].to_list()].to_list()

            func_abun_subset_melt['taxon_rel_abun'] = single_relabun.loc[func_abun_subset_melt['taxon'].to_list()].to_list()

            func_abun_subset_melt['taxon_function_abun'] = func_abun_subset_melt['genome_function_count'] * func_abun_subset_melt['taxon_abun']

            func_abun_subset_melt['taxon_rel_function_abun'] = func_abun_subset_melt['genome_function_count'] * func_abun_subset_melt['taxon_rel_abun']

            func_abun_subset_melt['norm_taxon_function_contrib'] = func_abun_subset_melt['taxon_function_abun'] / \
                                                                   func_abun_subset_melt.groupby("function").sum(numeric_only = True)["taxon_function_abun"][func_abun_subset_melt["function"]].to_numpy()

        # Collapse sequences identified as "rare" to the same category.
        rare_seqs = [r for r in rare_seqs if r in func_abun_subset_melt['taxon'].to_list()]

        if len(rare_seqs) > 0:
            func_abun_subset_melt.loc[func_abun_subset_melt['taxon'].isin(rare_seqs), 'taxon'] = 'RARE'
            func_abun_subset_melt = func_abun_subset_melt.groupby(['function', 'taxon'],
                                                                  as_index=False).sum()

        func_abun_subset_melt['sample'] = sample

        # Order column names.
        if skip_abun:

            func_abun_subset_melt = func_abun_subset_melt[['sample',
                                                           'function',
                                                           'taxon',
                                                           'genome_function_count']]
        else:
            func_abun_subset_melt = func_abun_subset_melt[['sample',
                                                           'function',
                                                           'taxon',
                                                           'taxon_abun',
                                                           'taxon_rel_abun',
                                                           'genome_function_count',
                                                           'taxon_function_abun',
                                                           'taxon_rel_function_abun',
                                                           'norm_taxon_function_contrib']]

            # Make sure there are no NaN values in final column.
            if func_abun_subset_melt['norm_taxon_function_contrib'].isna().sum() > 0:
                sys.exit("Error - NaN values are present in the norm_taxon_function_contrib column, which indicates that the calculation failed.")


        if s_i == 0:
            metagenome_contrib_out = func_abun_subset_melt.copy()
            s_i += 1
        else:
            metagenome_contrib_out = pd.concat([metagenome_contrib_out,
                                                func_abun_subset_melt],
                                                axis=0)

    return(metagenome_contrib_out)


def contrib_to_unstrat(contrib_table, sample_order=None):
    '''Take in metagenome contribution table and return wide-format sample
    (unstratified) metagenome table.'''

    # Return empty datafame.
    if contrib_table.shape[0] == 0:
        return(pd.DataFrame())

    contrib_table = contrib_table[['sample', 'function', 'taxon_function_abun']]

    contrib_table = pd.pivot_table(data=contrib_table, columns='sample',
                                   index='function',
                                   values='taxon_function_abun',
                                   aggfunc=np.sum, fill_value=0)

    contrib_table.index.name = None
    contrib_table.columns.name = None

    if sample_order:
        contrib_table = contrib_table.reindex(columns=sample_order,
                                              fill_value=0)

    return(contrib_table)
