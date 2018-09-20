#!/usr/bin/env python

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.0.0-b.8"

import argparse
from os import path
import sys
from tempfile import TemporaryDirectory
from picrust2.place_seqs import place_seqs_pipeline
from picrust2.wrap_hsp import castor_hsp_workflow
from picrust2.default import (default_fasta, default_tree, default_tables,
                              default_map, default_regroup_map,
                              default_pathway_map)
from picrust2.run_minpath import run_minpath_pipeline
from picrust2.metagenome_pipeline import run_metagenome_pipeline
from picrust2.util import (make_output_dir, check_files_exist,
                           add_descrip_col)

def full_pipeline(study_fasta,
                  input_table,
                  output_folder,
                  threads,
                  ref_msa,
                  tree,
                  in_traits,
                  custom_trait_tables,
                  marker_gene_table,
                  pathway_map,
                  no_pathways,
                  regroup_map,
                  no_regroup,
                  stratified,
                  max_nsti,
                  min_reads,
                  min_samples,
                  hsp_method,
                  calculate_NSTI,
                  confidence,
                  seed,
                  no_gap_fill,
                  per_sequence_contrib,
                  no_descrip,
                  verbose):
    '''Function that contains wrapper commands for full PICRUSt2 pipeline.
    Descriptions of all of these input arguments/options are given in the
    picrust2_pipeline.py script.'''

    # Check that input files exist.
    check_files_exist([study_fasta, input_table])

    # Make output folder.
    make_output_dir(output_folder)

    out_tree = path.join(output_folder, "out.tre")

    if custom_trait_tables is None:

        # Check that specified functional categories are allowed.
        FUNC_TRAIT_OPTIONS = ['COG', 'EC', 'KO', 'PFAM', 'TIGRFAM']
        funcs = in_traits.split(",")
        for func in funcs:
            if func not in FUNC_TRAIT_OPTIONS:
                sys.exit("Error - specified category " + func + " is not one of "
                         "the default categories.")

        # Add EC to this set if pathways are to be predicted.
        if "EC" not in funcs and not no_pathways:
            funcs.append("EC")

        rxn_func = "EC"

        func_tables = default_tables

    else:
        funcs = []
        func_tables = {}

        table_i = 0

        for custom in custom_trait_tables.split(","):

            func_id = path.splitext(path.basename(custom))[0]
            funcs.append(func_id)
            func_tables[func_id] = custom

            if table_i == 0:
                rxn_func = func_id
                table_i += 1

    # Append marker as well, since this also needs to be run.
    funcs.append("marker")
    func_tables["marker"] = marker_gene_table

    # Methods for discrete trait prediction with CI enabled.
    discrete_set = set(['emp_prob', 'mp'])

    if confidence and hsp_method in discrete_set:
        ci_setting = True
    else:
        ci_setting = False

    gap_fill_opt = not no_gap_fill

    with TemporaryDirectory() as temp_dir:

        if verbose:
            print("Placing sequences onto reference tree.")

        place_seqs_pipeline(study_fasta=study_fasta,
                            ref_msa=ref_msa,
                            tree=tree,
                            out_tree=out_tree,
                            threads=threads,
                            papara_output=None,
                            out_dir=temp_dir,
                            chunk_size=5000,
                            print_cmds=verbose)

        if verbose:
            print("Finished placing sequences on output tree: " + out_tree)

    # Get predictions for all specified functions and keep track of outfiles.
    predicted_funcs = {}

    for func in funcs:

        # Only output NSTI in 16S table.
        nsti_setting = False
        if func == "marker" and calculate_NSTI:
            nsti_setting = True

        if verbose:
            print("Running hidden-state prediction for " + func)

        hsp_table, ci_table = castor_hsp_workflow(tree_path=out_tree,
                                                  trait_table_path=func_tables[func],
                                                  hsp_method=hsp_method,
                                                  calc_nsti=nsti_setting,
                                                  calc_ci=ci_setting,
                                                  check_input=False,
                                                  num_proc=threads,
                                                  ran_seed=seed)

        count_outfile = path.join(output_folder, func + "_predicted.tsv")

        # Add "_nsti" to filename if output.
        if nsti_setting:
            count_outfile = path.join(output_folder, func + "_nsti_predicted.tsv")

        # Keep track of output file name for next step of pipeline.
        predicted_funcs[func] = count_outfile

        if verbose:
            print("Writing out predicted gene family abundances to " + count_outfile)

        hsp_table.to_csv(path_or_buf=count_outfile, index_label="sequence", sep="\t")

        # Output the CI file as well if option set.
        if ci_setting:
            ci_outfile = path.join(output_folder, func + "_predicted_ci.tsv")

            if verbose:
                print("Writing out predicted gene family CIs to " + ci_outfile)

            ci_table.to_csv(path_or_buf=ci_outfile, index_label="sequence",
                            sep="\t")

    marker_infile = predicted_funcs["marker"]

    # Inititalize dictionary of function names to pandas tables to return.
    func_dfs = {}

    # Each value will be a list of 2 elements corresponding to the unstratified
    # and stratified tables respectively (stratified will be None of not calculated).

    # Loop over each function again and run metagenome pipeline.
    for func in funcs:

        if func == "marker":
            continue

        func_infile = predicted_funcs[func]

        if verbose:
            print("Running metagenome pipeline for " + func)

        func_output_dir = path.join(output_folder, func + "_metagenome_out")

        # Infer metagenome abundances per-sample.
        with TemporaryDirectory() as temp_dir:

            # Pass arguments to key function and get predicted functions
            # stratified and unstratified by genomes.
            strat_pred, unstrat_pred = run_metagenome_pipeline(input_biom=input_table,
                                                               function=func_infile,
                                                               marker=marker_infile,
                                                               out_dir=func_output_dir,
                                                               max_nsti=max_nsti,
                                                               min_reads=min_reads,
                                                               min_samples=min_samples,
                                                               strat_out=stratified,
                                                               proc=threads,
                                                               output_normfile=True)

            unstrat_pred.index.name = "function"
            unstrat_pred.reset_index(inplace=True)

            if not no_descrip and custom_trait_tables is None:
                unstrat_pred = add_descrip_col(inputfile=unstrat_pred,
                                               mapfile=default_map[func],
                                               in_df=True)

            # Write out stratified table only if that option was specified.
            if stratified:
                strat_pred.reset_index(inplace=True)

                if not no_descrip and custom_trait_tables is None:
                    strat_pred = add_descrip_col(inputfile=strat_pred,
                                                 mapfile=default_map[func],
                                                 in_df=True)

            # Strat pred will be None when --stratified unset.
            func_dfs[func] = [unstrat_pred, strat_pred]

    # Infer pathway abundances and coverages unless --no_pathways set.
    if not no_pathways:

        with TemporaryDirectory() as temp_dir:

            predicted_rxn = path.join(temp_dir, "pred_metagenome_rxn.tsv")

            if stratified:
                func_dfs[rxn_func][1].to_csv(path_or_buf=predicted_rxn, sep="\t", index=False)
            else:
                func_dfs[rxn_func][0].to_csv(path_or_buf=predicted_rxn, sep="\t", index=False)

            if verbose:
                print("Inferring pathways from predicted " + rxn_func)

            # Set regrouping mapfile to be empty if no_regroup set.
            if no_regroup:
                regroup_map = None

            unstrat_abun, unstrat_cov, strat_abun, strat_cov = run_minpath_pipeline(
                                                                        inputfile=predicted_rxn,
                                                                        mapfile=pathway_map,
                                                                        regroup_mapfile=regroup_map,
                                                                        proc=threads,
                                                                        out_dir=temp_dir,
                                                                        gap_fill=gap_fill_opt,
                                                                        per_sequence_contrib=per_sequence_contrib,
                                                                        print_cmds=verbose)

        unstrat_abun.reset_index(inplace=True)

        if not no_descrip and custom_trait_tables is None:
            unstrat_abun = add_descrip_col(inputfile=unstrat_abun,
                                           mapfile=default_map["METACYC"],
                                           in_df=True)

        unstrat_cov.reset_index(inplace=True)

        if not no_descrip and custom_trait_tables is None:
            unstrat_cov = add_descrip_col(inputfile=unstrat_cov,
                                        mapfile=default_map["METACYC"],
                                        in_df=True)

        # Write stratified output only if something besides None was returned.
        if strat_abun is not None:

            if not no_descrip and custom_trait_tables is None:
                strat_abun = add_descrip_col(inputfile=strat_abun,
                                             mapfile=default_map["METACYC"],
                                             in_df=True)

        if strat_cov is not None:

            if not no_descrip and custom_trait_tables is None:
                strat_cov = add_descrip_col(inputfile=strat_cov,
                                            mapfile=default_map["METACYC"],
                                            in_df=True)
        return(func_dfs, unstrat_abun, unstrat_cov, strat_abun, strat_cov)

    else:

        return(func_dfs, None, None, None, None)


if __name__ == "__main__":
    main()
