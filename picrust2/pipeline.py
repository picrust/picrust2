#!/usr/bin/env python

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.0.1-b"

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
                  hmm,
                  in_traits,
                  custom_trait_tables,
                  marker_gene_table,
                  pathway_map,
                  no_pathways,
                  regroup_map,
                  no_regroup,
                  stratified,
                  alignment_tool,
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

    if path.exists(output_folder):
        sys.exit("Stopping - output directory " + output_folder +
                 " already exists.")

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

        no_descrip = True

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

    if verbose:
        print("Placing sequences onto reference tree.")


    # Define folders for intermediate files.
    intermediate_dir = path.join(output_folder, "intermediate")
    place_seqs_intermediate = path.join(intermediate_dir, "place_seqs")
    make_output_dir(intermediate_dir)
    make_output_dir(place_seqs_intermediate)

    place_seqs_pipeline(study_fasta=study_fasta,
                        ref_msa=ref_msa,
                        tree=tree,
                        hmm=hmm,
                        out_tree=out_tree,
                        alignment_tool=alignment_tool,
                        threads=threads,
                        out_dir=place_seqs_intermediate,
                        chunk_size=5000,
                        print_cmds=verbose)

    if verbose:
        print("Finished placing sequences on output tree: " + out_tree)

    # Get predictions for all specified functions and keep track of outfiles.
    predicted_funcs = {}

    for func in funcs:

        count_outfile = hsp_pipeline_steps(func=func,
                                           calculate_NSTI=calculate_NSTI,
                                           out_tree=out_tree,
                                           func_table_in=func_tables[func],
                                           hsp_method=hsp_method,
                                           ci_setting=ci_setting,
                                           threads=threads,
                                           seed=seed,
                                           output_folder=output_folder,
                                           verbose=verbose)

        # Keep track of output file name for next step of pipeline.
        predicted_funcs[func] = count_outfile

    marker_infile = predicted_funcs["marker"]

    # Inititalize dictionary of function names to output filenames to return.
    func_output = {}

    # Each value will be a list of 2 elements corresponding to the unstratified
    # and stratified tables respectively (stratified will be None of not calculated).

    # Loop over each function again and run metagenome pipeline.
    for func in funcs:

        if func == "marker":
            continue

        if verbose:
            print("Running metagenome pipeline for " + func)
        
        func_infile = predicted_funcs[func]

        func_output_dir = path.join(output_folder, func + "_metagenome_out")

        func_map = None

        if func in default_map:
            func_map = default_map[func]

        func_strat_out, func_unstrat_out = metagenome_pipeline_steps(input_table=input_table,
                                                                     func_infile=func_infile,
                                                                     marker_infile=marker_infile,
                                                                     func_output_dir=func_output_dir,
                                                                     no_descrip=no_descrip,
                                                                     max_nsti=max_nsti,
                                                                     min_reads=min_reads,
                                                                     min_samples=min_samples,
                                                                     stratified=stratified,
                                                                     threads=threads,
                                                                     func_map=func_map,
                                                                     verbose=verbose)
        if stratified:
            func_output[func] = func_strat_out
        else:
            func_output[func] = func_unstrat_out


    pathway_outfiles = None

    # Infer pathway abundances and coverages unless --no_pathways set.
    if not no_pathways:

        pathways_intermediate = path.join(intermediate_dir, "pathways")
        make_output_dir(pathways_intermediate)

        if verbose:
            print("Inferring pathways from predicted " + rxn_func)

        predicted_rxn = func_output[rxn_func]

        # Set regrouping mapfile to be empty if no_regroup set.
        if no_regroup:
            regroup_map = None

        unstrat_abun, unstrat_cov, strat_abun, strat_cov = run_minpath_pipeline(
                                                                    inputfile=predicted_rxn,
                                                                    mapfile=pathway_map,
                                                                    regroup_mapfile=regroup_map,
                                                                    proc=threads,
                                                                    out_dir=pathways_intermediate,
                                                                    gap_fill=gap_fill_opt,
                                                                    per_sequence_contrib=per_sequence_contrib,
                                                                    print_cmds=verbose)

        pathways_out = path.join(output_folder, "pathways_out")

        unstrat_abun.index.name = 'pathway'
        unstrat_cov.index.name = 'pathway'
        unstrat_abun.reset_index(inplace=True)
        unstrat_cov.reset_index(inplace=True)

        pathway_outfiles = {}

        if not no_descrip:
            unstrat_abun = add_descrip_col(inputfile=unstrat_abun,
                                           mapfile=default_map["METACYC"],
                                           in_df=True)
        if not no_descrip:
            unstrat_cov = add_descrip_col(inputfile=unstrat_cov,
                                          mapfile=default_map["METACYC"],
                                          in_df=True)

        if verbose:
            print("Writing predicted pathway abundances and coverages to " + pathways_out)

        make_output_dir(pathways_out)

        unstrat_abun_outfile = path.join(pathways_out, "path_abun_unstrat.tsv")
        unstrat_cov_outfile = path.join(pathways_out, "path_cov_unstrat.tsv")

        unstrat_abun.to_csv(path_or_buf=unstrat_abun_outfile,  sep="\t", index=False)
        unstrat_cov.to_csv(path_or_buf=unstrat_cov_outfile,  sep="\t", index=False)

        pathway_outfiles["unstrat_abun"] = unstrat_abun_outfile
        pathway_outfiles["unstrat_cov"] = unstrat_cov_outfile

        strat_abun_outfile = None
        strat_cov_outfile = None

        # Write stratified output only if something besides None was returned.
        if strat_abun is not None:

            if not no_descrip:
                strat_abun = add_descrip_col(inputfile=strat_abun,
                                             mapfile=default_map["METACYC"],
                                             in_df=True)
            strat_abun_outfile = path.join(pathways_out, "path_abun_strat.tsv")
            strat_abun.to_csv(path_or_buf=strat_abun_outfile,  sep="\t", index=False)

        if strat_cov is not None:

            if not no_descrip:
                strat_cov = add_descrip_col(inputfile=strat_cov,
                                            mapfile=default_map["METACYC"],
                                            in_df=True)

            strat_cov_outfile = path.join(pathways_out, "path_cov_strat.tsv")
            strat_cov.to_csv(path_or_buf=strat_cov_outfile,  sep="\t", index=False)

        pathway_outfiles["strat_abun"] = strat_abun_outfile
        pathway_outfiles["strat_cov"] = strat_cov_outfile

    return(func_output, pathway_outfiles)


def hsp_pipeline_steps(func, calculate_NSTI, out_tree, func_table_in,
                       hsp_method, ci_setting, threads, seed, output_folder,
                       verbose):
    '''HSP pipeline steps moved to separate function for improved garbage
    collection (i.e. so that large objects no longer needed are removed from
    memory).'''

    # Only output NSTI in 16S table.
    nsti_setting = False
    if func == "marker" and calculate_NSTI:
        nsti_setting = True

    if verbose:
        print("Running hidden-state prediction for " + func)

    hsp_table, ci_table = castor_hsp_workflow(tree_path=out_tree,
                                              trait_table_path=func_table_in,
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

    return(count_outfile)


def metagenome_pipeline_steps(input_table, func_infile, marker_infile,
                              func_output_dir, no_descrip, max_nsti, min_reads,
                              min_samples, stratified, threads, func_map,
                              verbose):
    '''Steps wraping metagenome pipeline moved to separate function to decrease
    memory usage.'''

    # Infer metagenome abundances per-sample.
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

    if not no_descrip and func_map:
        unstrat_pred = add_descrip_col(inputfile=unstrat_pred,
                                       mapfile=func_map,
                                       in_df=True)

    # Write out stratified table only if that option was specified.
    if stratified:
        strat_pred.reset_index(inplace=True)

        if not no_descrip and func_map:
            strat_pred = add_descrip_col(inputfile=strat_pred,
                                         mapfile=func_map,
                                         in_df=True)

    if verbose:
        print("Writing metagenome output files for " + func + " to: " +
              func_output_dir)

    unstrat_outfile = path.join(func_output_dir, "pred_metagenome_unstrat.tsv")
    unstrat_pred.to_csv(path_or_buf=unstrat_outfile, sep="\t", index=False)    
    
    strat_outfile = None
    if stratified:
        strat_outfile = path.join(func_output_dir, "pred_metagenome_strat.tsv")
        strat_pred.to_csv(path_or_buf=strat_outfile, sep="\t", index=False)

    # Return output filenames.
    return(strat_outfile, unstrat_outfile)


if __name__ == "__main__":
    main()
