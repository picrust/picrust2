#!/usr/bin/env python

from os import path
import sys
from picrust2.default import (default_ref_dir_bac, default_ref_dir_arc, default_tables_bac,
                              default_tables_arc, default_pathway_map)
from picrust2.place_seqs import identify_ref_files
from picrust2.util import (make_output_dir, check_files_exist, read_fasta,
                           read_fasta_ids, system_call_check, read_seqabun, prune_tree)
from picrust2.split_domains import get_lowest_nsti, combine_domain_predictions


def full_pipeline_split(study_fasta,
                  input_table,
                  output_folder,
                  processes,
                  placement_tool,
                  ref_dir1,
                  ref_dir2,
                  in_traits,
                  custom_trait_tables_ref1,
                  custom_trait_tables_ref2,
                  marker_gene_table_ref1,
                  marker_gene_table_ref2,
                  pathway_map,
                  rxn_func,
                  no_pathways,
                  regroup_map,
                  no_regroup,
                  stratified,
                  max_nsti,
                  min_reads,
                  min_samples,
                  hsp_method,
                  edge_exponent,
                  min_align,
                  skip_minpath,
                  no_gap_fill,
                  coverage,
                  per_sequence_contrib,
                  wide_table,
                  skip_norm,
                  remove_intermediate,
                  verbose):
    '''Function that contains wrapper commands for full PICRUSt2 pipeline.
    This is the version that incorporates separate placement steps for two 
    domains (by default, bacteria and archaea).
    Descriptions of all of these input arguments/options are given in the
    picrust2_pipeline_split.py script.'''

    # Throw warning if --per_sequence_contrib set but --stratified unset.
    if per_sequence_contrib and not stratified:
        print("\nThe option --per_sequence_contrib was set, but not the option "
              "--stratified. This means that a stratified pathway table will "
              "be output only (i.e. a stratified metagenome table will NOT "
              "be output).\n", file=sys.stderr)
    
    if ref_dir1 == default_ref_dir_bac:
        out_tree_ref1 = path.join(output_folder, "bac.tre")
        name_ref1 = 'bac'
    else:
        out_tree_ref1 = path.join(output_folder, "ref1.tre")
        name_ref1 = 'ref1'
    if ref_dir2 == default_ref_dir_arc:
        out_tree_ref2 = path.join(output_folder, "arc.tre")
        name_ref2 = 'arc'
    else:
        out_tree_ref2 = path.join(output_folder, "ref2.tre")
        name_ref2 = 'ref2'
    
    # Exit if only one set of custom trait tables has been gives
    if custom_trait_tables_ref1 is None and not custom_trait_tables_ref2 is None:
        sys.exit("You've set some custom trait tables for reference set 1 but not "
              "for reference set 2. Please set both of them.")
    elif not custom_trait_tables_ref1 is None and custom_trait_tables_ref2 is None:
        sys.exit("You've set some custom trait tables for reference set 2 but not "
              "for reference set 1. Please set both of them.")

    if custom_trait_tables_ref1 is None:

        # Check that specified functional categories are allowed.
        FUNC_TRAIT_OPTIONS = ['EC', 'KO', 'GO', 'PFAM', 'BIGG', 'CAZY', 'GENE_NAMES']
        funcs = in_traits.split(",")
        for func in funcs:
            if func not in FUNC_TRAIT_OPTIONS:
                sys.exit("Error - specified category " + func + " is not " +
                         "one of the default categories.")
                         
        funcs_ref1 = funcs
        funcs_ref2 = funcs

        func_tables_ref1 = default_tables_bac
        func_tables_ref2 = default_tables_arc

    else:
        # Split paths to input custom trait tables and take the basename to be
        # the function id.
        funcs_ref1, funcs_ref2 = [], []
        func_tables_ref1, func_tables_ref2 = {}, {}

        for custom in custom_trait_tables_ref1.split(","):

            func_id = path.splitext(path.basename(custom))[0]
            funcs_ref1.append(func_id)
            func_tables_ref1[func_id] = custom
            
        for custom in custom_trait_tables_ref2.split(","):

            func_id = path.splitext(path.basename(custom))[0]
            funcs_ref2.append(func_id)
            func_tables_ref2[func_id] = custom

    # Add reaction function to be in set of gene families if it is not already
    # and as long as pathways are also to be predicted.
    # Note that by default this is EC, so we would only be adding it if EC wasn't 
    # given as an option but we do want pathways to be predicted.
    if rxn_func not in funcs_ref1 and not no_pathways:
        orig_rxn_func = rxn_func
        rxn_func = path.splitext(path.basename(rxn_func))[0]
        funcs_ref1.append(rxn_func)
        funcs_ref2.append(rxn_func)

        if rxn_func not in func_tables_ref1:
            func_tables_ref1[rxn_func] = orig_rxn_func
        if rxn_func not in func_tables_ref2:
            func_tables_ref2[rxn_func] = orig_rxn_func

    # if not skip_norm:
    #     # Append marker as well, since this also needs to be run.
    #     funcs_ref1.append("marker")
    #     func_tables_ref1["marker"] = marker_gene_table
    # Append marker
    # funcs_ref1.append("marker")
    # funcs_ref2.append("marker")
    # func_tables_ref1["marker"] = marker_gene_table_ref1
    # func_tables_ref2["marker"] = marker_gene_table_ref2

    # Check that all input files exist.
    ref_msa_ref1, tree_ref1, hmm_ref1, model_ref1 = identify_ref_files(ref_dir1, placement_tool)
    ref_msa_ref2, tree_ref2, hmm_ref2, model_ref2 = identify_ref_files(ref_dir2, placement_tool)
    files2check = [study_fasta, input_table, ref_msa_ref1, tree_ref1, hmm_ref1, model_ref1, ref_msa_ref2, tree_ref2, hmm_ref2, model_ref2] + list(func_tables_ref1.values())  + list(func_tables_ref2.values())

    if not no_pathways:
        files2check.append(pathway_map)

        # Throw warning if default pathway mapfile used with non-default
        # reference files.
        if pathway_map == default_pathway_map and ref_dir1 != default_ref_dir_bac:
            print("Warning - non-default reference files specified with "
                  "default pathway mapfile of prokaryote-specific MetaCyc "
                  "pathways (--pathway_map option). This usage may be "
                  "unintended.", file=sys.stderr)

        if not no_regroup:
            files2check.append(regroup_map)

    # This will throw an error if any input files are not found.
    check_files_exist(files2check)

    # This will throw an error if any trait tables are entirely empty.
    check_files_exist(list(func_tables_ref1.values())  + list(func_tables_ref2.values()))

    # Check that sequence names in FASTA overlap with input table.
    check_overlapping_seqs(study_fasta, input_table, verbose)
    
    # Check that there are no duplicated sequence names in the FASTA.
    check_duplicated_seqnames(study_fasta, verbose)

    if path.exists(output_folder):
        sys.exit("Stopping since output directory " + output_folder +
                 " already exists.")

    # Make output folder.
    make_output_dir(output_folder)

    if verbose:
        print("Placing sequences onto reference tree", file=sys.stderr)

    # Define folders for intermediate files (unless --remove_intermediate set).
    if remove_intermediate:
        place_seqs_intermediate_ref1 = ""
        place_seqs_intermediate_ref2 = ""
        pathways_intermediate = ""
    else:
        intermediate_dir = path.join(output_folder, "intermediate")
        make_output_dir(intermediate_dir)
        place_seqs_intermediate_ref1 = path.join(intermediate_dir, "place_seqs_"+name_ref1)
        place_seqs_intermediate_ref2 = path.join(intermediate_dir, "place_seqs_"+name_ref2)
        pathways_intermediate = path.join(intermediate_dir, "pathways")

    # Run place_seqs.py.
    place_seqs_cmd_ref1 = ["place_seqs.py",
                      "--study_fasta", study_fasta,
                      "--ref_dir", ref_dir1,
                      "--out_tree", out_tree_ref1,
                      "--processes", str(processes),
                      "--intermediate", place_seqs_intermediate_ref1,
                      "--min_align", str(min_align),
                      "--chunk_size", str(5000),
                      "--placement_tool", placement_tool]
                      
    place_seqs_cmd_ref2 = ["place_seqs.py",
                      "--study_fasta", study_fasta,
                      "--ref_dir", ref_dir2,
                      "--out_tree", out_tree_ref2,
                      "--processes", str(processes),
                      "--intermediate", place_seqs_intermediate_ref2,
                      "--min_align", str(min_align),
                      "--chunk_size", str(5000),
                      "--placement_tool", placement_tool]
                      

    if verbose:
        place_seqs_cmd_ref1.append("--verbose")
        place_seqs_cmd_ref2.append("--verbose")

    system_call_check(place_seqs_cmd_ref1, print_command=verbose,
                      print_stdout=verbose, print_stderr=True)

    if verbose:
        print("Finished placing sequences on output tree for reference 1: " + out_tree_ref1,
              file=sys.stderr)
              
    system_call_check(place_seqs_cmd_ref2, print_command=verbose,
                      print_stdout=verbose, print_stderr=True)

    if verbose:
        print("Finished placing sequences on output tree for reference 2: " + out_tree_ref2,
              file=sys.stderr)
              
    # Get predictions for all specified functions and keep track of outfiles.
    predicted_funcs_split = {}
              
    # Get predictions for marker gene sets for each reference
    reference_sets = [[ref_dir1, out_tree_ref1, marker_gene_table_ref1, name_ref1], [ref_dir2, out_tree_ref2, marker_gene_table_ref2, name_ref2]]
    for ref in reference_sets:
        ref_dir = ref[0]
        out_tree = ref[1]
        marker_table = ref[2]
        name_ref = ref[3]
        
        hsp_outfile = path.join(output_folder, name_ref + "_marker_predicted_and_nsti.tsv.gz")
        
        predicted_funcs_split[name_ref + "_marker"] = hsp_outfile
        
        # Run hsp.py for both reference sets
        hsp_cmd = ["hsp.py",
                   "--tree", out_tree,
                   "--output", hsp_outfile,
                   "--observed_trait_table", marker_table,
                   "--hsp_method", hsp_method,
                   "--edge_exponent", str(edge_exponent),
                   "--seed", "100",
                   "--calculate_NSTI",
                   "--processes", "1"]
                   
        if verbose:
            hsp_cmd.append("--verbose")

        system_call_check(hsp_cmd, print_command=verbose,
                          print_stdout=verbose, print_stderr=True)
                          
    if verbose:
        print("Finished getting marker and NSTI predictions for both domains: " +  predicted_funcs_split[name_ref1 + "_marker"] + ", " + predicted_funcs_split[name_ref2 + "_marker"] +
              "\nNow finding the best one for each sequence.",
              file=sys.stderr)
    
    # Choose which of the reference sets is best for each sequence
    nsti_lowest_df, nsti_dom1_df, nsti_dom2_df = get_lowest_nsti(predicted_funcs_split[name_ref1 + "_marker"],
                                                                 predicted_funcs_split[name_ref2 + "_marker"],
                                                                 name_ref1,
                                                                 name_ref2,
                                                                 verbose=verbose)
    
    # Write the new files
    nsti_dom1_df.to_csv(path_or_buf=path.join(output_folder, name_ref1 + "_reduced_marker_predicted_and_nsti.tsv.gz"), sep="\t")
    nsti_dom2_df.to_csv(path_or_buf=path.join(output_folder, name_ref2 + "_reduced_marker_predicted_and_nsti.tsv.gz"), sep="\t")
    nsti_lowest_df.to_csv(path_or_buf=path.join(output_folder, "combined_marker_predicted_and_nsti.tsv.gz"), sep="\t")
    
    # Check whether we still have both domains present still
    if nsti_dom1_df.shape[0] == 0:
        if verbose:
            print("Don't have any " + name_ref1 + " in the study sequences. Continuing with only " + name_ref2,
                  file=sys.stderr)
        
        # Update the predicted_funcs_split dictionary
        predicted_funcs_split[name_ref1 + "_marker"] = ""
        predicted_funcs_split[name_ref2 + "_marker"] = path.join(output_folder, name_ref2 + "_reduced_marker_predicted_and_nsti.tsv.gz")
    elif nsti_dom2_df.shape[0] == 0:
        if verbose:
            print("Don't have any " + name_ref2 + " in the study sequences. Continuing with only " + name_ref1,
                  file=sys.stderr)
        
        # Update the predicted_funcs_split dictionary
        predicted_funcs_split[name_ref1 + "_marker"] = path.join(output_folder, name_ref1 + "_reduced_marker_predicted_and_nsti.tsv.gz")
        predicted_funcs_split[name_ref2 + "_marker"] = ""
    else:
        # Update the predicted_funcs_split dictionary
        predicted_funcs_split[name_ref1 + "_marker"] = path.join(output_folder, name_ref1 + "_reduced_marker_predicted_and_nsti.tsv.gz")
        predicted_funcs_split[name_ref2 + "_marker"] = path.join(output_folder, name_ref2 + "_reduced_marker_predicted_and_nsti.tsv.gz")
    
    # Now prune the trees to contain only the sequences of interest for each domain
    out_tree_ref1_red = path.join(output_folder, name_ref1 + "_reduced.tre")
    out_tree_ref2_red = path.join(output_folder, name_ref2 + "_reduced.tre")
    
    ref_d1_seqs = read_fasta_ids(ref_msa_ref1)
    ref_d2_seqs = read_fasta_ids(ref_msa_ref2)
    
    prune_tree(list(nsti_dom1_df.index.values)+ref_d1_seqs, out_tree_ref1, out_tree_ref1_red)
    prune_tree(list(nsti_dom2_df.index.values)+ref_d2_seqs, out_tree_ref2, out_tree_ref2_red)
    
    if verbose:
        print("Finished getting the best domain match for each sequence: " + path.join(output_folder, "combined_marker_predicted_and_nsti.tsv.gz") +
              " Now running hsp.py for the reduced reference sets.",
              file=sys.stderr)
    
    # Run hsp.py for each function database for each of the reference sets
    reference_sets = [[name_ref1, funcs_ref1, out_tree_ref1_red, func_tables_ref1, nsti_dom1_df], [name_ref2, funcs_ref2, out_tree_ref2_red, func_tables_ref2, nsti_dom2_df]]
    for ref in reference_sets:
        # First get the files we'll be using
        name = ref[0]
        funcs = ref[1]
        out_tree = ref[2]
        func_tables = ref[3]
        nsti = ref[4]
        if nsti.shape[0] == 0: 
            continue
        
        for func in funcs:
            hsp_outfile = path.join(output_folder, name + "_" + func + "_predicted.tsv.gz")
        
            predicted_funcs_split[name + "_" + func] = hsp_outfile
            
            hsp_cmd = ["hsp.py",
                       "--tree", out_tree,
                       "--output", hsp_outfile,
                       "--observed_trait_table", func_tables[func],
                       "--hsp_method", hsp_method,
                       "--edge_exponent", str(edge_exponent),
                       "--seed", "100",
                       "--processes", str(processes)]
                       
            if verbose:
                hsp_cmd.append("--verbose")

            system_call_check(hsp_cmd, print_command=verbose,
                          print_stdout=verbose, print_stderr=True)
                          
    if verbose:
        print("Finished getting functional predictions for all traits.",
              file=sys.stderr)
        
    # Get a list of the predictions to be joined together and join them for each domain
    # Note that this checks the names so that it is compatible with custom trait ables
    functions_to_combine = []
    for func in predicted_funcs_split:
        if "marker" not in func:
            functions_to_combine.append(func.replace(name_ref1+'_', "").replace(name_ref2+'_', ""))
    
    # Get only unique values
    functions_to_combine = set(functions_to_combine)
    
    # Join the predictions for each domain
    predicted_funcs = {}
    for func in functions_to_combine:
      
        combining = []
        
        for pred in predicted_funcs_split:
            if func in pred:
                  combining.append(pred)
        
        if len(combining) == 1:
            predicted_funcs[func] = predicted_funcs_split[combining[0]]
            print("Warning: There was only one file for the function: "+ func + "\n"
                  "Maybe that's fine if you used custom traits or there were no sequences "
                  "matching one of the domains.", file=sys.stderr)
        
        elif len(combining) > 2:
            sys.exit("More than two files were available for the function: " + func +"\n"
                     "If you are using your own custom trait files, please check that you "
                     "have not given the same table more than once for a domain.")
        
        else:
            print(combining, predicted_funcs_split)
            out_file = path.join(output_folder, "combined_" + func + "_predicted.tsv.gz")
            combine_domain_predictions(predicted_funcs_split[combining[0]], 
                                       predicted_funcs_split[combining[1]], 
                                       out_file,
                                       verbose=verbose)
            predicted_funcs[func] = out_file
    
    if verbose:
        print("Finished joining together all trait tables with the same trait name for both domains.",
              file=sys.stderr)
    
    predicted_funcs["marker"] = path.join(output_folder, "combined_marker_predicted_and_nsti.tsv.gz")

    # Now run metagenome pipeline commands.
    # Inititalize dictionary of function names --> metagenome output files.
    func_output = {}

    # Loop over each function again and run metagenome pipeline.
    for func in predicted_funcs:
      
        if func == "marker":
            continue

        if verbose:
            print("Running metagenome pipeline for " + func, file=sys.stderr)

        func_output_dir = path.join(output_folder, func + "_metagenome_out")

        metagenome_pipeline_cmd = ["metagenome_pipeline.py",
                                   "--input", input_table,
                                   "--function", predicted_funcs[func],
                                   "--min_reads", str(min_reads),
                                   "--min_samples", str(min_samples),
                                   "--out_dir", func_output_dir]

        # Initialize two-element list as value for each function.
        # First value will be unstratified output and second will be
        # stratified output.
        func_output[func] = [None, None]

        func_output[func][0] = path.join(func_output_dir,
                                         "pred_metagenome_unstrat.tsv.gz")

        if wide_table:
            metagenome_pipeline_cmd.append("--wide_table")

        metagenome_pipeline_cmd += ["--max_nsti", str(max_nsti)]

        if skip_norm:
            metagenome_pipeline_cmd.append("--skip_norm")
        else:
            metagenome_pipeline_cmd += ["--marker", predicted_funcs["marker"]]

        if stratified:
            metagenome_pipeline_cmd.append("--strat_out")

            if wide_table:
                func_output[func][1] = path.join(func_output_dir,
                                                 "pred_metagenome_strat.tsv.gz")
            else:
                func_output[func][1] = path.join(func_output_dir,
                                                 "pred_metagenome_contrib.tsv.gz")

        system_call_check(metagenome_pipeline_cmd, print_command=verbose,
                          print_stdout=verbose, print_stderr=True)

    # Now infer pathway abundances and coverages unless --no_pathways set.
    pathway_outfiles = None

    if not no_pathways:

        path_output_dir = path.join(output_folder, "pathways_out")

        if verbose:
            print("Inferring pathways from predicted " + rxn_func)

        # Determine whether stratified or unstratified table should be input.
        if not stratified or per_sequence_contrib:
            rxn_input_metagenome = func_output[rxn_func][0]
        else:
            rxn_input_metagenome = func_output[rxn_func][1]

        pathway_pipeline_cmd = ["pathway_pipeline.py",
                                "--input", rxn_input_metagenome,
                                "--out_dir", path_output_dir,
                                "--map", pathway_map,
                                "--intermediate", pathways_intermediate,
                                "--proc", str(processes)]

        if no_gap_fill:
            pathway_pipeline_cmd.append("--no_gap_fill")

        if skip_minpath:
            pathway_pipeline_cmd.append("--skip_minpath")

        if coverage:
            pathway_pipeline_cmd.append("--coverage")

        if no_regroup:
            pathway_pipeline_cmd.append("--no_regroup")
        else:
            pathway_pipeline_cmd += ["--regroup_map", regroup_map]

        if wide_table:
            pathway_pipeline_cmd.append("--wide_table")

        if per_sequence_contrib:
            pathway_pipeline_cmd.append("--per_sequence_contrib")

            if skip_norm:
                norm_sequence_abun = input_table
            else:
                norm_sequence_abun = path.join(output_folder,
                                               rxn_func + "_metagenome_out",
                                               "seqtab_norm.tsv.gz")

            pathway_pipeline_cmd += ["--per_sequence_abun", norm_sequence_abun]

            pathway_pipeline_cmd += ["--per_sequence_function",
                                     predicted_funcs[rxn_func]]

        if verbose:
            pathway_pipeline_cmd.append("--verbose")

        system_call_check(pathway_pipeline_cmd, print_command=verbose,
                          print_stdout=False, print_stderr=True)

        if verbose:
            print("Wrote predicted pathway abundances and coverages to " +
                  path_output_dir, file=sys.stderr)

        # Keep track of output filenames if this function is being used in
        # a non-default way (e.g. with a QIIME2 plugin).
        pathway_outfiles = {}

        pathway_outfiles["unstrat_abun"] = path.join(path_output_dir,
                                                     "path_abun_unstrat.tsv.gz")
        pathway_outfiles["unstrat_cov"] = path.join(path_output_dir,
                                                    "path_cov_unstrat.tsv.gz")

        pathway_outfiles["strat_abun"] = None
        pathway_outfiles["strat_cov"] = None

        if stratified or per_sequence_contrib:
            if wide_table:
                pathway_outfiles["strat_abun"] = path.join(path_output_dir,
                                                           "path_abun_strat.tsv.gz")

                if per_sequence_contrib:
                    pathway_outfiles["strat_cov"] = path.join(path_output_dir,
                                                              "path_cov_strat.tsv.gz")

            else:
                pathway_outfiles["strat_abun"] = path.join(path_output_dir,
                                                           "path_abun_contrib.tsv.gz")
                if per_sequence_contrib:
                    pathway_outfiles["strat_cov"] = path.join(path_output_dir,
                                                              "path_cov_contrib.tsv.gz")

    return(func_output, pathway_outfiles)


def check_overlapping_seqs(in_seq, in_tab, verbose):
    '''Check that ASV ids overlap between the input FASTA and sequence
    abundance table. Will throw an error if none overlap and will otherwise
    print number of overlapping ids to STDERR. Also throw warning if input
    ASV table contains a column called taxonomy'''

    FASTA_ASVs = set(read_fasta(in_seq).keys())

    in_table = read_seqabun(in_tab)

    table_ASVs = set(in_table.index.values)

    num_ASV_overlap = len(table_ASVs.intersection(FASTA_ASVs))

    if 'taxonomy' in in_table.columns:
        print("Warning - column named \"taxonomy\" in abundance table - if "
              "this corresponds to taxonomic labels this should be removed "
              "before running this pipeline.", file=sys.stderr)

    # Throw error if 0 ASVs overlap between the two files.
    if num_ASV_overlap == 0:
        sys.exit("Stopping - no ASV ids overlap between input FASTA and "
                 "sequence abundance table")

    # Otherwise print to STDERR how many ASVs overlap between the two files
    # if verbose set.
    if verbose:
        print(str(num_ASV_overlap) + " of " + str(len(table_ASVs)) +
              " sequence ids overlap between input table and FASTA.\n",
              file=sys.stderr)
              
def check_duplicated_seqnames(in_seq, verbose):
    '''Check that ASV ids are not duplicated in the input FASTA.
    Will throw an error if any of the ids are duplicated, and will otherwise
    print number of ids to STDERR.'''
    
    FASTA_ASVs = read_fasta_ids(in_seq)
    
    unique_ids = set(FASTA_ASVs)
    
    duplicated_ids = [i for i in unique_ids if FASTA_ASVs.count(i) > 1]
    
    if len(FASTA_ASVs) != len(set(FASTA_ASVs)):
      sys.exit("Stopping - there are duplicated ASV ids in the input FASTA.\n" +
              "These are the duplicated IDs: " + ", ".join(duplicated_ids))
    
    if verbose:
        print(str(len(FASTA_ASVs)) + " ASVs in the input FASTA. None of the sequence ids are duplicated.\n",
              file=sys.stderr)
