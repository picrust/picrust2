#!/usr/bin/env python

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.1.0-b"

from os import path
import sys
import biom
from picrust2.default import default_tables
from picrust2.place_seqs import identify_ref_files
from picrust2.util import (make_output_dir, check_files_exist, read_fasta,
                           system_call_check)


def full_pipeline(study_fasta,
                  input_table,
                  output_folder,
                  threads,
                  ref_dir,
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
                  skip_nsti,
                  skip_minpath,
                  no_gap_fill,
                  coverage,
                  per_sequence_contrib,
                  verbose):
    '''Function that contains wrapper commands for full PICRUSt2 pipeline.
    Descriptions of all of these input arguments/options are given in the
    picrust2_pipeline.py script.'''

    if path.exists(output_folder):
        sys.exit("Stopping since output directory " + output_folder +
                 " already exists.")

    # Make output folder.
    make_output_dir(output_folder)

    # Throw warning if --per_sequence_contrib set but --stratified unset.
    if per_sequence_contrib and not stratified:
        print("\nThe option --per_sequence_contrib was set, but not the option "
              "--stratified. This means that a stratified pathway table will "
              "be output only (i.e. a stratified metagenome table will NOT "
              "be output).\n", file=sys.stderr)

    out_tree = path.join(output_folder, "out.tre")

    if custom_trait_tables is None:

        # Check that specified functional categories are allowed.
        FUNC_TRAIT_OPTIONS = ['COG', 'EC', 'KO', 'PFAM', 'TIGRFAM']
        funcs = in_traits.split(",")
        for func in funcs:
            if func not in FUNC_TRAIT_OPTIONS:
                sys.exit("Error - specified category " + func + " is not " +
                         "one of the default categories.")

        # Add EC to this set if pathways are to be predicted.
        if "EC" not in funcs and not no_pathways:
            funcs.append("EC")

        rxn_func = "EC"

        func_tables = default_tables

    else:

        # Split paths to input custom trait tables and take the basename to be
        # the function id. The first table specified is assumed to be used
        # for inferring pathways.
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

    # Check that all input files exist. 
    ref_msa, tree, hmm, model = identify_ref_files(ref_dir)
    files2check = [study_fasta, input_table, ref_msa, tree, hmm, model] + list(func_tables.values())

    if not no_pathways:
        files2check.append(pathway_map)

        if not no_regroup:
            files2check.append(regroup_map)

    # This will throw an error if any input files are not found.
    check_files_exist(files2check)

    # Check that sequence names in FASTA overlap with input table.
    check_overlapping_seqs(study_fasta, input_table)

    if verbose:
        print("Placing sequences onto reference tree", file=sys.stderr)

    # Define folders for intermediate files.
    intermediate_dir = path.join(output_folder, "intermediate")
    make_output_dir(intermediate_dir)

    place_seqs_intermediate = path.join(intermediate_dir, "place_seqs")

    # Run place_seqs.py.
    place_seqs_cmd = ["place_seqs.py",
                      "--study_fasta", study_fasta,
                      "--ref_dir", ref_dir,
                      "--out_tree", out_tree,
                      "--threads", str(threads),
                      "--intermediate", place_seqs_intermediate,
                      "--chunk_size", str(5000)]

    if verbose:
        place_seqs_cmd.append("--print_cmds")

    system_call_check(place_seqs_cmd, print_out=verbose)

    if verbose:
        print("Finished placing sequences on output tree: " + out_tree,
              file=sys.stderr)

    # Get predictions for all specified functions and keep track of outfiles.
    predicted_funcs = {}

    for func in funcs:

        # Change output filename for NSTI and non-NSTI containing files.
        hsp_outfile = path.join(output_folder, func + "_predicted")

        if func == "marker" and not skip_nsti:
            hsp_outfile = hsp_outfile + "_and_nsti.tsv"
        else:
            hsp_outfile = hsp_outfile + ".tsv"

        # Keep track of output filename for next step of pipeline.
        predicted_funcs[func] = hsp_outfile

        # Run hsp.py for each function database.
        hsp_cmd = ["hsp.py",
                   "--tree", out_tree,
                   "--output", hsp_outfile,
                   "--observed_trait_table", func_tables[func],
                   "--hsp_method", hsp_method,
                   "--processes", str(threads),
                   "--seed", "100"]

        # Add flags to command if specified.
        if func == "marker" and not skip_nsti:
            hsp_cmd.append("--calculate_NSTI")

        system_call_check(hsp_cmd, print_out=verbose)

    # Now run metagenome pipeline commands.
    # Inititalize dictionary of function names --> metagenome output files.
    func_output = {}

    # Loop over each function again and run metagenome pipeline.
    for func in funcs:

        if func == "marker":
            continue

        if verbose:
            print("Running metagenome pipeline for " + func, file=sys.stderr)

        func_output_dir = path.join(output_folder, func + "_metagenome_out")

        metagenome_pipeline_cmd = ["metagenome_pipeline.py",
                                   "--input", input_table,
                                   "--function", predicted_funcs[func],
                                   "--marker", predicted_funcs["marker"],
                                   "--min_reads", str(min_reads),
                                   "--min_samples", str(min_samples),
                                   "--out_dir", func_output_dir]

        # Initialize 2-element list as value for each function.
        # First value will be unstratified output and second will be
        # stratified output.
        func_output[func] = [None, None]

        func_output[func][0] = path.join(func_output_dir,
                                         "pred_metagenome_unstrat.tsv")

        if not skip_nsti:
            metagenome_pipeline_cmd += ["--max_nsti", str(max_nsti)]

        if stratified:
            metagenome_pipeline_cmd.append("--strat_out")
            func_output[func][1] = path.join(func_output_dir,
                                             "pred_metagenome_strat.tsv")

        # Note that STDERR is printed for this command since it outputs how
        # many ASVs were above the NSTI cut-off (if specified).
        system_call_check(metagenome_pipeline_cmd, print_out=verbose,
                          print_stderr=True)

    # Now infer pathway abundances and coverages unless --no_pathways set.
    pathway_outfiles = None

    if not no_pathways:

        pathways_intermediate = path.join(intermediate_dir, "pathways")
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
                                "--proc", str(threads)]

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

        if per_sequence_contrib:
            pathway_pipeline_cmd.append("--per_sequence_contrib")

            norm_sequence_abun = path.join(output_folder,
                                           rxn_func + "_metagenome_out",
                                           "seqtab_norm.tsv")

            pathway_pipeline_cmd += ["--per_sequence_abun", norm_sequence_abun]

            pathway_pipeline_cmd += ["--per_sequence_function",
                                      predicted_funcs[rxn_func]]

        if verbose:
            pathway_pipeline_cmd.append("--print_cmds")

        system_call_check(pathway_pipeline_cmd, print_out=verbose)

        if verbose:
            print("Wrote predicted pathway abundances and coverages to " +
                  path_output_dir, file=sys.stderr)

        # Keep track of output filenames if this function is being used in
        # a non-default way (e.g. with a QIIME2 plugin).
        pathway_outfiles = {}

        pathway_outfiles["unstrat_abun"] = path.join(path_output_dir,
                                                     "path_abun_unstrat.tsv")
        pathway_outfiles["unstrat_cov"] = path.join(path_output_dir,
                                                    "path_cov_unstrat.tsv")

        if stratified:
            pathway_outfiles["strat_abun"] = path.join(path_output_dir,
                                                       "path_abun_strat.tsv")
            pathway_outfiles["strat_cov"] = path.join(path_output_dir,
                                                      "path_cov_strat.tsv")
        else:
            pathway_outfiles["strat_abun"] = None
            pathway_outfiles["strat_cov"] = None

    return(func_output, pathway_outfiles)


def check_overlapping_seqs(in_seq, in_tab):
    '''Check that ASV ids overlap between the input FASTA and sequence
    abundance table. Will throw an error if none overlap and will otherwise
    print number of overlapping ids to STDERR.'''

    FASTA_ASVs = set(read_fasta(in_seq).keys())

    BIOM_ASVs = set(biom.load_table(in_tab).ids(axis='observation'))

    num_ASV_overlap = len(BIOM_ASVs.intersection(FASTA_ASVs))

    # Throw error if 0 ASVs overlap between the two files.
    if num_ASV_overlap == 0:
        sys.exit("Stopping - no ASV ids overlap between input FASTA and "
                 "sequence abundance table")

    # Otherwise print to STDER how many ASVs overlap between the two files.
    print(str(num_ASV_overlap) + " of " + str(len(BIOM_ASVs)) + " sequence "
          "ids overlap between input table and FASTA.\n", file=sys.stderr)
