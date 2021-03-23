#!/usr/bin/env python

__copyright__ = "Copyright 2018-2021, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.4.0-b"

from os import path
import sys
from picrust2.default import (default_tables, default_pathway_map,
                              default_ref_dir)
from picrust2.place_seqs import identify_ref_files
from picrust2.util import (make_output_dir, check_files_exist, read_fasta,
                           system_call_check, read_seqabun)


def full_pipeline(study_fasta,
                  input_table,
                  output_folder,
                  processes,
                  placement_tool,
                  ref_dir,
                  in_traits,
                  custom_trait_tables,
                  marker_gene_table,
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
                  skip_nsti,
                  skip_minpath,
                  no_gap_fill,
                  coverage,
                  per_sequence_contrib,
                  wide_table,
                  skip_norm,
                  remove_intermediate,
                  verbose):
    '''Function that contains wrapper commands for full PICRUSt2 pipeline.
    Descriptions of all of these input arguments/options are given in the
    picrust2_pipeline.py script.'''

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

        func_tables = default_tables

    else:
        # Split paths to input custom trait tables and take the basename to be
        # the function id.
        funcs = []
        func_tables = {}

        for custom in custom_trait_tables.split(","):

            func_id = path.splitext(path.basename(custom))[0]
            funcs.append(func_id)
            func_tables[func_id] = custom

    # Add reaction function to be in set of gene families if it is not already
    # and as long as pathways are also to be predicted.
    if rxn_func not in funcs and not no_pathways:
        orig_rxn_func = rxn_func
        rxn_func = path.splitext(path.basename(rxn_func))[0]
        funcs.append(rxn_func)

        if rxn_func not in func_tables:
            func_tables[rxn_func] = orig_rxn_func

    if not skip_norm:
        # Append marker as well, since this also needs to be run.
        funcs.append("marker")
        func_tables["marker"] = marker_gene_table

    # Check that all input files exist.
    ref_msa, tree, hmm, model = identify_ref_files(ref_dir, placement_tool)
    files2check = [study_fasta, input_table, ref_msa, tree, hmm, model] + list(func_tables.values())

    if not no_pathways:
        files2check.append(pathway_map)

        # Throw warning if default pathway mapfile used with non-default
        # reference files.
        if pathway_map == default_pathway_map and ref_dir != default_ref_dir:
            print("Warning - non-default reference files specified with "
                  "default pathway mapfile of prokaryote-specific MetaCyc "
                  "pathways (--pathway_map option). This usage may be "
                  "unintended.", file=sys.stderr)

        if not no_regroup:
            files2check.append(regroup_map)

    # This will throw an error if any input files are not found.
    check_files_exist(files2check)

    # Check that sequence names in FASTA overlap with input table.
    check_overlapping_seqs(study_fasta, input_table, verbose)

    if path.exists(output_folder):
        sys.exit("Stopping since output directory " + output_folder +
                 " already exists.")

    # Make output folder.
    make_output_dir(output_folder)

    if verbose:
        print("Placing sequences onto reference tree", file=sys.stderr)

    # Define folders for intermediate files (unless --remove_intermediate set).
    if remove_intermediate:
        place_seqs_intermediate = ""
        pathways_intermediate = ""
    else:
        intermediate_dir = path.join(output_folder, "intermediate")
        make_output_dir(intermediate_dir)
        place_seqs_intermediate = path.join(intermediate_dir, "place_seqs")
        pathways_intermediate = path.join(intermediate_dir, "pathways")

    # Run place_seqs.py.
    place_seqs_cmd = ["place_seqs.py",
                      "--study_fasta", study_fasta,
                      "--ref_dir", ref_dir,
                      "--out_tree", out_tree,
                      "--processes", str(processes),
                      "--intermediate", place_seqs_intermediate,
                      "--min_align", str(min_align),
                      "--chunk_size", str(5000),
                      "--placement_tool", placement_tool]

    if verbose:
        place_seqs_cmd.append("--verbose")

    system_call_check(place_seqs_cmd, print_command=verbose,
                      print_stdout=verbose, print_stderr=True)

    if verbose:
        print("Finished placing sequences on output tree: " + out_tree,
              file=sys.stderr)

    # Get predictions for all specified functions and keep track of outfiles.
    predicted_funcs = {}

    if not skip_norm:
        # Make sure marker database is first in the list. This is because this will
        # be run on a single core and so will be easier to identify any errors
        # if the program exits when working on this function type.
        funcs.insert(0, funcs.pop(funcs.index("marker")))

    for func in funcs:
        # Change output filename for NSTI and non-NSTI containing files.
        hsp_outfile = path.join(output_folder, func + "_predicted")

        if (func == "marker" and not skip_nsti) or (skip_norm and not skip_nsti):
            hsp_outfile = hsp_outfile + "_and_nsti.tsv.gz"
        else:
            hsp_outfile = hsp_outfile + ".tsv.gz"

        # Keep track of output filename for next step of pipeline.
        predicted_funcs[func] = hsp_outfile

        # Run hsp.py for each function database.
        hsp_cmd = ["hsp.py",
                   "--tree", out_tree,
                   "--output", hsp_outfile,
                   "--observed_trait_table", func_tables[func],
                   "--hsp_method", hsp_method,
                   "--edge_exponent", str(edge_exponent),
                   "--seed", "100"]

        # Add flags to command if specified.
        if (func == "marker" and not skip_nsti) or (skip_norm and not skip_nsti):
            hsp_cmd.append("--calculate_NSTI")

        # Run marker on only 1 processor.
        if func == "marker":
            hsp_cmd += ["--processes", "1"]
        else:
            hsp_cmd += ["--processes", str(processes)]

        if verbose:
            hsp_cmd.append("--verbose")

        system_call_check(hsp_cmd, print_command=verbose,
                          print_stdout=verbose, print_stderr=True)

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

        if not skip_nsti:
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
