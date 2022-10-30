#!/usr/bin/env python

__copyright__ = "Copyright 2018-2022, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.5.1"

import sys
from os import path, chdir, getcwd
import json
import numpy as np
from picrust2.util import (system_call_check, make_output_dir, read_fasta,
                           read_phylip, write_fasta, write_phylip,
                           read_stockholm)


def place_seqs_pipeline(study_fasta,
                        ref_dir,
                        placement_tool,
                        out_tree,
                        threads,
                        out_dir,
                        min_align,
                        chunk_size,
                        verbose):
    '''Full pipeline for running sequence placement.'''

    # Throw error if there is a space in the study FASTA filepath.
    if " " in study_fasta:
        sys.exit("Stopping - remove the space from the input FASTA filepath.")

    # Check header-lines of FASTA: if they have multiple whitespace-delimited fields
    # this can sometimes be an issue for hmmalign.
    check_fasta_headers(study_fasta)


    # Identify reference files to use. Note that model file will be different
    # depending on which placement pipeline is indicated.
    ref_msa, tree, hmm, model = identify_ref_files(ref_dir, placement_tool)

    # Run hmmalign to align study sequences with reference MSA.
    # This is necessary for EPA-ng and is used as a check for sequences that
    # should be excluded for both the EPA-ng and SEPP pipelines (based on
    # sequences that align very little and are likely problematic).
    out_stockholm = path.join(out_dir, "query_align.stockholm")

    system_call_check("hmmalign --trim --dna --mapali " +
                      ref_msa + " --informat FASTA -o " +
                      out_stockholm + " " + hmm + " " + study_fasta,
                      print_command=verbose, print_stdout=verbose,
                      print_stderr=verbose)

    hmmalign_out = read_stockholm(out_stockholm, clean_char=True)

    ref_seqnames = set(list(read_fasta(ref_msa).keys()))
    study_seqs = read_fasta(study_fasta)

    ref_hmmalign_subset = {seq: hmmalign_out[seq] for seq in ref_seqnames}

    study_hmmalign_subset = check_alignments(raw_seqs=study_seqs,
                                             aligned_seqs=hmmalign_out,
                                             min_align=min_align,
                                             verbose=verbose)
    if placement_tool == "epa-ng":

        # Write out subsetted alignment and reference subset.
        study_msa_fastafile = path.join(out_dir, "study_seqs_hmmalign.fasta")
        ref_msa_fastafile = path.join(out_dir, "ref_seqs_hmmalign.fasta")

        write_fasta(ref_hmmalign_subset, ref_msa_fastafile)
        write_fasta(study_hmmalign_subset, study_msa_fastafile)

        # Run EPA-ng to place input sequences and output JPLACE file.
        epa_out_dir = path.join(out_dir, "epa_out")

        run_epa_ng(tree=tree,
                   ref_msa_fastafile=ref_msa_fastafile,
                   study_msa_fastafile=study_msa_fastafile,
                   model=model,
                   chunk_size=chunk_size,
                   threads=threads,
                   out_dir=epa_out_dir,
                   print_cmds=verbose)

        jplace_outfile = path.join(epa_out_dir, "epa_result_parsed.jplace")

    elif placement_tool == "sepp":

        # For SEPP, filter out any funky sequences from original FASTA before
        # running.
        study_seqs_subset = {seq: study_seqs[seq] for seq in study_hmmalign_subset}
        study_fasta_filt = path.join(out_dir, "study_seqs_filtered.fasta")
        write_fasta(study_seqs_subset, study_fasta_filt)

        sepp_out_dir = path.join(out_dir, "sepp_out")

        run_sepp(tree=tree,
                 ref_msa_fastafile=ref_msa,
                 study_msa_fastafile=study_fasta_filt,
                 raxml_model=model,
                 threads=threads,
                 out_dir=sepp_out_dir,
                 print_cmds=verbose)

        jplace_outfile = path.join(sepp_out_dir, "output_placement.json")

    else:
        sys.exit("Option placement_tool needs to be either \"epa-ng\" or \
                 \"sepp\". It was set to this instead: \"" + placement_tool +
                 "\".")

    gappa_jplace_to_newick(jplace_file=jplace_outfile,
                           outfile=out_tree,
                           print_cmds=verbose)


def run_papara(tree: str, ref_msa: dict, study_fasta: str, out_dir: str,
               threads=1, print_cmds=False):
    '''Run PaPaRa to place study sequences into reference multiple-sequence
    alignment (MSA). Will return dictionary of the the output MSA (sequence ids
    as keys). Expects path to tree and study FASTA as strings. Expects
    reference MSA as a dictionary output by read_fasta. This MSA will be
    converted to phylip format before running PaPaRa.'''

    # Get absolute paths to input files.
    tree = path.abspath(tree)
    study_fasta = path.abspath(study_fasta)

    # Change working directory to out directory (but keep track of original).
    # This is necessary because PaPaRa outputs into the current working
    # directory.
    orig_wd = getcwd()
    chdir(out_dir)

    # Convert ref sequences from MSA FASTA to phylip.
    write_phylip(ref_msa, "ref_seqs.phylip")

    # Make call to papara to place sequences (outputs phylip format).
    system_call_check("papara -t " + tree + " -s ref_seqs.phylip " +
                      "-q " + study_fasta + " -j " + str(threads) +
                      " -n out", print_command=print_cmds,
                      print_stdout=print_cmds, print_stderr=print_cmds)

    # Change back to original working directory.
    chdir(orig_wd)

    # Read in papara phylip output and return.
    return(read_phylip(path.join(out_dir, "papara_alignment.out"),
                       check_input=True))

def split_ref_study_papara(papara_out: dict, ref_seqnames: set, ref_fasta: str,
                           study_fasta: str):
    '''Split PaPaRa phylip output into FASTA MSA files of study sequences and
    reference sequences separately. Expects PaPaRa output already read in
    as dictionary. Takes in the PaPaRa output as a dictionary, a set that
    contains all sequence ids in reference MSA, and the output FASTA
    filenames.'''

    # Determine study sequence id based on those found in the all and ref sets.
    all_seqnames = set(list(papara_out.keys()))
    study_seqnames = all_seqnames.difference(ref_seqnames)

    # Get subsets of PaPaRa output MSA of reference study sequences only.
    ref_papara_subset = {seq: papara_out[seq] for seq in ref_seqnames}
    study_papara_subset = {seq: papara_out[seq] for seq in study_seqnames}

    write_fasta(ref_papara_subset, ref_fasta)
    write_fasta(study_papara_subset, study_fasta)

def run_epa_ng(tree: str, ref_msa_fastafile: str, study_msa_fastafile: str,
               model: str, out_dir: str, chunk_size=5000,
               threads=1, print_cmds=False):
    '''Run EPA-ng on specified tree, reference MSA, and study sequence MSA.
    Will output a .jplace file in out_dir.'''

    make_output_dir(out_dir)

    epa_ng_command = ["epa-ng", "--tree", tree,
                      "--ref-msa", ref_msa_fastafile,
                      "--query", study_msa_fastafile,
                      "--chunk-size", str(chunk_size),
                      "-T", str(threads),
                      "-m", model,
                      "-w", out_dir,
                      "--filter-acc-lwr", "0.99",
                      "--filter-max", "100"]

    system_call_check(epa_ng_command, print_command=print_cmds,
                      print_stdout=print_cmds, print_stderr=print_cmds)

    # Parse jplace file so that output is reprodicible.
    jplace_orig = path.join(out_dir, "epa_result.jplace")
    jplace_parsed = path.join(out_dir, "epa_result_parsed.jplace")
    parse_jplace(jplace_orig, jplace_parsed)

def gappa_jplace_to_newick(jplace_file: str, outfile: str, print_cmds=False):
    '''System call to gappa binary to convert jplace object to newick
    treefile (with specified filename).'''

    gappa_out_dir = path.dirname(jplace_file)

    # Run gappa to convert jplace to newick.
    system_call_check("gappa examine graft --jplace-path " + jplace_file +
                      " --fully-resolve --out-dir " + gappa_out_dir,
                      print_command=print_cmds, print_stdout=print_cmds,
                      print_stderr=print_cmds)

    # Expected name of output newick file.
    jplace_ext = path.splitext(jplace_file)[1]
    newick_file = jplace_file.replace(jplace_ext, ".newick")

    # Rename newick file to be specified outfile.
    system_call_check("mv " + newick_file + " " + outfile,
                      print_command=print_cmds)

def identify_ref_files(in_dir, placement_method):
    '''Given a directory will check whether the four required reference files
    are present and will return the path to each file in a list in the order:
    FASTA, TREE, HMM, MODEL. Will return paths to different model files
    depending if the EPA-ng or SEPP placement methods are specified.'''

    # Remove any trailing slashes.
    in_dir = in_dir.rstrip('/')

    base_path = path.basename(in_dir)

    # Expected filenames to be in this directory.
    # Note that only one of the possible FASTA filenames should be in the
    # directory. An error will be thrown if multiple files are found.
    possible_fasta = [path.join(in_dir, base_path + ".fna.gz"),
                      path.join(in_dir, base_path + ".fasta.gz"),
                      path.join(in_dir, base_path + ".fna"),
                      path.join(in_dir, base_path + ".fasta")]

    path2return = []

    # Flag for whether FASTA file is present.
    missing_fasta = False

    # Identify which possible FASTA files are present.
    for poss in possible_fasta:
        if path.isfile(poss):
            path2return.append(poss)

    # Throw error if multiple FASTAs identified.
    if len(path2return) == 0:
        missing_fasta = True
    elif len(path2return) > 1:
        sys.exit("Found multiple FASTA files in specified reference directory. "
                 "Only one should be present. Files found: \n" +
                 "\n".join(path2return))

    # Check whether other reference files exist as well.
    # List to keep track of which required files are missing.
    missing_files = []

    if placement_method == "epa-ng":
        expected_modelfile = path.join(in_dir, base_path + ".model")
    elif placement_method == "sepp":
        expected_modelfile = path.join(in_dir, base_path + ".raxml_info")
    else:
        sys.exit("Option placement_method needs to be either \"epa-ng\" or \
                 \"sepp\". It was set to this instead: \"" + placement_method +
                 "\".")

    other_expected = [path.join(in_dir, base_path + ".tre"),
                      path.join(in_dir, base_path + ".hmm"),
                      expected_modelfile]

    for other in other_expected:
        if path.isfile(other):
            path2return.append(other)
        else:
            missing_files.append(other)
    
    if missing_fasta:
        print("No FASTA file found in specified directory. Expected to find "
              "one of:\n" + "\n".join(possible_fasta) + "\n\n", file=sys.stderr)

        if len(missing_files) > 0:
            print("In addition to the missing FASTA file, the following "
                  "expected file(s) could not be found:\n" +
                  "\n".join(missing_files) + "\n\n", file=sys.stderr)
    elif len(missing_files) > 0:
            if len(missing_files) > 0:
                print("The following expected file(s) could not be found:\n" +
                      "\n".join(missing_files) + "\n\n", file=sys.stderr)

    if missing_fasta or len(missing_files) > 0:
        sys.exit("Error - missing at least one of the four reference files in "
                 "this specified directory: " + in_dir)

    return(path2return)

def parse_jplace(jplace_in, jplace_out):
    '''Parse jplace file to retain only a single placement per ASV (with the
    highest likelihood and highest edge number). Also, will order that the
    ASVs are reported in the file alphanumerically. This is needed so that
    GAPPA constructs a newick treefile consistently. Requires input and output
    jplace filenames to be specified.'''

    # Read in jplace file as JSON object.
    with open(jplace_in, 'r') as f:
        datastore = json.load(f)

    # Loop through all ASV placements and parse out ASV names and info.
    placement_names = []

    for i in range(len(datastore["placements"])):

        placement_names.append(datastore["placements"][i]["n"][0])

        asv_placements = datastore["placements"][i]["p"]
        
        # For ASVs with multiple placements only retain a single placement.
        if len(asv_placements) == 1:
            next
        else:
            placement_info = np.array(asv_placements)

            # First identify all indices matching the max likelihood.
            max_lik = np.amax(placement_info[:, 1])

            max_lik_i = np.where(placement_info[:, 1] == max_lik)

            # Then figure out the highest edge number of all placements with
            # this likelihood and retain only that placement.
            max_edge_num = np.amax(placement_info[max_lik_i, 0])

            placement2keep = np.where(placement_info[:, 0] == max_edge_num)[0]

            if len(placement2keep) != 1:
                sys.exit("Multiple ASVs on same edge in jplace file.")

            datastore["placements"][i]["p"] = [datastore["placements"][i]["p"][placement2keep[0]]]

    # Sort placements by ASV names.
    placement_names_sorted_i = np.argsort(placement_names)
    datastore["placements"] = [datastore["placements"][i] for i in placement_names_sorted_i]

    # Write out sorted and parsed jplace object.
    with open(jplace_out, 'w') as f:
        json.dump(datastore, f, indent=4, sort_keys=False)

def check_alignments(raw_seqs, aligned_seqs, min_align, verbose):
    '''Check that all study sequences are at least a high % of their original
    length in alignment. If not then remove sequences below this cut-off and
    throw detailed warning. Throw critical error if all sequences are below
    cut-off. This cut-off is specified by --min_align option. Also keep
    track of input sequence lengths and report range of lengths (only if 
    verbose set).'''

    poorly_aligned = []
    passing = {}

    min_study_seq_length = None
    max_study_seq_length = None

    raw_seqnames = set(raw_seqs.keys())

    for seq_id in raw_seqnames:
        orig_seq = raw_seqs[seq_id]
        aligned_seq = aligned_seqs[seq_id]

        orig_seq = orig_seq.replace("-", "")
        orig_seq = orig_seq.replace(".", "")

        aligned_seq = aligned_seq.replace("-", "")
        aligned_seq = aligned_seq.replace(".", "")

        orig_seq_length = len(orig_seq)

        if len(aligned_seq) < orig_seq_length * min_align:
            poorly_aligned.append(seq_id)
        else:
            passing[seq_id] = aligned_seqs[seq_id]

        if verbose:
            if not min_study_seq_length or orig_seq_length < min_study_seq_length:
                min_study_seq_length = orig_seq_length

            if not max_study_seq_length or orig_seq_length > max_study_seq_length:
                max_study_seq_length = orig_seq_length

    if len(poorly_aligned) == len(raw_seqnames):
        sys.exit("Stopping - all " + str(len(raw_seqnames)) + " input "
                 "sequences aligned poorly to reference sequences (--min_align "
                 "option specified a minimum proportion of " + str(min_align) +
                 " aligning to reference sequences).")

    elif len(poorly_aligned) > 0:
        print("Warning - " + str(len(poorly_aligned)) + " input sequences "
              "aligned poorly to reference sequences (--min_align option "
              "specified a minimum proportion of " + str(min_align) +
              " aligning to reference sequences). These input sequences will "
              "not be placed and will be excluded from downstream steps.\n\n"
              "This is the set of poorly aligned input sequences to be "
              "excluded: " + ", ".join(poorly_aligned) + "\n", file=sys.stderr)

    if verbose:
        if min_study_seq_length == max_study_seq_length:
            print("All raw input sequences were the same length (" +
                  str(min_study_seq_length) + ")\n", file=sys.stderr)
        else:
            print("Raw input sequences ranged in length from " +
                  str(min_study_seq_length) + " to " +
                  str(max_study_seq_length) + "\n", file=sys.stderr)

    return(passing)


def run_sepp(tree: str, ref_msa_fastafile: str, study_msa_fastafile: str,
             raxml_model: str, out_dir: str, threads=1, print_cmds=False,
             set_seed=297834):
    '''Run SEPP on specified tree, reference MSA, and study sequence MSA.
    Will output a .jplace file in out_dir.'''

    sepp_command = ["run_sepp.py",
                    "--tree", tree,
                    "--raxml", raxml_model,
                    "--cpu", str(threads),
                    "--molecule", "dna",
                    "--outdir", out_dir,
                    "-seed", str(set_seed),
                    "--alignment", ref_msa_fastafile,
                    "--fragment", study_msa_fastafile]

    # Give error if either reference or query FASTA gzipped.
    ref_fasta_basename = path.basename(ref_msa_fastafile)
    ref_fasta_ext = path.splitext(ref_fasta_basename)[1]
    
    study_fasta_basename = path.basename(study_msa_fastafile)
    study_fasta_ext = path.splitext(study_fasta_basename)[1]

    if ref_fasta_ext == ".gz" or study_fasta_ext == ".gz":

        if ref_fasta_ext == ".gz":

            print("\nTo place sequences with SEPP all input FASTAs must be "
                  "decompressed. Please run gunzip on this reference file "
                  "before re-running: " + ref_msa_fastafile + "\n",
                  file=sys.stderr)
        
        if study_fasta_ext == ".gz":

            print("\nTo place sequences with SEPP all input FASTAs must be "
                  "decompressed. Please run gunzip on the query FASTA file "
                  "before re-running: " + study_msa_fastafile + "\n",
                  file=sys.stderr)
        
        sys.exit("\nStopped running due to at least one input FASTA being "
                 "gzipped (which SEPP does not allow).\n")

    system_call_check(sepp_command, print_command=print_cmds,
                      print_stdout=print_cmds, print_stderr=print_cmds)


def check_fasta_headers(filename : str):
    '''Check that headers of specified FASTA do not contain whitespace.'''

    if filename[-3:] == ".gz":
        fasta_in = gzip.open(filename, "rt")
    else:
        fasta_in = open(filename, "r")

    for line in fasta_in:

        line = line.rstrip()

        if len(line) == 0:
            continue

        if line[0] == ">":

            if len(line.split()) > 1:
                sys.exit("\nStopping - input FASTA file header lines should not contain whitespace. Please alter and re-run. \n\nHeader with whitespace:\n\n" + line + "\n")

