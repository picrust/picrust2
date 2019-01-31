#!/usr/bin/env python

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.0.4-b"

import sys
from os import path, chdir, getcwd
from picrust2.util import (system_call_check, make_output_dir, read_fasta,
                           read_phylip, write_fasta, write_phylip,
                           read_stockholm)


def place_seqs_pipeline(study_fasta,
                        ref_dir,
                        out_tree,
                        threads,
                        out_dir,
                        chunk_size,
                        print_cmds):
    '''Full pipeline for running sequence placement.'''

    # Identify reference files to use.
    ref_msa, tree, hmm, model = identify_ref_files(ref_dir)

    # Run hmmalign to place study sequences into reference MSA.
    out_stockholm = path.join(out_dir, "query_align.stockholm")

    system_call_check("hmmalign  --trim --dna --mapali " + ref_msa + " --informat FASTA -o " +
                      out_stockholm + " " + hmm + " " + study_fasta,
                      print_out=print_cmds)

    hmmalign_out = read_stockholm(out_stockholm, clean_char=True)

    # Specify split FASTA files to be created.
    study_msa_fastafile = path.join(out_dir, "study_seqs_hmmalign.fasta")
    ref_msa_fastafile = path.join(out_dir, "ref_seqs_hmmalign.fasta")

    ref_seqnames = set(list(read_fasta(ref_msa).keys()))

    study_seqnames = set(read_fasta(study_fasta).keys())

    ref_hmmalign_subset = {seq: hmmalign_out[seq] for seq in ref_seqnames}
    study_hmmalign_subset = {seq: hmmalign_out[seq] for seq in study_seqnames}

    write_fasta(ref_hmmalign_subset, ref_msa_fastafile)
    write_fasta(study_hmmalign_subset, study_msa_fastafile)

    # Run EPA-NG to output .jplace file.
    epa_out_dir = path.join(out_dir, "epa_out")

    run_epa_ng(tree=tree, ref_msa_fastafile=ref_msa_fastafile,
               study_msa_fastafile=study_msa_fastafile, model=model,
               chunk_size=chunk_size, threads=threads, out_dir=epa_out_dir,
               print_cmds=print_cmds)

    jplace_outfile = path.join(epa_out_dir, "epa_result.jplace")

    gappa_jplace_to_newick(jplace_file=jplace_outfile, outfile=out_tree,
                           print_cmds=print_cmds)


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
                      " -n out", print_out=print_cmds)

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
    '''Run EPA-NG on specified tree, reference MSA, and study sequence MSA.
    Will output a .jplace file in out_dir.'''

    make_output_dir(out_dir)

    epa_ng_command = ["epa-ng", "--tree", tree, "--ref-msa", ref_msa_fastafile,
                      "--query", study_msa_fastafile, "--chunk-size",
                      str(chunk_size), "-T", str(threads), "-m", model, "-w",
                      out_dir]

    system_call_check(epa_ng_command, print_out=print_cmds)

def gappa_jplace_to_newick(jplace_file: str, outfile: str, print_cmds=False):
    '''System call to gappa binary to convert jplace object to newick
    treefile (with specified filename).'''

    gappa_out_dir = path.dirname(jplace_file)

    # Run gappa to convert jplace to newick.
    system_call_check("gappa analyze graft --jplace-path " + jplace_file +
                      " --fully-resolve --out-dir " + gappa_out_dir,
                      print_out=print_cmds)

    # Expected name of output newick file.
    newick_file = jplace_file.replace(".jplace", ".newick")

    # Rename newick file to be specified outfile.
    system_call_check("mv " + newick_file + " " + outfile,
                      print_out=print_cmds)

def identify_ref_files(in_dir):
    '''Given a directory will check whether the four required reference files
    are present and will return the path to each file in a list in the order:
    FASTA, TREE, HMM, MODEL.'''

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

    other_expected = [path.join(in_dir, base_path + ".tre"),
                      path.join(in_dir, base_path + ".hmm"),
                      path.join(in_dir, base_path + ".model")]

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
