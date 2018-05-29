#!/usr/bin/env python

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.0.0-b.2"

from os import path, chdir, getcwd
from picrust2.util import (system_call_check, make_output_dir, read_fasta,
                           read_phylip, write_fasta, write_phylip)


def place_seqs_pipeline(study_fasta,
                        ref_msa,
                        tree,
                        out_tree,
                        threads,
                        papara_output,
                        out_dir,
                        chunk_size,
                        print_cmds):
    '''Full pipeline for running sequence placement.'''

    # Read in ref seqs FASTA as a dict.
    ref_msa = read_fasta(ref_msa)

    # Either read in PaPaRa output or run it.
    if papara_output:
        # Read in PaPaRa output if already done.
        papara_out = read_phylip(papara_output, check_input=True)

    else:
        # Run PaPaRa to place study sequences and read in Phylip file.
        papara_out = run_papara(tree=tree, ref_msa=ref_msa,
                                study_fasta=study_fasta, out_dir=out_dir,
                                threads=threads, print_cmds=print_cmds)

    # Specify split FASTA files to be created.
    study_msa_fastafile = path.join(out_dir, "study_seqs_papara.fasta")
    ref_msa_fastafile = path.join(out_dir, "ref_seqs_papara.fasta")

    # Split PaPaRa output into two FASTA files containging study and reference
    # sequences respectively.
    split_ref_study_papara(papara_out=papara_out,
                           ref_seqnames=set(list(ref_msa.keys())),
                           study_fasta=study_msa_fastafile,
                           ref_fasta=ref_msa_fastafile)

    # Run EPA-NG to output .jplace file.
    epa_out_dir = path.join(out_dir, "epa_out")

    run_epa_ng(tree=tree, ref_msa_fastafile=ref_msa_fastafile,
               study_msa_fastafile=study_msa_fastafile, chunk_size=chunk_size,
               threads=threads, out_dir=epa_out_dir, print_cmds=print_cmds)

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
               out_dir: str, chunk_size=5000, threads=1, print_cmds=False):
    '''Run EPA-NG on specified tree, reference MSA, and study sequence MSA.
    Will opath.joinutput a .jplace file in out_dir.'''

    make_output_dir(out_dir)

    system_call_check("epa-ng --tree " + tree + " --ref-msa " +
                      ref_msa_fastafile + " --query " + study_msa_fastafile +
                      " --chunk-size " + str(chunk_size) + " -T " +
                      str(threads) + " -w " + out_dir, print_out=print_cmds)


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
