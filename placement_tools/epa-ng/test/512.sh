#!/bin/bash

OUT="lou"
mkdir -p $OUT
rm -rf $OUT/*
./leave_one_out.py ../../standard-RAxML/raxmlHPC ../bin/epa data/lucas/tree.newick data/lucas/msa.fasta $OUT $1
