#!/bin/bash

OUT="lou"

mkdir -p $OUT
rm -rf $OUT/*
./leave_one_out.py ../../standard-RAxML/raxmlHPC ../bin/epa-ng data/ref.tre data/range_combined.fasta $OUT $1
