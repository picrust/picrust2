#!/bin/bash

OUT="lou/rerun"
EPA="../bin/epa-ng"
LUCAS="data/lucas"
BASEDIR=$(dirname $0)

mkdir -p "$OUT/$1"

echo "RUNNING EPA"
$EPA -t "$1/tree.newick" -s "$LUCAS/msa.fasta" -O -w "$OUT/$1"

echo "OLD RESULTS:"
./jplace_compare.py -v "$1/RAxML_portableTree.leave_one_out.jplace" "$1/epa_result.jplace"

echo "NEW RESULTS:"
./jplace_compare.py -v "$1/RAxML_portableTree.leave_one_out.jplace" "$OUT/$1/epa_result.jplace"
