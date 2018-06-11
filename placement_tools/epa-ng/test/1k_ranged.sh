#!/bin/bash

ABSPATH=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd)
OUT=$ABSPATH/1k_ranged
LOG=$OUT/log
TREE=$ABSPATH/data/lucas/20k.newick
MSA=$ABSPATH/data/lucas/1k.fasta
mkdir -p $OUT
rm -rf $OUT/*
touch $LOG

echo "RUNNING RAXML" >> $LOG
(time ../../standard-RAxML/raxmlHPC-SSE3 -f v -s $MSA -t $TREE -n 1k -m GTRGAMMA -w $OUT) &>> $LOG

echo "RUNNING EPA" > $LOG
(time ./epa_ranged -t $TREE -s $MSA -Or -w $OUT) &>> $LOG

echo "RUNNING JPLACE_COMPARE" >> $LOG
(./jplace_compare.py -v $OUT/RAxML_portableTree.1k.jplace $OUT/epa_result.jplace) &>> $LOG
