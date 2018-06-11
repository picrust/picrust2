#!/bin/bash

ABSPATH=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd)
OUT=$ABSPATH/1k
LOG=$OUT/log
TREE=$ABSPATH/data/lucas/20k.newick
MSA=$ABSPATH/data/lucas/1k.fasta
REF_MSA=$ABSPATH/data/lucas/1k_reference.fasta
QRY_MSA=$ABSPATH/data/lucas/1k_query_100.fasta
mkdir -p $OUT
rm -rf $OUT/*
touch $LOG

#echo "RUNNING RAXML" >> $LOG
#(time ../../standard-RAxML/raxmlHPC-SSE3 -f v -H -s $MSA -t $TREE -n 1k -m GTRGAMMA -w $OUT) &>> $LOG

#echo "RUNNING EPA - DUMPING BINARY" >> $LOG
#(time ../bin/epa-ng -t $TREE -s $REF_MSA -OB -w $OUT) &>> $LOG

#echo "RUNNING EPA - FROM BINARY" >> $LOG
#(time ../bin/epa-ng -b $OUT/epa_binary_file -q $QRY_MSA -O -w $OUT) &>> $LOG

echo "RUNNING EPA" >> $LOG
export OMP_NUM_THREADS=4
valgrind --tool=helgrind ../bin/epa-ng -t $TREE -s $REF_MSA -q $QRY_MSA -g 0.95 -w $OUT

#echo "RUNNING PPLACER" >> $LOG
#(time ./pplacer -t data/lucas/20k.newick -s 1k/RAxML_info.1k --out-dir $OUT data/lucas/1k.fasta) &>> $LOG

#echo "RUNNING JPLACE_COMPARE" >> $LOG
#(./jplace_compare.py -v $OUT/RAxML_portableTree.1k.jplace $OUT/epa_result.jplace) &>> $LOG
