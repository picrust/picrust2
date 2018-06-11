#!/bin/bash

ABSPATH=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd)
OUT=$ABSPATH/runtime_thorough
LOG=$OUT/log
TREE=$ABSPATH/data/lucas/20k.newick
MSA=$ABSPATH/data/lucas/1k_100.fasta
REF_MSA=$ABSPATH/data/lucas/1k_reference.fasta
QRY_MSA=$ABSPATH/data/lucas/1k_query_100.fasta
mkdir -p $OUT
rm -rf $OUT/*
touch $LOG

echo "EPA" >> $LOG
for i in 1
do
	{ time ../bin/epamk -t $TREE -s $REF_MSA -q $QRY_MSA  -w $OUT > /dev/null 2>&1 ; } 2>> $LOG
done


echo "RAXML" >> $LOG
for i in 1
do
	rm -rf $OUT/RAxML*
        { time ../../standard-RAxML/raxmlHPC-AVX -f v -H -s $MSA -t $TREE -n 1k -m GTRGAMMA -w $OUT > /dev/null 2>&1 ; } 2>> $LOG
done

echo "PPLACER" >> $LOG
for i in 1
do
        { time ./pplacer -t $TREE -s 1k/RAxML_info.1k --no-pre-mask --max-strikes 0 -j 0 --out-dir $OUT $MSA > /dev/null 2>&1 ; } 2>> $LOG
done

