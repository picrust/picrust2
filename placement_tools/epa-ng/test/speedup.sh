#!/bin/bash

ABSPATH=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd)
OUT=$ABSPATH/speedup
LOG=$OUT/log
TREE=$ABSPATH/data/lucas/20k.newick
MSA=$ABSPATH/data/lucas/1k.fasta
REF_MSA=$ABSPATH/data/lucas/1k_reference.fasta
QRY_MSA=$ABSPATH/data/lucas/1k_query_100.fasta
mkdir -p $OUT
rm -rf $OUT/*
touch $LOG

export OMP_PROC_BIND=true

for i in 20 24 28 32 36 40 44 48
do
	export OMP_NUM_THREADS=$i
	echo "$i THREADS:" >> $LOG
	for j in {1..5}
	do
		{ time ../bin/epamk -t $TREE -s $REF_MSA -q $QRY_MSA  -w $OUT > /dev/null 2>&1 ; } 2>> $LOG
	done
done
