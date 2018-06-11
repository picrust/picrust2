#!/bin/bash
#$ -S /bin/bash
#$ -pe mvapich16 96
#$ -binding linear:16
#$ -q bridge.q
#$ -cwd
#$ -j y 
#$ -l h_rt=00:30:00

#@ job_type = MPICH
#@ class = micro
#@ node = 20
#@ tasks_per_node = 16
#@ island_count = 1
#@ energy_policy_tag = epa
#@ minimize_time_to_solution = yes
#@ wall_clock_limit = 04:00:00
#@ job_name = pepa_TRACE_VT
#@ network.MPI = sn_all,not_shared,us
#@ initialdir = /home/barberpe/epa/test
#@ output = out_$(jobid).log
#@ error  = err_$(jobid).log
#@ notification=always
#@ notify_user=pierre.barbera@h-its.org
#@ queue


source /etc/profile.d/modules.sh
module load openmpi/gcc
module unload gcc
module load gcc

#export VT_MAX_FLUSHES=0
#export VT_PROFM_LDIR=$(ws_allocate vt_trace_tmp 5)
export OMP_PROC_BIND=true
export OMP_NUM_THREADS=16

ABSPATH=/home/barberpe/epa/test
OUT=$ABSPATH/trace
LOG=$OUT/log
TREE=$ABSPATH/data/lucas/20k.newick
MSA=$ABSPATH/data/lucas/1k.fasta
REF_MSA=$ABSPATH/data/lucas/1k_reference.fasta
QRY_MSA=$ABSPATH/data/lucas/1k_query.fasta
mkdir -p $OUT
rm -rf $OUT/*
touch $LOG

#mpiexec ../bin/epa-ng -t $TREE -s $REF_MSA -q $QRY_MSA -B -w $OUT
mpiexec ../bin/epa-ng -b $ABSPATH/binfile -q $QRY_MSA -w $OUT #-g 0.99

