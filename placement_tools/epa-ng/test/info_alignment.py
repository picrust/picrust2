#!/usr/bin/python

from Bio import SeqIO
from os import path
import sys
import numpy as np

def help():
    print "Script to deterimine some characteristics of a MSA file"
    print "USAGE:\t<path to fasta file>"
def wrng(msg):
    print msg
    help()
    exit()
def err(msg):
    print msg
    print "Aborting"
    exit()

def get_valid_span(seq):
    upper = len(seq)
    lower = 0
    while seq[lower] == '-':
        lower += 1
    while seq[upper - 1] == '-':
        upper -= 1
    return upper - lower

if len(sys.argv) < 2 or len(sys.argv) > 2:
    print "incorrect number of arguments"
    help()
    exit()

fasta_file = path.abspath(sys.argv[1])
if not path.isfile(fasta_file):
    wrng("File doesn't exist or isn't a file")

valid_fractions = []
count = 0
for record in SeqIO.parse(fasta_file, "fasta"):
    seq = str(record.seq).upper()
    valid_fractions.append(get_valid_span(seq)/float(len(seq)))
    count +=1

print "Number of sequences parsed: " + str(count)
print "Fraction of sequence that is delimited by gap-only regions"
print "\tMax: " + str(np.max(valid_fractions))
print "\tMin: " + str(np.min(valid_fractions))
print "\tMean: " + str(np.mean(valid_fractions))
