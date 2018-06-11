#!/usr/bin/python
import sys
import os
from dendropy import *

tree = Tree.get(path='data/lucas/20k.newick', schema="newick")
msa = DnaCharacterMatrix.get(path='data/lucas/1k.fasta', schema="fasta")

exception_list = ["928f713140c528eca38f36a03c4ba0cbfba31ca7", "6d40160742d6426f2e6160745a816d556206ddb1",
"23526fb1ece3ff6533674ea7c6012bcf68b50e01"]

remove_list = []
lbl = ""
num = 900

for s in msa:
    lbl = s.label.replace('_', ' ')
    if (tree.find_node_with_taxon_label(lbl) == None and num > 0):# and lbl not in exception_list):
        remove_list += [s]
        num -= 1

print len(remove_list)

msa.discard_sequences(remove_list)

with open("data/lucas/problem_set.fasta", "wb") as outfile:
    msa.write(file=outfile, schema='fasta')
