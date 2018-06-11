#!/usr/bin/python

import sys
import os
from dendropy import *
from subprocess import call
import glob

def help():
    print "USAGE:\traxml_path epa_path tree_file reference_MSA_file query_MSA_file output_dir [number of runs]"
def wrng(msg):
    print msg
    help()
    exit()
def err(msg):
    print msg
    print "Aborting"
    exit()

if len(sys.argv) < 7 or len(sys.argv) > 7:
    print "incorrect number of arguments"
    help()
    exit()

first_x = False
runs=0

raxml = sys.argv[1]
epa = sys.argv[2]
tree_file = os.path.abspath(sys.argv[3])
ref_MSA_file = os.path.abspath(sys.argv[4])
query_MSA_file = os.path.abspath(sys.argv[5])
output_dir = os.path.abspath(sys.argv[6])
if len(sys.argv) == 8:
    runs = int(sys.argv[7])
    first_x = True

if not os.path.isfile(raxml):
    wrng("raxml doesn't exist or isn't a file")
if not os.path.isfile(epa):
    wrng("epa doesn't exist or isn't a file")
if not os.path.isfile(tree_file):
    wrng("tree_file doesn't exist or isn't a file")
if not os.path.isfile(ref_MSA_file):
    wrng("reference_MSA_file doesn't exist or isn't a file")
if not os.path.isfile(query_MSA_file):
    wrng("query_MSA_file doesn't exist or isn't a file")
if not os.path.isdir(output_dir):
    wrng("output_dir doesn't exist or isn't a directory")

# read in tree and MSA files
tree = Tree.get(path=tree_file, schema="newick", rooting="force-unrooted")
msa = DnaCharacterMatrix.get(path=ref_MSA_file, schema="fasta")

num_failed = 0
num_run = 0

# count the tips
num_tips = 0.0
for n in tree.leaf_node_iter():
    num_tips += 1
progress = 0.0

with open(os.path.join(output_dir, "results.log"), 'wb') as log_file:
    # for every tip:
    for node in tree.leaf_node_iter():
        # trim taxon from the cloned tree
        lou_tree = tree.clone()
        to_prune = TaxonNamespace()
        to_prune.add_taxon(node.taxon)
        lou_tree.prune_taxa(to_prune)
        lou_tree.deroot()

        # write tree to tmp folder
        cur_outdir = os.path.join(output_dir, node.taxon.label)
        # print cur_outdir
        if not os.path.exists(cur_outdir):
            os.makedirs(cur_outdir)

        cur_treefile = os.path.join(cur_outdir, "tree.newick")
        lou_tree.write(path=cur_treefile, schema="newick", suppress_rooting=True)

        # call raxml with trimmed files
        # print "calling raxml:"
        params = [raxml, "-f", "v", "-s", ref_MSA_file, "-t", cur_treefile, "-n", "leave_one_out", "-m","GTRGAMMA",
        "-w", cur_outdir]
        # print params
        ret = call(params, stdout=open(os.devnull, 'wb'))

        # call epa with trimmed files
	params = [epa, "-t", cur_treefile, "-s", ref_MSA_file,"-q", ref_MSA_file , "-w", cur_outdir]        

	# print "calling epa to dump binary: "
        #params = [epa, "-t", cur_treefile,"-s", ref_MSA_file, "-OB", "-w", cur_outdir]
        #ret = call(params, stdout=open(os.devnull, 'wb'))

        # print "calling epa from binary file: "
        #binfile = os.path.join(cur_outdir, "epa_binary_file")
        #params = [epa, "-b", binfile, "-q", ref_MSA_file, "-w", cur_outdir]
        # print params
        ret = call(params, stdout=open(os.devnull, 'wb'))


        # call validation script on both jplace files, log to log file
        jplace_files = glob.glob(os.path.join(cur_outdir, "*.jplace"))
        assert len(jplace_files) == 2

        raxml_jplace =  os.path.join(cur_outdir, "RAxML_portableTree.leave_one_out.jplace")
        epa_jplace = os.path.join(cur_outdir, "epa_result.jplace")

        # print "calling compare script:"
        params = ["./pquery_emd", raxml_jplace, epa_jplace]
        # # print params
        ret = call(params, stdout=log_file)

        # print "calling compare script:"
       	#params = ["./pquery_emd", raxml_jplace, epa_jplace]
        # print params
       	#log_file.flush()
        #ret = call(params, stdout=log_file)

        log_file.write("\n")

        num_failed += ret
        num_run += 1

        if first_x and num_run == runs:
            break

        progress_old = progress
        progress = (num_run / num_tips)*100

        if  (progress - progress_old) > 1:
            print str(progress) + "%"

    failed_string = "Failed " + str(num_failed) + " out of " + str(num_run) + " tests"
    log_file.write(failed_string + "\n")

print failed_string
