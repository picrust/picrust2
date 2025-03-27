import os
import argparse

parser = argparse.ArgumentParser(

            description="",

                epilog='''
                Example for running:
                python count_hmm_hits.py -hf hmmsearch_archaea -g alkB -o archaea_alkB.txt -e -25
                In this example, we'd count the number of hits within each file of the folder hmmsearch_archaea that have e-value below
                10 to the power of -25. These would be saved in tab-delimited format in the archaea_alkB.txt file, with the column header alkB
                
                ''', formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-hf', '--hit_folder', required=True,
                            type=str, help='output folder containing hmmsearch results for all genomes')
parser.add_argument('-g', '--gene_name', required=True,
                            type=str, help='Name of gene required for output')
parser.add_argument('-o', '--out_file', required=True,
                            type=str, help='Name of output file, required for writing counts of hits to')
parser.add_argument('-e', '--e_value', default='-25',
                            type=str, help='Threshold to include hits below this limit (default: %(default)s)')

args = parser.parse_args()
hit_folder, e_value, gene_name, out_file = args.hit_folder, args.e_value, args.gene_name, args.out_file

e_value = 10**(int(e_value))

files = os.listdir(hit_folder)
all_hits = ['assembly\t'+gene_name+'\n']
for f in files:
    hits = 0
    for row in open(hit_folder+'/'+f, 'r'):
        if row[0] == '#': continue
        row = list(filter(None, row.replace('\n', '').split(' ')))
        if float(row[4]) <= e_value:
            hits += 1
    all_hits.append(f.replace('.out', '')+'\t'+str(hits)+'\n')

all_hits[-1] = all_hits[-1].replace('\n', '')

with open(out_file, 'w') as f:
    for h in all_hits:
        w = f.write(h)
