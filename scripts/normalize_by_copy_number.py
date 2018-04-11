#!/usr/bin/env python

from __future__ import division

__license__ = "GPL"
__version__ = "2-alpha.7"

import argparse
from biom import load_table, Table
from picrust2.predict_metagenomes import transfer_observation_metadata
from os import path
from os.path import join
from picrust2.util import get_picrust_project_dir, make_output_dir_for_file, write_biom_table
import gzip

parser = argparse.ArgumentParser(

    description="Normalize an OTU table by marker gene copy number",

    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-i', '--input', metavar='PATH', required=True, type=str,
                    help='The input study table in biom format')

parser.add_argument('-o', '--output', metavar='PATH', required=True, type=str,
                    help='The output normalized table in biom format')

parser.add_argument('-c', '--input_count', metavar='PATH', required=True,
                    type=str,
                    help='Table of predicted marker gene copy numbers for ' +
                         'study sequences (can be gzipped)')

parser.add_argument('--metadata_identifer', default='CopyNumber',
                    help='Identifier for copy number entry as observation ' +
                         'metadata')

parser.add_argument('--load_precalc_file_in_biom', default=False,
                    action="store_true",
                    help='Instead of loading the precalculated file in ' +
                         'tab-delimited format (with sequence ids as row ' +
                         'ids and traits as columns) load the data in BIOM ' +
                         'format (with sequences as SampleIds and traits as ' +
                         'ObservationIds)')

parser.add_argument('--verbose', default=False, action="store_true",
                    help='Details will be printed to screen while program ' +
                         'runs')

def main():

    args = parser.parse_args()

    otu_table = load_table(args.input)

    ids_to_load = otu_table.ids(axis='observation')

    input_count_table = args.input_count

    if args.verbose:
        print("Loading trait table: " + input_count_table)

    count_table = load_table(input_count_table).transpose()

    #Need to only keep data relevant to our otu list
    ids = []
    for x in otu_table.iter(axis='observation'):
        ids.append(str(x[1]))

    ob_id = count_table.ids(axis='observation')[0]

    filtered_otus=[]
    filtered_values=[]
    for x in ids:
        if count_table.exists(x, axis='sample'):
            filtered_otus.append(x)
            filtered_values.append(otu_table.data(x, axis='observation'))

    filtered_otu_table = Table(filtered_values, filtered_otus, otu_table.ids())

    copy_numbers_filtered={}
    for x in filtered_otus:
        value = count_table.get_value_by_ids(ob_id,x)
        try:
            #data can be floats so round them and make them integers
            value = int(round(float(value)))

        except ValueError:
            raise ValueError("Invalid type passed as copy number for OTU ID %s. Must be int-able." % (value))
        if value < 1:
            raise ValueError("Copy numbers must be greater than or equal to 1.")

        copy_numbers_filtered[x]={args.metadata_identifer:value}

    filtered_otu_table.add_metadata(copy_numbers_filtered, axis='observation')

    def metadata_norm(v, i, md):
        return v / float(md[args.metadata_identifer])
    normalized_table = filtered_otu_table.transform(metadata_norm, axis='observation')
 
    #move Observation Metadata from original to filtered OTU table
    normalized_table = transfer_observation_metadata(otu_table, normalized_table, 'observation')

    make_output_dir_for_file(args.output)
    write_biom_table(normalized_table, args.output)


if __name__ == "__main__":
    main()
