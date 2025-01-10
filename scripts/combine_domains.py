#!/usr/bin/env python

import argparse
from importlib.metadata import version
from picrust2.split_domains import combine_domain_predictions

parser = argparse.ArgumentParser(

    description="This script combined the functional predictions for two different domains into a single table. "
                "Typically this script would be called after performing hidden state prediction separately "
                "for two domains. The output can be used for the metagenome_pipeline.",
    epilog='''
Usage example:
combine_domains.py --table_dom1 arc_ec_predicted.tsv.gz --table_dom2 bac_ec_predicted.tsv.gz -o combined_ec_predicted.tsv.gz'
''', formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('--table_dom1', metavar='PATH', required=True, type=str,
                    help='The input file of predictions for domain1 in tab-delimited format.')

parser.add_argument('--table_dom2', metavar='PATH', required=True, type=str,
                    help='The input file of predictions for domain2 in tab-delimited format.')
                         
parser.add_argument('-o', '--output', metavar='PATH', required=True, type=str,
                    help='Output table with predicted abundances per study sequence in input files. '
                        'If the extension \".gz\" is added the table will automatically be gzipped.')
                        
parser.add_argument('--verbose', default=False, action='store_true',
                    help='If specified, print out wrapped commands and other '
                         'details to screen.')

def main():

    args = parser.parse_args()
    
    # Combine the domains
    combine_domain_predictions(args.table_dom1, 
                              args.table_dom2, 
                              args.output,
                              verbose=args.verbose)
    

if __name__ == "__main__":
    main()
