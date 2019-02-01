#!/usr/bin/env python

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.1.0-b"

import argparse
from os import path
import sys
import time
from picrust2.default import (default_ref_dir, default_tables, default_map,
                              default_regroup_map, default_pathway_map)
from picrust2.util import make_output_dir
from picrust2.pipeline import full_pipeline

HSP_METHODS = ['mp', 'emp_prob', 'pic', 'scp', 'subtree_average']

parser = argparse.ArgumentParser(

    description="Wrapper for full PICRUSt2 pipeline. Will output "
                "predictions for E.C. numbers, KEGG orthologs, and MetaCyc "
                "pathway abundances and coverages by default. Note that this "
                "is a convenience script and you can also run each step "
                "individually. Descriptions of gene families and pathways "
                "will be added automatically to output files unless custom "
                "trait tables are input.",

    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-s', '--study_fasta', metavar='PATH', required=True,
                    type=str, help='FASTA of unaligned study sequences')

parser.add_argument('-i', '--input', metavar='PATH', required=True, type=str,
                    help='Input table of sequence abundances (BIOM or TSV '
                         'format)')

parser.add_argument('-o', '--output', metavar='PATH', required=True,
                    type=str, help='Output folder for final files')

parser.add_argument('--threads', type=int, default=1,
                    help='Number of threads to use (default: %(default)d).')

parser.add_argument('-r', '--ref_dir', metavar='PATH', type=str,
                    default=default_ref_dir,
                    help='Directory containing reference sequence files '
                         '(default: %(default)s). Please see the online '
                         'documentation for how to name the files in this '
                         'directory.')

parser.add_argument('--in_traits', type=str.upper, default='EC,KO',
                    help='Comma-delimited list (with no spaces) of which gene '
                         'families to predict from this set: COG, EC, KO, '
                         'PFAM, TIGRFAM. Note that E.C. numbers will always '
                         'be predicted unless --no_pathways is set '
                         '(default: %(default)s).')

parser.add_argument('--custom_trait_tables', type=str, metavar='PATH', default=None,
                    help='Optional path to custom trait tables with gene families '
                         'as columns and genomes as rows (overrides '
                         '--in_traits setting) to be used for hidden-state '
                         'prediction. Multiple tables can be specified by '
                         'delimiting filenames by commas. Importantly, the '
                         'first custom table specified will be used for '
                         'inferring pathway abundances. Typically this command '
                         'would be used with a custom marker gene table '
                         '(--marker_gene_table) as well')

parser.add_argument('--marker_gene_table', type=str, metavar='PATH', 
                    default=default_tables["16S"],
                    help='Path to marker gene copy number table (16S copy '
                         'numbers by default).')

parser.add_argument('--pathway_map', metavar='MAP', type=str,
                    default=default_pathway_map,
                    help='MinPath mapfile. The default mapfile maps MetaCyc '
                         'reactions to prokaryotic pathways '
                         '(default: %(default)s).')

parser.add_argument('--no_pathways', default=False, action='store_true',
                    help='Flag to indicate that pathways should NOT be '
                         'inferred (otherwise they will be inferred by '
                         'default). Predicted E.C. number abundances are used '
                         'to infer pathways when default reference files are '
                         'used.')

parser.add_argument('--regroup_map', metavar='ID_MAP',
                    default=default_regroup_map, type=str,
                    help='Mapfile of ids to regroup gene families to before '
                         'running MinPath. The default mapfile is for '
                         'regrouping E. C. numbers to MetaCyc reactions '
                         '(default: %(default)s).')

parser.add_argument('--no_regroup', default=False, action="store_true",
                    help='Do not regroup input gene families to reactions '
                         'as specified in the regrouping mapfile. This option '
                         'should only be used if you are using custom '
                         'reference and/or mapping files.')

parser.add_argument('--stratified', default=False, action='store_true',
                    help='Flag to indicate that stratified tables should be '
                         'generated at all steps (will increase run-time)')

parser.add_argument('--max_nsti', metavar='FLOAT', type=float, default=2.0,
                    help='Sequences with NSTI values above this value will '
                         'be excluded (default: %(default)d).')

parser.add_argument('--min_reads', metavar='INT', type=int, default=1,
                    help='Minimum number of reads across all samples for '
                         'each input ASV. ASVs below this cut-off will be '
                         'counted as part of the \"RARE\" category in the '
                         'stratified output (default: %(default)d).')

parser.add_argument('--min_samples', metavar='INT', type=int, default=1,
                    help='Minimum number of samples that an ASV needs to be '
                         'identfied within. ASVs below this cut-off will be '
                         'counted as part of the \"RARE\" category in the '
                         'stratified output (default: %(default)d).')

parser.add_argument('-m', '--hsp_method', default='mp',
                    choices=HSP_METHODS,
                    help='HSP method to use.' +
                    '"mp": predict discrete traits using max parsimony. '
                    '"emp_prob": predict discrete traits based on empirical '
                    'state probabilities across tips. "subtree_average": '
                    'predict continuous traits using subtree averaging. '
                    '"pic": predict continuous traits with phylogentic '
                    'independent contrast. "scp": reconstruct continuous '
                    'traits using squared-change parsimony (default: '
                    '%(default)s).')

parser.add_argument('--skip_nsti', default=False, action='store_true',
                    help='Do not calculate nearest-sequenced taxon index '
                    '(NSTI).')

parser.add_argument('--skip_minpath', default=False, action="store_true",
                    help='Do not run MinPath to identify which pathways are '
                         'present as a first pass (on by default).')

parser.add_argument('--no_gap_fill', default=False, action="store_true",
                    help='Do not perform gap filling before predicting '
                         'pathway abundances (Gap filling is on otherwise by '
                         'default.')

parser.add_argument('--coverage', default=False, action="store_true",
                    help='Calculate pathway coverages as well as abundances, '
                         'which are experimental and only useful for '
                         'advanced users.')

parser.add_argument('--per_sequence_contrib', default=False,
                    action="store_true",
                    help='Flag to specify that MinPath is run on the genes '
                    'contributed by each sequence (i.e. a predicted '
                    'genome) individually. Note this will greatly increase '
                    'the runtime. The output will be the predicted pathway '
                    'abundance contributed by each individual sequence. This '
                    'is in contrast to the default stratified output, which '
                    'is the contribution to the community-wide pathway '
                    'abundances. Pathway coverage stratified by contributing '
                    'sequence will also be output when --coverage is set '
                    '(default: %(default)s).')

parser.add_argument('--verbose', default=False, action='store_true',
                    help='If specified, print out wrapped commands to screen')

parser.add_argument('-v', '--version', default=False, action='version',
                    version="%(prog)s " + __version__)


def main():

    # Get start time.
    start_time = time.time()

    args = parser.parse_args()

    func_outfiles, pathway_outfiles = full_pipeline(study_fasta=args.study_fasta,
                                                    input_table=args.input,
                                                    output_folder=args.output,
                                                    threads=args.threads,
                                                    ref_dir=args.ref_dir,
                                                    in_traits=args.in_traits,
                                                    custom_trait_tables=args.custom_trait_tables,
                                                    marker_gene_table=args.marker_gene_table,
                                                    pathway_map=args.pathway_map,
                                                    no_pathways=args.no_pathways,
                                                    regroup_map=args.regroup_map,
                                                    skip_minpath=args.skip_minpath,
                                                    no_regroup=args.no_regroup,
                                                    coverage=args.coverage,
                                                    stratified=args.stratified,
                                                    max_nsti=args.max_nsti,
                                                    min_reads=args.min_reads,
                                                    min_samples=args.min_samples,
                                                    hsp_method=args.hsp_method,
                                                    skip_nsti=args.skip_nsti,
                                                    no_gap_fill=args.no_gap_fill,
                                                    per_sequence_contrib=args.per_sequence_contrib,
                                                    verbose=args.verbose)

    # Print out elapsed time.
    elapsed_time = time.time() - start_time
    print("Completed PICRUSt2 pipeline in " + "%.2f" % elapsed_time + " seconds.")

if __name__ == "__main__":
    main()
