#!/usr/bin/env python

import argparse
from importlib.metadata import version
import sys
import time
from picrust2.default import (default_ref_dir, default_tables,
                              default_regroup_map, default_pathway_map)
from picrust2.util import restricted_float
from picrust2.pipeline import full_pipeline

HSP_METHODS = ['mp', 'emp_prob', 'pic', 'scp', 'subtree_average']

parser = argparse.ArgumentParser(

    description="Wrapper for full PICRUSt2 pipeline. Run sequence placement "
                "with EPA-NG and GAPPA to place study sequences (i.e. OTUs "
                "and ASVs) into a reference tree. Then runs hidden-state "
                "prediction with the castor R package to predict genome for "
                "each study sequence. Metagenome profiles are then generated, "
                "which can be optionally stratified by the contributing "
                "sequence. Finally, pathway abundances are predicted based on "
                "metagenome profiles. By default, output files include "
                "predictions for Enzyme Commission (EC) numbers, "
                "KEGG Orthologs (KOs), and MetaCyc "
                "pathway abundances. However, this script enables users to "
                "use custom reference and trait tables to customize analyses.",

    epilog='''
Run full default pipeline with 10 cores (only unstratified output):
picrust2_pipeline.py -s study_seqs.fna -i seqabun.biom -o picrust2_out --processes 10

Run full default pipeline with 10 cores with stratified output (including pathway stratified output based on per-sequence contributions):
picrust2_pipeline.py -s study_seqs.fna -i seqabun.biom -o picrust2_out --processes 10 --stratified --per_sequence_contrib''',

    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-s', '--study_fasta', metavar='PATH', required=True,
                    type=str, help='FASTA of unaligned study sequences (e.g. '
                                   'ASVs). The headerline should be only one field (i.e. no additional whitespace-delimited fields).')

parser.add_argument('-i', '--input', metavar='PATH', required=True, type=str,
                    help='Input table of sequence abundances (BIOM, TSV or '
                         'mothur shared file format).')

parser.add_argument('-o', '--output', metavar='PATH', required=True,
                    type=str, help='Output folder for final files.')

parser.add_argument('-p', '--processes', type=int, default=1,
                    help='Number of processes to run in parallel (default: '
                         '%(default)d).')

parser.add_argument('-t', '--placement_tool', metavar='epa-ng|sepp',
                    choices=['epa-ng', 'sepp'], default="epa-ng",
                    help='Placement tool to use when placing sequences into '
                         'reference tree. One of \"epa-ng\" or \"sepp\" '
                         'must be input (default: %(default)s)')

parser.add_argument('-r', '--ref_dir', metavar='PATH', type=str,
                    default=default_ref_dir,
                    help='Directory containing reference sequence files '
                         '(default: %(default)s). Please see the online '
                         'documentation for how to name the files in this '
                         'directory.')

parser.add_argument('--in_traits', type=str.upper, default='EC,KO',
                    help='Comma-delimited list (with no spaces) of which gene '
                         'families to predict from this set: COG, EC, KO, '
                         'PFAM, TIGRFAM. Note that EC numbers will always '
                         'be predicted unless --no_pathways is set '
                         '(default: %(default)s).')

parser.add_argument('--custom_trait_tables', type=str, metavar='PATH',
                    default=None,
                    help='Optional path to custom trait tables with gene '
                         'families as columns and genomes as rows (overrides '
                         '--in_traits setting) to be used for hidden-state '
                         'prediction. Multiple tables can be specified by '
                         'delimiting filenames by commas. Importantly, the '
                         'first custom table specified will be used for '
                         'inferring pathway abundances. Typically this '
                         'command would be used with a custom marker gene '
                         'table (--marker_gene_table) as well.')

parser.add_argument('--marker_gene_table', type=str, metavar='PATH',
                    default=default_tables["16S"],
                    help='Path to marker gene copy number table (16S copy '
                         'numbers by default).')

parser.add_argument('--pathway_map', metavar='MAP', type=str,
                    default=default_pathway_map,
                    help='MinPath mapfile. The default mapfile maps MetaCyc '
                         'reactions to prokaryotic pathways '
                         '(default: %(default)s).')

parser.add_argument('--reaction_func', metavar='MAP', type=str, default="EC",
                    help='Functional database to use as reactions for '
                         'inferring pathway abundances (default: '
                         '%(default)s). This should be either the short-form '
                         'of the database as specified in --in_traits, or the '
                         'path to the file as would be specified for '
                         '--custom_trait_tables. Note that when functions '
                         'besides the default EC numbers are used typically '
                         'the --no_regroup option would also be set.')

parser.add_argument('--no_pathways', default=False, action='store_true',
                    help='Flag to indicate that pathways should NOT be '
                         'inferred (otherwise they will be inferred by '
                         'default). Predicted EC number abundances are used '
                         'to infer pathways when the default reference files are '
                         'used.')

parser.add_argument('--regroup_map', metavar='ID_MAP',
                    default=default_regroup_map, type=str,
                    help='Mapfile of ids to regroup gene families to before '
                         'running MinPath. The default mapfile is for '
                         'regrouping EC numbers to MetaCyc reactions '
                         '(default: %(default)s).')

parser.add_argument('--no_regroup', default=False, action="store_true",
                    help='Do not regroup input gene families to reactions '
                         'as specified in the regrouping mapfile. This option '
                         'should only be used if you are using custom '
                         'reference and/or mapping files.')

parser.add_argument('--stratified', default=False, action='store_true',
                    help='Flag to indicate that stratified tables should be '
                         'generated at all steps (will increase run-time).')

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

parser.add_argument('-e', '--edge_exponent', default=0.5, type=float,
                    help='Setting for maximum parisomony hidden-state '
                          'prediction. Specifies weighting transition costs '
                          'by the inverse length of edge lengths. If 0, then '
                          'edge lengths do not influence predictions. Must be '
                          'a non-negative real-valued number (default: '
                          '%(default)f).')

parser.add_argument('--min_align', type=restricted_float, default=0.8,
                    help='Proportion of the total length of an input query '
                         'sequence that must align with reference sequences. '
                         'Any sequences with lengths below this value after '
                         'making an alignment with reference sequences will '
                         'be excluded from the placement and all subsequent '
                         'steps. (default: %(default)d).')

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
                         'genome) individually. Note this will greatly '
                         'increase the runtime. The output will be the '
                         'predicted pathway abundance contributed by each '
                         'individual sequence. This is in contrast to the '
                         'default stratified output, which is the contribution '
                         'to the community-wide pathway abundances. Pathway '
                         'coverage stratified by contributing sequence will '
                         'also be output when --coverage is set (default: '
                         '%(default)s).')

parser.add_argument('--wide_table', default=False, action='store_true',
                    help='Output wide-format stratified table of metagenome '
                         'and pathway predictions when \"--stratified\" is '
                         'set. This is the deprecated method of generating '
                         'stratified tables since it is extremely memory '
                         'intensive. The stratified filenames contain '
                         '\"strat\" rather than \"contrib\" when this option '
                         'is used.')

parser.add_argument('--skip_norm', default=False, action='store_true',
                    help='Skip normalizing sequence abundances by predicted '
                         'marker gene copy numbers (typically 16S rRNA '
                         'genes). This step will be performed automatically '
                         'unless this option is specified.')

parser.add_argument('--remove_intermediate', default=False,
                    action='store_true',
                    help='Remove the intermediate outfiles of the sequence '
                         'placement and pathway inference steps.')

parser.add_argument('--verbose', default=False, action='store_true',
                    help='Print out details as commands are run.')

parser.add_argument('-v', '--version', default=False, action='version',
                    version="PICRUSt2 " + version('PICRUSt2'))


def main():

    start_time = time.time()

    args = parser.parse_args()

    func_outfiles, pathway_outfiles = full_pipeline(study_fasta=args.study_fasta,
                                                    input_table=args.input,
                                                    output_folder=args.output,
                                                    processes=args.processes,
                                                    placement_tool=args.placement_tool,
                                                    ref_dir=args.ref_dir,
                                                    in_traits=args.in_traits,
                                                    custom_trait_tables=args.custom_trait_tables,
                                                    marker_gene_table=args.marker_gene_table,
                                                    pathway_map=args.pathway_map,
                                                    rxn_func=args.reaction_func,
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
                                                    edge_exponent=args.edge_exponent,
                                                    min_align=args.min_align,
                                                    skip_nsti=args.skip_nsti,
                                                    no_gap_fill=args.no_gap_fill,
                                                    per_sequence_contrib=args.per_sequence_contrib,
                                                    wide_table=args.wide_table,
                                                    skip_norm=args.skip_norm,
                                                    remove_intermediate=args.remove_intermediate,
                                                    verbose=args.verbose)

    if args.verbose:
        elapsed_time = time.time() - start_time
        print("Completed PICRUSt2 pipeline in " + "%.2f" % elapsed_time +
              " seconds.", file=sys.stderr)


if __name__ == "__main__":
    main()
