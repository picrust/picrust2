#!/usr/bin/env python

__copyright__ = "Copyright 2018-2021, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.4.2"

import argparse
from picrust2.pathway_pipeline import pathway_pipeline
from picrust2.util import (make_output_dir, check_files_exist,
                           TemporaryDirectory)
from picrust2.default import default_regroup_map, default_pathway_map
from os import path

parser = argparse.ArgumentParser(

    description="Script to infer the presence and abundances of pathways "
    "based on gene family abundances in a sample. By default, this script "
    "expects a table of EC number abundances (as output by PICRUSt2). "
    "However, alternative reaction to pathways mapping files can also be "
    "specified. By default, EC numbers are first regrouped to MetaCyc "
    "reactions, which are then linked to MetaCyc pathways through the default "
    "database.\n\n\n"

    "Stratified output will only be output if a stratified metagenome is "
    "input (or if --per_sequence_contrib is set). "
    "Please note that by default stratified abundances are based on "
    "how much predicted genomes (e.g. sequences) contribute to the "
    "community-wide abundance, not the abundance of the pathway based on the "
    "predicted genes in that genome alone. In other words, a predicted genome "
    "might be contributing greatly to the community-wide pathway abundance "
    "simply because one required gene for that pathway is at extremely high "
    "abundance in that genome even though no other required genes for that "
    "pathway are present. In contrast, the --per_sequence_contrib option "
    "should be used to get the predicted abundance and coverage of each "
    "pathway based on the predicted gene families within each genome. Note "
    "that using the --per_sequence_contrib option can greatly increase "
    "runtime.",
    epilog='''
    Usage examples:

    Default mapping of predicted EC number abundances to MetaCyc pathways:
    pathway_pipeline.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -o pathways_out --processes 1

    Mapping predicted KO abundances to legacy KEGG pathways (with stratified output that represents contributions to community-wide abundances):
    pathway_pipeline.py -i KO_metagenome_out/pred_metagenome_strat.tsv.gz -o KEGG_pathways_out --no_regroup --map picrust2/picrust2/default_files/pathway_mapfiles/KEGG_pathways_to_KO.tsv

    Map EC numbers to MetaCyc pathways and get stratified output corresponding to contribution of predicted gene family abundances within each predicted genome:
    pathway_pipeline.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -o pathways_out_per_seq --per_sequence_contrib --per_sequence_abun EC_metagenome_out/seqtab_norm.tsv.gz --per_sequence_function EC_predicted.tsv.gz

''', formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-i', '--input', metavar='IN_TABLE', required=True,
                    type=str,
                    help='Input TSV table of gene family abundances (either '
                         'the unstratified or stratified output of '
                         'metagenome_pipeline.py).')

parser.add_argument('-o', '--out_dir', metavar='DIRECTORY', required=True,
                    type=str, help='Output folder for pathway abundance output.')

parser.add_argument('-m', '--map', metavar='MAP', type=str,
                    default=default_pathway_map,
                    help='Mapping of pathways to reactions. The default '
                          'mapfile maps MetaCyc reactions to prokaryotic '
                          'MetaCyc pathways (default: %(default)s).')

parser.add_argument('--skip_minpath', default=False, action="store_true",
                    help='Do not run MinPath to identify which pathways are '
                         'present as a first pass (on by default).')

parser.add_argument('--no_gap_fill', default=False, action="store_true",
                    help='Do not perform gap filling before predicting '
                         'pathway abundances (Gap filling is on otherwise by '
                         'default.')

parser.add_argument('--intermediate', metavar='DIR', type=str, default=None,
                    help='Output folder for intermediate files (will be '
                         'deleted otherwise).')

parser.add_argument('-p', '--processes', default=1, type=int,
                    help='Number of processes to run in parallel '
                         '(default: %(default)d).')

parser.add_argument('--no_regroup', default=False, action="store_true",
                    help='Do not regroup input gene families to reactions '
                         'as specified in the regrouping mapfile.')

parser.add_argument('--coverage', default=False, action="store_true",
                    help='Calculate pathway coverages as well as abundances, '
                         'which are experimental and only useful for '
                         'advanced users.')

parser.add_argument('-r', '--regroup_map', metavar='ID_MAP',
                    default=default_regroup_map, type=str,
                    help='Mapfile of ids to regroup gene families to before '
                         'running MinPath. The default mapfile is for '
                         'regrouping EC numbers to MetaCyc reactions '
                         '(default: %(default)s).')

parser.add_argument('--per_sequence_contrib', default=False,
                    action="store_true",
                    help='Flag to specify that MinPath is run on the genes '
                    'contributed by each sequence (i.e. a predicted '
                    'genome) individually. Note this will greatly increase '
                    'the runtime. The output will be the predicted pathway '
                    'abundance contributed by each individual sequence. This '
                    'is in contrast to the default stratified output, which '
                    'is the contribution to the community-wide pathway '
                    'abundances. Experimental pathway coverage stratified by contributing '
                    'sequence will also be output when --coverage is set. '
                    'Options --per_sequence_contrib and '
                    '--per_sequence_function need to be set when this option '
                    'is used (default: %(default)s).')

parser.add_argument('--per_sequence_abun', metavar='PATH',
                    default=None, type=str,
                    help='Path to table of sequence abundances across samples '
                         'normalized by marker copy number (typically the normalized '
                         'sequence abundance table output at the metagenome '
                         'pipeline step). This input is required when the '
                         '--per_sequence_contrib option is set. (default: '
                         '%(default)s).')

parser.add_argument('--per_sequence_function', metavar='PATH',
                    default=None, type=str,
                    help='Path to table of function abundances per sequence, '
                         'which was outputted at the hidden-state prediction '
                         'step. This input is required when the '
                         '--per_sequence_contrib option is set. Note that '
                         'this file should be the same input table as used '
                         'for the metagenome pipeline step (default: '
                         '%(default)s).')

parser.add_argument('--wide_table', default=False,
                    action="store_true",
                    help='Flag to specify that wide-format stratified table '
                         'should be output rather than metagenome '
                         'contribution table. This is the deprecated method '
                         'of generating stratified tables since it is '
                         'extremely memory intensive (default: %(default)s).')

parser.add_argument('--verbose', default=False, action='store_true',
                    help='If specified, print out wrapped commands and other '
                         'details to screen.')

parser.add_argument('-v', '--version', default=False, action='version',
                    version="%(prog)s " + __version__)


def main():

    args = parser.parse_args()

    # Check that input files exist.
    check_files_exist([args.input, args.map])

    gap_fill_opt = not args.no_gap_fill

    run_minpath_opt = not args.skip_minpath

    # If intermediate output directory set then create and output there.
    # Otherwise make a temporary directory for the intermediate files.
    if args.intermediate:

        make_output_dir(args.intermediate)

        unstrat_abun, \
            unstrat_cov, \
            strat_abun, \
            strat_cov, \
            path_abun_by_seq, \
            path_cov_by_seq, \
            unstrat_abun_per_seq = pathway_pipeline(
                            inputfile=args.input,
                            mapfile=args.map,
                            regroup_mapfile=args.regroup_map,
                            proc=args.processes,
                            out_dir=args.intermediate,
                            run_minpath=run_minpath_opt,
                            coverage=args.coverage,
                            gap_fill_on=gap_fill_opt,
                            no_regroup=args.no_regroup,
                            per_sequence_contrib=args.per_sequence_contrib,
                            per_sequence_abun=args.per_sequence_abun,
                            per_sequence_function=args.per_sequence_function,
                            wide_table=args.wide_table,
                            verbose=args.verbose)
    else:
        with TemporaryDirectory() as temp_dir:
            unstrat_abun, \
                unstrat_cov, \
                strat_abun, \
                strat_cov, \
                path_abun_by_seq, \
                path_cov_by_seq, \
                unstrat_abun_per_seq = pathway_pipeline(
                            inputfile=args.input,
                            mapfile=args.map,
                            regroup_mapfile=args.regroup_map,
                            proc=args.processes,
                            out_dir=temp_dir,
                            run_minpath=run_minpath_opt,
                            coverage=args.coverage,
                            gap_fill_on=gap_fill_opt,
                            no_regroup=args.no_regroup,
                            per_sequence_contrib=args.per_sequence_contrib,
                            per_sequence_abun=args.per_sequence_abun,
                            per_sequence_function=args.per_sequence_function,
                            wide_table=args.wide_table,
                            verbose=args.verbose)

    make_output_dir(args.out_dir)

    # Write output files. The unstratified abundance table will always be
    # written, but the other files will only be written if applicable.
    unstrat_abun_outfile = path.join(args.out_dir, "path_abun_unstrat.tsv.gz")
    unstrat_abun.to_csv(path_or_buf=unstrat_abun_outfile, sep="\t",
                        index_label="pathway", compression="gzip")

    if args.coverage:
        unstrat_cov_outfile = path.join(args.out_dir,
                                        "path_cov_unstrat.tsv.gz")
        unstrat_cov.to_csv(path_or_buf=unstrat_cov_outfile, sep="\t",
                           index_label="pathway", compression="gzip")

    if strat_abun is not None:

        if args.wide_table:
            strat_abun_outfile = path.join(args.out_dir,
                                           "path_abun_strat.tsv.gz")
        else:
            strat_abun_outfile = path.join(args.out_dir,
                                           "path_abun_contrib.tsv.gz")

        strat_abun.to_csv(path_or_buf=strat_abun_outfile, sep="\t",
                          index=False, compression="gzip")

    if args.coverage and strat_cov is not None:
        if args.wide_table:
            strat_cov_outfile = path.join(args.out_dir,
                                          "path_cov_strat.tsv.gz")
        else:
            strat_cov_outfile = path.join(args.out_dir,
                                          "path_cov_contrib.tsv.gz")

        strat_cov.to_csv(path_or_buf=strat_cov_outfile, sep="\t",
                         index=False, compression="gzip")

    if path_abun_by_seq is not None:
        genome_path_abun_outfile = path.join(args.out_dir,
                                             "path_abun_predictions.tsv.gz")
        path_abun_by_seq.to_csv(path_or_buf=genome_path_abun_outfile, sep="\t",
                                index=True, compression="gzip",
                                index_label="sequence")

    if args.coverage and path_cov_by_seq is not None:
        genome_path_cov_outfile = path.join(args.out_dir,
                                            "path_cov_predictions.tsv.gz")
        path_cov_by_seq.to_csv(path_or_buf=genome_path_cov_outfile, sep="\t",
                               index=True, compression="gzip",
                               index_label="sequence")

    if unstrat_abun_per_seq is not None:
        unstrat_abun_per_seq_outfile = path.join(args.out_dir,
                                                 "path_abun_unstrat_per_seq.tsv.gz")
        unstrat_abun_per_seq.to_csv(path_or_buf=unstrat_abun_per_seq_outfile,
                                    sep="\t", index_label="pathway",
                                    compression="gzip")


if __name__ == "__main__":
    main()
