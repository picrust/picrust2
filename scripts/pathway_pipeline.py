#!/usr/bin/env python

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.1.0-b"

import argparse
from picrust2.pathway_pipeline import pathway_pipeline
from picrust2.util import make_output_dir, check_files_exist, TemporaryDirectory
from picrust2.default import default_regroup_map, default_pathway_map
from os import path

parser = argparse.ArgumentParser(

    description=
    "Wrapper for MinPath to infer which pathways are present "
    "given gene family abundances in a sample and to calculate the abundance "
    "and coverage of these pathways. By default this script expects a table "
    "of E.C. number abundances (as output by PICRUSt2). By default, these E.C. "
    "numbers will first be regrouped to MetaCyc reactions, which are then "
    "linked to MetaCyc pathways through the default database.\n\n\n"

    "Stratified output will only be output if a stratified metagenome is "
    "input. Note that by default, STRATIFIED "
    "ABUNDANCES ARE BASED ON HOW MUCH THAT PREDICTED GENOME "
    "(E.G. SEQUENCE) CONTRIBUTES TO THE COMMUNITY-WIDE ABUNDANCE, NOT THE "
    "ABUNDANCE OF THE PATHWAY BASED ON THE PREDICTED GENES IN "
    "THAT GENOME ALONE. In other words, a predicted genome might "
    "be contributing a lot to the community-wide pathway "
    "abundance simply because one required gene for that pathway "
    "is at extremely high abundance in that genome even though no "
    "other required genes for that pathway are present. In contrast, if you "
    "use the --per_sequence_contrib option you will get the predicted abundance and "
    "coverage of each pathway based on the predicted gene families WITHIN each "
    "genome. Note that using this option will likely greatly increase runtime.",

    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-i', '--input', metavar='IN_TABLE', required=True,
                    type=str,
                    help='Input TSV table of gene family abundances')

parser.add_argument('-o', '--out_dir', metavar='DIRECTORY', required=True,
                    type=str, help='Output folder for pathway abundance output')

parser.add_argument('-m', '--map', metavar='MAP', type=str,
                    default=default_pathway_map, 
                    help='MinPath mapfile. The default mapfile maps MetaCyc '
                         'reactions to prokaryotic pathways '
                         '(default: %(default)s).')

parser.add_argument('--skip_minpath', default=False, action="store_true",
                    help='Do not run MinPath to identify which pathways are '
                         'present as a first pass (on by default).')

parser.add_argument('--no_gap_fill', default=False, action="store_true",
                    help='Do not perform gap filling before predicting ' +
                         'pathway abundances (Gap filling is on otherwise by ' +
                         'default.')

parser.add_argument('--intermediate', metavar='DIR', type=str, default=None,
                    help='Output folder for intermediate files (wont be ' +
                         'kept unless this option is set.')

parser.add_argument('-p', '--proc', default=1, type=int,
                    help='Number of processes to run (default: %(default)d).')

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
                         'regrouping E.C. numbers to MetaCyc reactions '
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
                    'abundances. Pathway coverage stratified by contributing '
                    'sequence will also be output when --coverage is set. '
                    'Options --per_sequence_contrib and '
                    '--per_sequence_function need to be set when this option '
                    'is used (default: %(default)s).')

parser.add_argument('--per_sequence_abun', metavar='PATH',
                    default=None, type=str,
                    help='Path to table of sequence abundances across samples '
                         'normalized by marker copy number (i.e. normalized '
                         'sequence abundance table outputted by metagenome '
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

parser.add_argument('--print_cmds', default=False, action="store_true",
                    help='If specified, print out wrapped commands to screen')

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

        unstrat_abun, unstrat_cov, strat_abun, strat_cov = pathway_pipeline(
                                                      inputfile=args.input,
                                                      mapfile=args.map,
                                                      regroup_mapfile=args.regroup_map,
                                                      proc=args.proc,
                                                      out_dir=args.intermediate,
                                                      run_minpath=run_minpath_opt,
                                                      coverage=args.coverage,
                                                      gap_fill_on=gap_fill_opt,
                                                      no_regroup=args.no_regroup,
                                                      per_sequence_contrib=args.per_sequence_contrib,
                                                      per_sequence_abun=args.per_sequence_abun,
                                                      per_sequence_function=args.per_sequence_function,
                                                      print_cmds=args.print_cmds)
    else:
        with TemporaryDirectory() as temp_dir:
            unstrat_abun, unstrat_cov, strat_abun, strat_cov = pathway_pipeline(
                                                            inputfile=args.input,
                                                            mapfile=args.map,
                                                            regroup_mapfile=args.regroup_map,
                                                            proc=args.proc,
                                                            out_dir=temp_dir,
                                                            run_minpath=run_minpath_opt,
                                                            coverage=args.coverage,
                                                            gap_fill_on=gap_fill_opt,
                                                            no_regroup=args.no_regroup,
                                                            per_sequence_contrib=args.per_sequence_contrib,
                                                            per_sequence_abun=args.per_sequence_abun,
                                                            per_sequence_function=args.per_sequence_function,
                                                            print_cmds=args.print_cmds)

    make_output_dir(args.out_dir)

    # Write output files.
    unstrat_abun_outfile = path.join(args.out_dir, "path_abun_unstrat.tsv")
    unstrat_abun.to_csv(path_or_buf=unstrat_abun_outfile,  sep="\t",
                       index_label="pathway")

    if args.coverage:
        unstrat_cov_outfile = path.join(args.out_dir, "path_cov_unstrat.tsv")
        unstrat_cov.to_csv(path_or_buf=unstrat_cov_outfile,  sep="\t",
                           index_label="pathway")

    # Write stratified output only if something besides None was returned.
    if strat_abun is not None:
        strat_abun_outfile = path.join(args.out_dir, "path_abun_strat.tsv")
        strat_abun.to_csv(path_or_buf=strat_abun_outfile,  sep="\t",
                          index=False)

    if args.coverage and strat_cov is not None:
        strat_cov_outfile = path.join(args.out_dir, "path_cov_strat.tsv")
        strat_cov.to_csv(path_or_buf=strat_cov_outfile,  sep="\t",
                         index=False)

if __name__ == "__main__":
    main()
