#!/usr/bin/env python

__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2.0.0-b.6"

import argparse
from picrust2.run_minpath import run_minpath_pipeline
from tempfile import TemporaryDirectory
from picrust2.util import make_output_dir, check_files_exist
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

parser.add_argument('-r', '--regroup_map', metavar='ID_MAP',
                    default=default_regroup_map, type=str,
                    help='Mapfile of ids to regroup gene families to before '
                         'running MinPath. The default mapfile is for '
                         'regrouping E.C. numbers to MetaCyc reactions '
                         '(default: %(default)s).')

parser.add_argument('--per_sequence_contrib', default=False, action="store_true",
                    help='Run MinPath on the gene families contributed by '
                    'each sequence (i.e. a predicted genome) individually. '
                    'This will only matter when a stratified table is input. '
                    'Note this will GREATLY increase the runtime, but will '
                    'output the predicted pathway abundance contributed by the '
                    'predicted gene families in each predicted genome alone '
                    '(i.e. not the contribution to the community-wide '
                    'abundance). Pathway coverage stratified by contributing '
                    'sequence will also be output when this option is set '
                    '(default: %(default)d).')

parser.add_argument('--print_cmds', default=False, action="store_true",
                    help='If specified, print out wrapped commands to screen')

parser.add_argument('-v', '--version', default=False, action='version',
                    version="%(prog)s " + __version__)

def main():

    args = parser.parse_args()

    # Check that input files exist.
    check_files_exist([args.input, args.map])

    gap_fill_opt = not args.no_gap_fill

    # If no regrouping flag set then set input regrouping mapfile to be None.
    if args.no_regroup:
        args.regroup_map = None

    # If intermediate output directory set then create and output there.
    # Otherwise make a temporary directory for the intermediate files.
    if args.intermediate:

        make_output_dir(args.intermediate)

        unstrat_abun, unstrat_cov, strat_abun, strat_cov = run_minpath_pipeline(
                                                      inputfile=args.input,
                                                      mapfile=args.map,
                                                      regroup_mapfile=args.regroup_map,
                                                      proc=args.proc,
                                                      out_dir=args.intermediate,
                                                      gap_fill=gap_fill_opt,
                                                      per_sequence_contrib=args.per_sequence_contrib,
                                                      print_cmds=args.print_cmds)
    else:
        with TemporaryDirectory() as temp_dir:
            unstrat_abun, unstrat_cov, strat_abun, strat_cov = run_minpath_pipeline(
                                                            inputfile=args.input,
                                                            mapfile=args.map,
                                                            regroup_mapfile=args.regroup_map,
                                                            proc=args.proc,
                                                            out_dir=temp_dir,
                                                            gap_fill=gap_fill_opt,
                                                            per_sequence_contrib=args.per_sequence_contrib,
                                                            print_cmds=args.print_cmds)

    make_output_dir(args.out_dir)

    # Write output files.
    unstrat_abun_outfile = path.join(args.out_dir, "path_abun_unstrat.tsv")
    unstrat_abun.to_csv(path_or_buf=unstrat_abun_outfile,  sep="\t",
                       index_label="pathway")

    unstrat_cov_outfile = path.join(args.out_dir, "path_cov_unstrat.tsv")
    unstrat_cov.to_csv(path_or_buf=unstrat_cov_outfile,  sep="\t",
                       index_label="pathway")

    # Write stratified output only if something besides None was returned.
    if strat_abun is not None:
        strat_abun_outfile = path.join(args.out_dir, "path_abun_strat.tsv")
        strat_abun.to_csv(path_or_buf=strat_abun_outfile,  sep="\t",
                          index=False)

    if strat_cov is not None:
        strat_cov_outfile = path.join(args.out_dir, "path_cov_strat.tsv")
        strat_cov.to_csv(path_or_buf=strat_cov_outfile,  sep="\t",
                         index=False)

if __name__ == "__main__":
    main()
