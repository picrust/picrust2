#!/usr/bin/env python

__author__ = "Gavin Douglas"
__copyright__ = "Copyright 2018, The PICRUSt Project"
__credits__ = ["Gavin Douglas", "Morgan Langille"]
__license__ = "GPL"
__version__ = "2-alpha.2"

from cogent.util.option_parsing import parse_command_line_parameters, make_option
from picrust.wrap_hsp import castor_hsp_loocv_wrapper

script_info = {}
script_info['brief_description'] = "Leave-one-genome out validation of HSP"

script_info['script_description'] = "Given a set of tip names in tree " +\
                                   "will run HSP when leaving that tip out." +\
                                   "Expected and predicted values will be " +\
                                   "output along with the Spearman rho " +\
                                   "between the expected and observed and " +\
                                   "the NSTI value for that tip. Currently " +\
                                   "only HSP with maximum parsimony is enabled"

script_info['required_options'] = [

  make_option('-i', '--observed_trait_table', type="existing_filepath",
              help='the input trait table describing directly observed ' +
                   'traits (e.g. sequenced genomes) in tab-delimited format'),

  make_option('-t', '--tree', type="existing_filepath",
              help='the full reference tree, in newick format'),

  make_option('-n', '--names', type="existing_filepath",
              help='File with tip names to leave out - one per line')
]

script_info['optional_options'] = [

  make_option('--exp_out', type="new_filepath",
              default='expected_traits.tsv',
              help='the output filepath for expected trait values ' +
                   '[default: %default]'),

  make_option('--pred_out', type="new_filepath",
              default='expected_traits.tsv',
              help='the output filepath for predicted trait values ' +
                   '[default: %default]'),

  make_option('-m', '--metrics_out', type="new_filepath",
              default='rho_nsti_out.tsv',
              help='the output filepath for Rho and NSTI values per left out genome ' +
                   '[default: %default]'),

  make_option('-p', '--processes', default=1, type="int",
              help='Number of processes to run in parallel.' +
                   '[default: %default]'),

  make_option('--debug', default=False, action="store_true",
              help='Flag to specify run in debugging mode. ' +
                   '[default: %default]')
]

script_info['version'] = __version__


def main():

    option_parser, opts, args = parse_command_line_parameters(**script_info)

    castor_hsp_loocv_wrapper(tree_path=opts.tree,
                             trait_table_path=opts.observed_trait_table,
                             tips_path=opts.names,
                             hsp_method="mp",
                             expected_out_path=opts.exp_out,
                             predicted_out_path=opts.pred_out,
                             metrics_out_path=opts.metrics_out,
                             num_cores=opts.processes,
                             HALT_EXEC=opts.debug)

if __name__ == "__main__":
    main()
