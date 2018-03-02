#!/usr/bin/env python
# File created on 15 Jan 2017
from __future__ import division

__author__ = "Gavin Douglas"
__copyright__ = "Copyright 2011-2017, The PICRUSt Project"
__credits__ = ["Gavin Douglas", "Morgan Langille"]
__license__ = "GPL"
__version__ = "2-alpha.1"
__maintainer__ = "Gavin Douglas"
__email__ = "gavinmdouglas@gmail.com"
__status__ = "Development"


from cogent.util.option_parsing import parse_command_line_parameters, make_option
from picrust.wrap_hsp import castor_hsp_wrapper
from picrust.util import make_output_dir_for_file

script_info = {}
script_info['brief_description'] = "Given a tree and a set of known " +\
                                   "character states will output " +\
                                   "predictions for unobserved character " +\
                                   "states. Note this does not require a " +\
                                   "separate ancestral state " +\
                                   "reconstruction step to be run."
script_info['script_description'] = "This script performs hidden state " +\
                                    "prediction on tips in the input tree " +\
                                    "with unknown trait values. Typically " +\
                                    "this script is used to predict the " +\
                                    "abundance of gene families present in " +\
                                    "each taxon, given a tree and a set of " +\
                                    "known trait values. This script " +\
                                    "outputs a table of trait predictions. " +\
                                    "Note that this script assumes that the input " +\
                                    "trait values will include \"0\" counts."

HSP_METHODS = ['emp_prob', 'mk_model', 'mp', 'pic', 'scp', 'subtree_average']

# Define command-line interface
script_info['output_description'] = "Output is a tab-delimited table of " +\
                                   "predicted character states"
script_info['required_options'] = [

  make_option('-i', '--observed_trait_table', type="existing_filepath",
              help='the input trait table describing directly observed ' +
                   'traits (e.g. sequenced genomes) in tab-delimited format'),

  make_option('-t', '--tree', type="existing_filepath",
              help='the full reference tree, in newick format')
]

script_info['optional_options'] = [

  make_option('-o', '--output_trait_table', type="new_filepath",
              default='predicted_traits.tsv',
              help='the output filepath for trait predictions ' +
                   '[default: %default]'),

  make_option('--ci_out', type="new_filepath",
              default='predicted_traits_ci.tsv',
              help='the output filepath for confidence intervals trait ' +
                   'predictions (if -c option set) [default: %default]'),

  make_option('-m', '--hsp_method', default='emp_prob', choices=HSP_METHODS,
              help='HSP method to use, options: ' +
                   ", ".join(HSP_METHODS) + '. "emp_prob": ' +
                   'predict discrete traits based on empirical state ' +
                   'probabilities across tips. "mk_model": predict ' +
                   'discrete traits based on fixed-rates continuous time ' +
                   'Markov model. "mp": predict discrete traits using max ' +
                   'parsimony. "subtree_average": predict continuous traits ' +
                   'using subtree averaging. "pic": predict continuous traits '+
                   'with phylogentic independent contrast. "scp": ' +
                   'reconstruct continuous traits using squared-change ' +
                   'parsimony [default: %default]'),

  make_option('-n', '--calculate_NSTI', default=False,
              action="store_true",
              help='if specified, calculate NSTI and add to output file ' +
                   '[default: %default]'),

  make_option('-c', "--confidence", default=False, action="store_true",
              help='if specified, output 95% confidence intervals (only ' +
                   'possible for mk_model, emp_prob, and mp settings) ' +
                   '[default: %default]'),

  make_option('--check', default=False, action="store_true",
              help='if specified, check input trait table before hsp ' +
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

    # Methods for discrete trait prediction with CI enabled.
    discrete_set = set(['emp_prob', 'mk_model', 'mp'])

    if opts.confidence and opts.hsp_method in discrete_set:
      ci_setting = True
    else:
      ci_setting = False

    hsp_table, ci_table = castor_hsp_wrapper(tree_path=opts.tree,
                                             trait_table_path=opts.observed_trait_table,
                                             hsp_method=opts.hsp_method,
                                             calc_nsti=opts.calculate_NSTI,
                                             calc_ci=ci_setting,
                                             check_input=opts.check,
                                             num_cores=opts.processes,
                                             HALT_EXEC=opts.debug)

    # Output the table to file.
    make_output_dir_for_file(opts.output_trait_table)
    hsp_table.writeToFile(opts.output_trait_table, sep='\t')

    # Output the CI file as well if option set.
    if (ci_setting):
        make_output_dir_for_file(opts.ci_out)
        ci_table.writeToFile(opts.ci_out, sep='\t')


if __name__ == "__main__":
    main()
