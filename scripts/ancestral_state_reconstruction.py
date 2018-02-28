#!/usr/bin/env python
# File created on 1 Feb 2012
from __future__ import division

__author__ = "Morgan Langille"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
__credits__ = ["Morgan Langille"]
__license__ = "GPL"
__version__ = "1.1.3"
__maintainer__ = "Morgan Langille"
__email__ = "morgan.g.i.langille@gmail.com"
__status__ = "Development"


from cogent.util.option_parsing import parse_command_line_parameters, make_option
from picrust.count import wagner_for_picrust
from picrust.wrap_asr import ape_ace_wrapper, castor_asr_pic_wrapper
from picrust.ancestral_state_reconstruction import run_asr_in_parallel
from picrust.util import make_output_dir_for_file, make_output_dir

script_info = {}
script_info['brief_description'] = "Runs ancestral state reconstruction " +\
                                    "given a tree and trait table"

script_info['script_description'] = "Provides a common interface for " +\
                                    "running various ancestral state " +\
                                    " reconstruction methods (e.g. PIC)."

script_info['script_usage'] = [("Example 1", 
                                "Provide a tree file and trait table file:",
                                "%prog -i trait_table.tab -t " +\
                                "pruned_tree.newick -o asr_counts.tab " +\
                                "-c asr_ci.tab")]

script_info['output_description']= "A table containing trait information " +\
                                   "for internal nodes of the tree."

script_info['required_options'] = [
    make_option('-t', '--input_tree_fp', type="existing_filepath", 
                help='the tree to use for ASR'),
    make_option('-i', '--input_trait_table_fp', type="existing_filepath", 
                help='the trait table to use for ASR')]

asr_method_choices = ['castor_pic', 'ace_ml', 'ace_reml', 'ace_pic', 'wagner']
parallel_method_choices = ['multithreaded', 'sge', 'torque']

script_info['optional_options'] = [
make_option('-m', '--asr_method', type='choice',
            help='Method for ancestral state reconstruction. Valid choices '+\
                 'are: ' + ', '.join(asr_method_choices) +\
                 ' [default: %default]',
            choices=asr_method_choices, default='ace_pic'),

make_option('-o', '--output_fp' ,type="new_filepath",
            help='output trait table [default:%default]', 
            default='asr_counts.tab'),

make_option('-c', '--output_ci_fp', type="new_filepath",
            help='output table containing 95% confidence intervals, loglik, ' +\
                 'and brownian motion parameters for each asr prediction ' +\
                 '[default:%default]', default='asr_ci.tab'),

make_option('--no_castor_ci', action="store_true",
            help='Do not calculate CI when running castor_pic (saves time)' +\
                 '[default:%default]', default=False),

make_option('-p', '--parallel', action="store_true",
            help='allow parallelization of asr', default=False),

make_option('-j','--parallel_method', type='choice',
            help='Method for parallelization. Valid choices are: ' +\
                 ', '.join(parallel_method_choices) + ' [default: %default]',
            choices=parallel_method_choices, default='sge'),

make_option('-n', '--num_jobs', action='store', type='int',
            help='Number of jobs to be submitted (if --parallel). ' +\
                 '[default: %default]',
            default=100),

make_option('-d', '--debug', action="store_true", 
            help='To aid with debugging; get the command that the app ' +\
                 'controller is going to run',
            default=False)
]

script_info['version'] = __version__
       

def main():
    option_parser, opts, args =\
                   parse_command_line_parameters(**script_info)
    
    if(opts.parallel):
        tmp_dir='jobs/'
        make_output_dir(tmp_dir)
        asr_table, ci_table =run_asr_in_parallel(tree=opts.input_tree_fp,
                                                 table=opts.input_trait_table_fp,
                                                 asr_method=opts.asr_method,
                                                 parallel_method=opts.parallel_method,
                                                 num_jobs=opts.num_jobs,
                                                 tmp_dir=tmp_dir,
                                                 verbose=opts.verbose)
    else:
        # Call the appropiate ASR app controller.
        if(opts.asr_method == 'castor_pic'):

            if opts.no_castor_ci:
                calc_castor_ci = False
            else:
                calc_castor_ci = True

            asr_table, ci_table = castor_asr_pic_wrapper(opts.input_tree_fp,
                                                         opts.input_trait_table_fp,
                                                         calc_ci=calc_castor_ci,
                                                         HALT_EXEC=opts.debug)
        if(opts.asr_method == 'wagner'):
            asr_table = wagner_for_picrust(opts.input_tree_fp,
                                           opts.input_trait_table_fp,
                                           HALT_EXEC=opts.debug)
        elif(opts.asr_method == 'ace_ml'):
            asr_table, ci_table = ape_ace_wrapper(opts.input_tree_fp,
                                                  opts.input_trait_table_fp,
                                                  'ML',
                                                  HALT_EXEC=opts.debug)
        elif(opts.asr_method == 'ace_pic'):
            asr_table, ci_table = ape_ace_wrapper(opts.input_tree_fp,
                                                  opts.input_trait_table_fp,
                                                  'pic',
                                                  HALT_EXEC=opts.debug)
        elif(opts.asr_method == 'ace_reml'):
            asr_table, ci_table = ape_ace_wrapper(opts.input_tree_fp,
                                                  opts.input_trait_table_fp,
                                                  'REML',
                                                  HALT_EXEC=opts.debug)

    # Output the table to file.
    make_output_dir_for_file(opts.output_fp)
    asr_table.writeToFile(opts.output_fp, sep='\t')

    # Output the CI file (unless the method is wagner 
    # or castor_pic and calc_castor_ci set to False).
    if not (opts.asr_method == 'wagner' or 
           (opts.asr_method == 'castor_pic' and not calc_castor_ci)):
        make_output_dir_for_file(opts.output_ci_fp)
        ci_table.writeToFile(opts.output_ci_fp, sep='\t')
        

if __name__ == "__main__":
    main()
