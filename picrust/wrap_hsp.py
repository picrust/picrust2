#!/usr/bin/env python

""" Application controller for HSP functions within the R package 'castor'.

"""

__author__ = "Gavin Douglas"
__copyright__ = "Copyright 2018, The PICRUSt Project"
__license__ = "GPL"
__version__ = "2-alpha.2"

from __future__ import division
from os import remove
from cogent.app.util import CommandLineApplication, get_tmp_filename
from cogent import LoadTable
from picrust.util import get_picrust_project_dir
from os.path import join

class Castor_hsp(CommandLineApplication):
    """ Application controller for HSP functions within 
    the 'castor' R package."""

    castor_hsp_script_fp = join(get_picrust_project_dir(), 'picrust', 
        'support_files', 'R', 'castor_hsp.R')

    _command = castor_hsp_script_fp
    _input_handler = '_input_as_string'
    _suppress_stdout = False
    _suppress_stderr = False

    # Overridden to call script with R rather than directly - this is useful
    # because permisssions on the script are set to 644 when PICRUSt is
    # installed with setup.py. This is fine if we're executing it with R, but
    # not if we're trying to execute it directly.
    def _get_base_command(self):
        """ Returns the full command string

            input_arg: the argument to the command which represents the input
                to the program, this will be a string, either
                representing input or a filename to get input from
         """
        command_parts = []
        # Append a change directory to the beginning of the command to change
        # to self.WorkingDir before running the command
        # WorkingDir should be in quotes -- filenames might contain spaces
        cd_command = ''.join(['cd ', str(self.WorkingDir), ';'])
        if self._command is None:
            raise ApplicationError, '_command has not been set.'
        command = self._command
        parameters = self.Parameters

        command_parts.append(cd_command)
        command_parts.append("R")
        command_parts.append("-f")
        command_parts.append(command)
        command_parts.append("--args")
        command_parts.append(self._command_delimiter.join(filter(
            None, (map(str, parameters.values())))))

        return self._command_delimiter.join(command_parts).strip()
    BaseCommand = property(_get_base_command)


class Castor_hsp_LOOCV(CommandLineApplication):
    """ Application controller for HSP functions within 
    the 'castor' R package when leaving genomes out."""

    castor_loocv_hsp_script_fp = join(get_picrust_project_dir(), 'picrust', 
        'support_files', 'R', 'castor_hsp_loocv.R')

    _command = castor_loocv_hsp_script_fp
    _input_handler = '_input_as_string'
    _suppress_stdout = False
    _suppress_stderr = False

    # Overridden to call script with R rather than directly - this is useful
    # because permisssions on the script are set to 644 when PICRUSt is
    # installed with setup.py. This is fine if we're executing it with R, but
    # not if we're trying to execute it directly.
    def _get_base_command(self):
        """ Returns the full command string

            input_arg: the argument to the command which represents the input
                to the program, this will be a string, either
                representing input or a filename to get input from
         """
        command_parts = []
        # Append a change directory to the beginning of the command to change
        # to self.WorkingDir before running the command
        # WorkingDir should be in quotes -- filenames might contain spaces
        cd_command = ''.join(['cd ', str(self.WorkingDir), ';'])
        if self._command is None:
            raise ApplicationError, '_command has not been set.'
        command = self._command
        parameters = self.Parameters

        command_parts.append(cd_command)
        command_parts.append("R")
        command_parts.append("-f")
        command_parts.append(command)
        command_parts.append("--args")
        command_parts.append(self._command_delimiter.join(filter(
            None, (map(str, parameters.values())))))

        return self._command_delimiter.join(command_parts).strip()
    BaseCommand = property(_get_base_command)


def castor_hsp_wrapper(tree_path,
                       trait_table_path,
                       hsp_method,
                       calc_nsti=False,
                       calc_ci=False,
                       check_input=False,
                       num_cores=1,
                       ran_seed=None,
                       HALT_EXEC=False):
                       
    '''Runs the Castor_hsp application controller given path of tree and 
    trait table and returns a Table'''

    # Initialize Castor_hsp app controller
    castor_hsp = Castor_hsp(HALT_EXEC=HALT_EXEC)

    tmp_output_count_path = get_tmp_filename()
    tmp_output_ci_path = get_tmp_filename()

    # Quote file names
    tree_path = '"{0}"'.format(tree_path)
    trait_table_path = '"{0}"'.format(trait_table_path)

    # Need to format boolean setting as string for R to read in as argument.
    if calc_nsti:
        calc_nsti_setting = "TRUE"
    else:
        calc_nsti_setting = "FALSE"

    if calc_ci:
        calc_ci_setting = "TRUE"
    else:
        calc_ci_setting = "FALSE"

    if check_input:
        check_input_setting = "TRUE"
    else:
        check_input_setting = "FALSE"

    as_string = " ".join([tree_path,
                          trait_table_path,
                          hsp_method,
                          calc_nsti_setting,
                          calc_ci_setting,
                          check_input_setting,
                          str(num_cores),
                          tmp_output_count_path,
                          tmp_output_ci_path,
                          str(ran_seed)])

    # Run castor_hsp here
    result = castor_hsp(data=as_string)

    # Load the output into Table objects
    try:
        asr_table = LoadTable(filename=tmp_output_count_path, header=True,
                              sep='\t')
    except IOError:
        raise RuntimeError,\
         ("R reported an error on stderr:" +\
          " %s" % "\n".join(result["StdErr"].readlines()))

    remove(tmp_output_count_path)

    if calc_ci:
        asr_ci_table = LoadTable(filename=tmp_output_ci_path, header=True,
                                 sep='\t')
        remove(tmp_output_ci_path)
    else:
        asr_ci_table = None

    return asr_table, asr_ci_table


def castor_hsp_loocv_wrapper(tree_path,
                             trait_table_path,
                             tips_path,
                             hsp_method,
                             expected_out_path,
                             predicted_out_path,
                             metrics_out_path,
                             num_cores=1,
                             HALT_EXEC=False):
                       
    '''Runs the Castor_hsp_LOOCV application controller given path of tree and 
    trait table and returns a Table'''

    # Initialize Castor_hsp_LOOCV app controller
    castor_hsp_loocv = Castor_hsp_LOOCV(HALT_EXEC=HALT_EXEC)

    # Quote file names
    tree_path = '"{0}"'.format(tree_path)
    trait_table_path = '"{0}"'.format(trait_table_path)
    tips_path = '"{0}"'.format(tips_path)
    expected_out_path = '"{0}"'.format(expected_out_path)
    predicted_out_path = '"{0}"'.format(predicted_out_path)
    metrics_out_path = '"{0}"'.format(metrics_out_path)

    as_string = " ".join([tree_path,
                          trait_table_path,
                          tips_path,
                          hsp_method,
                          expected_out_path,
                          predicted_out_path,
                          metrics_out_path,
                          str(num_cores)])

    # Run castor_hsp_loocv here
    result = castor_hsp_loocv(data=as_string)

    # Check if error returned and print out if so:
    result_stderr = result["StdErr"].readlines()
    if result_stderr:
        raise RuntimeError,\
         ("R reported an error on stderr:" +\
          " %s" % "\n".join(result_stderr))

