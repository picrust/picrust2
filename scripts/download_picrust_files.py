#!/usr/bin/env python
from __future__ import absolute_import, print_function
from future.moves.urllib.parse import urljoin
from future.moves.urllib.request import urlopen, HTTPError


__author__ = "The PICRUSt Development Team"
__copyright__ = "Copyright 2017, PICRUSt Project"
__credits__ = ["Morgan Langille", "Jesse Zaneveld", "Greg Caporaso",
               "Daniel McDonald", "Dan Knights", "Joshua Reyes",
               "Jose Clemente", "Rob Knight", "Rob Beiko",
               "Curtis Huttenhower"]
__license__ = "GPL"
__url__ = "http://picrust.github.com"
__version__ = "1.1.3"
__maintainer__ = "Morgan Langille"
__email__ = "morgan.g.i.langille@gmail.com"


from cogent.util.option_parsing import (
    make_option,
    parse_command_line_parameters,
)
import os
import picrust
from picrust.util import atomic_write


PROJECT_ROOT = os.path.dirname(os.path.abspath(picrust.__file__))
DATA_DIR = os.path.join(PROJECT_ROOT, 'data')

BASE_URL = (
    'http://kronos.pharmacology.dal.ca/public_files/picrust/'
    'picrust_precalculated_v1.1.2/'
)

type_of_prediction_choices = ['ko', 'cog', 'rfam']
gg_version_choices = ['13_5', '18may2012']

FILES = {
    ('16s', '13_5'): '16S_13_5_precalculated.tab.gz',
    ('ko', '13_5'): 'ko_13_5_precalculated.tab.gz',
    ('cog', '13_5'): 'cog_13_5_precalculated.tab.gz',
    ('rfam', '13_5'): 'rfam_13_5_precalculated.tab.gz',
    ('16s_ci', '13_5'): '16S_13_5_precalculated_variances.tab.gz',
    ('ko_ci', '13_5'): 'ko_13_5_precalculated_variances.tab.gz',
    ('cog_ci', '13_5'): 'cog_13_5_precalculated_variances.tab.gz',
    ('rfam_ci', '13_5'): 'rfam_13_5_precalculated_variances.tab.gz',
    ('16s', '18may2012'): '16S_18may2012_precalculated.tab.gz',
    ('ko', '18may2012'): 'ko_18may2012_precalculated.tab.gz',
    ('cog', '18may2012'): 'cog_18may2012_precalculated.tab.gz',
}

script_info = {
    'brief_description': 'Downloads PICRUSt pre-calculated files.',
    'script_description': (
        'Downloads PICRUSt pre-calculated files to the data directory'
        ' ({}).'.format(DATA_DIR)),
    'script_usage': [
        ('', 'Download default pre-calculated files:', '%prog')],
    'output_description': (
        'Prints the result of the download attempt to the screen (STDOUT).'),
    'required_options': [],
    'optional_options': [
        make_option(
            '-t', '--type_of_prediction',
            default=type_of_prediction_choices[0],
            type="choice",
            choices=type_of_prediction_choices,
            help='Type of functional predictions. Valid choices are:'
                 ' {choices} [default: %default]'.format(
                    choices=', '.join(type_of_prediction_choices))
        ),
        make_option(
            '-g', '--gg_version',
            default=gg_version_choices[0],
            type="choice",
            choices=gg_version_choices,
            help='Version of GreenGenes that was used for OTU picking. Valid'
                 ' choices are: {choices} [default: %default]'.format(
                    choices=', '.join(gg_version_choices))
        ),
        make_option(
            '--with_confidence',
            default=False,
            action="store_true",
            help='Download confidence interval files (only available for'
                 ' GreenGenes 13_5) [default: %default]'
        ),
        make_option(
            '--force',
            action="store_true",
            default=False,
            help=('Force download of files (i.e. overwrite existing) [default:'
                  ' %default]'))
    ],
    'help_on_no_arguments': False,
    'version': __version__,
}


def download_picrust_files(
        output_path, with_confidence, gg_version, type_of_prediction, force):
    if not os.path.exists(output_path):
        print('Creating output directory {path}...'.format(path=output_path))
        os.mkdir(output_path)
    elif not os.path.isdir(output_path):
        print('Error: {path} already exists but is not a directory!'.format(
            path=output_path))
        exit(1)

    try:
        files_to_download = [
            FILES[('16s', gg_version)],
            FILES[(type_of_prediction, gg_version)]
        ]
    except KeyError:
        msg = (
            'Error: {type_of_prediction} predictions not available for'
            ' GreenGenes {gg_version}')
        print(msg.format(
            gg_version=gg_version, type_of_prediction=type_of_prediction))
        exit(1)

    if with_confidence:
        try:
            ci_files = [
                FILES[('16s_ci', gg_version)],
                FILES[(type_of_prediction + '_ci', gg_version)],
            ]
        except KeyError:
            msg = 'Error: CI files not available for {gg_version}'
            print(msg.format(gg_version=gg_version))
            exit(1)

        files_to_download.extend(ci_files)

    print('Downloading files to {path}...'.format(path=output_path))

    for filename in files_to_download:
        if 'variances' not in filename:
            url = urljoin(BASE_URL, gg_version + '/' + filename)
        else:
            url = urljoin(BASE_URL, gg_version + '/variances/' + filename)

        file_path = os.path.join(DATA_DIR, filename)

        if os.path.exists(file_path) and not force:
            msg = 'Skipping download of {filename}: file already exists.'
            print(msg.format(filename=filename))
            continue
        else:
            print('Downloading {filename}...'.format(filename=filename))

        try:
            response = urlopen(url)
            with atomic_write(file_path) as f:
                f.write(response.read())
        except HTTPError as e:
            msg = 'There was a problem downloading {filename}: {exception}.'
            print(msg.format(filename=filename, exception=e))
            exit(1)

    print('Done.')


def main():
    _, opts, _ = parse_command_line_parameters(**script_info)
    download_picrust_files(
        output_path=DATA_DIR,
        with_confidence=opts.with_confidence,
        gg_version=opts.gg_version,
        type_of_prediction=opts.type_of_prediction,
        force=opts.force,
    )


if __name__ == '__main__':
    main()
