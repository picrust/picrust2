#!/usr/bin/env python

from __future__ import division

__author__ = "Morgan Langille"
__copyright__ = "Copyright 2015, The PICRUSt Project"
__credits__ = ["Morgan Langille"]
__license__ = "GPL"
__version__ = "2-alpha.1"
__maintainer__ = "Morgan Langille"
__email__ = "morgan.g.i.langille@gmail.com"
__status__ = "Development"


from os import remove
from cogent import LoadTable
from cogent.util.table import Table
from cogent.app.util import get_tmp_filename
from picrust.util import get_picrust_project_dir
from os.path import join

from picrust.parallel import submit_jobs, wait_for_output_files


def combine_asr_tables(output_files,verbose=False):
    """ Combine all tables coming from asr output. Cuts 2nd column out and joins them together into single table.
    Assumes all output files have same row identifiers and that these are in the same order.
    """

    #Going to store an array of arrays here
    combined_table=[]

    #load in the first column (containing row ids). File doesn't matter since they should all have identical first columns.
    table=LoadTable(filename=output_files[0],header=True,sep='\t')
    row_ids = table.getRawData(columns=[table.Header[0]])
    combined_table.append([table.Header[0]])
    for row_id in row_ids:
        combined_table.append([row_id])

    #Now add the rest of the files to the table
    for i,output_file in enumerate(output_files):
        if verbose:
            print "Combining file {0} of {1}: {2}".format(i,len(output_files),output_file)
        #pull out the second column (first column with actual preditions)
        table=LoadTable(filename=output_file,header=True,sep='\t')
        predictions = table.getRawData(columns=[table.Header[1]])

        #Add the header for our column to the list of headers
        combined_table[0].append(table.Header[1])

        #Add rest of values in the column
        j=1
        for prediction in predictions:
            combined_table[j].append(prediction)
            j+=1

    return combined_table

def run_asr_in_parallel(tree, table, asr_method, parallel_method='sge',tmp_dir='jobs/',num_jobs=100, verbose=False):
    '''Runs the ancestral state reconstructions in parallel'''

    asr_script_fp = join(get_picrust_project_dir(),'scripts','ancestral_state_reconstruction.py')

    if(parallel_method=='sge'):
        cluster_jobs_fp=join(get_picrust_project_dir(),'scripts','start_parallel_picrust_jobs_sge.py')
    elif(parallel_method=='multithreaded'):
        cluster_jobs_fp=join(get_picrust_project_dir(),'scripts','start_parallel_picrust_jobs.py')
    elif(parallel_method=='torque'):
        cluster_jobs_fp=join(get_picrust_project_dir(),'scripts','start_parallel_picrust_jobs_torque.py')
    else:
        raise RuntimeError

    if(verbose):
        print "Loading trait table..."

    #foreach trait in the table, create a new tmp file with just that trait, and create the job command and add it a tmp jobs file
    table=LoadTable(filename=table, header=True, sep='\t')

    #get dimensions of the table
    dim=table.Shape

    created_tmp_files=[]
    output_files=[]
    ci_files=[]

    #create a tmp file to store the job commands (which we will pass to our parallel script to run)
    jobs_fp=get_tmp_filename(tmp_dir=tmp_dir,prefix='jobs_asr_')
    jobs=open(jobs_fp,'w')
    created_tmp_files.append(jobs_fp)

    if(verbose):
        print "Creating temporary input files in: ",tmp_dir

    #iterate over each column
    for i in range(1,dim[1]):
        #create a new table with only a single trait
        single_col_table=table.getColumns([0,i])

        #write the new table to a tmp file
        single_col_fp=get_tmp_filename(tmp_dir=tmp_dir,prefix='in_asr_')
        single_col_table.writeToFile(single_col_fp,sep='\t')
        created_tmp_files.append(single_col_fp)

        #create tmp output files
        tmp_output_fp=get_tmp_filename(tmp_dir=tmp_dir,prefix='out_asr_')
        output_files.append(tmp_output_fp)
        tmp_ci_fp=get_tmp_filename(tmp_dir=tmp_dir,prefix='out_asr_ci_')
        ci_files.append(tmp_ci_fp)

        #create the job command
        cmd= "{0} -i {1} -t {2} -m {3} -o {4} -c {5}".format(asr_script_fp, single_col_fp, tree, asr_method, tmp_output_fp, tmp_ci_fp)

        #add job command to the the jobs file
        jobs.write(cmd+"\n")

    jobs.close()
    created_tmp_files.extend(output_files)
    created_tmp_files.extend(ci_files)

    if(verbose):
        print "Launching parallel jobs."

    #run the job command
    job_prefix='asr'
    submit_jobs(cluster_jobs_fp ,jobs_fp,job_prefix,num_jobs=num_jobs)

    if(verbose):
        print "Jobs are now running. Will wait until finished."

    #wait until all jobs finished (e.g. simple poller)
    wait_for_output_files(output_files)

    if(verbose):
        print "Jobs are done running. Now combining all tmp files."
    #Combine output files
    combined_table=combine_asr_tables(output_files)
    combined_ci_table=combine_asr_tables(ci_files)

    #create a Table object
    combined_table=Table(header=combined_table[0],rows=combined_table[1:])
    combined_ci_table=Table(header=combined_ci_table[0],rows=combined_ci_table[1:])

    #clean up all tmp files
    for file in created_tmp_files:
        remove(file)

    #return the combined table
    return combined_table,combined_ci_table
