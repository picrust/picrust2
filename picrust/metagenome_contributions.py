#!/usr/bin/env python
# File created on 22 Feb 2012
from __future__ import division

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2015, The PICRUST project"
__credits__ = ["Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "1.1.3"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"


from numpy import dot, array, around
from biom.table import Table
from picrust.predict_metagenomes import get_overlapping_ids,extract_otu_and_genome_data


def make_pathway_filter_fn(ok_values,metadata_key='KEGG_Pathways',\
  search_only_pathway_level=None):
    """return a filter function that filters observations by pathway

    ok_values -- valid pathway category names (e.g. ['Metabolism'] in the KEGG
      pathway categories.

      By default, if the pathway is present at any level, keep the observation.
      Matches can be limited to a specific level of the pathway
      using search_only_pathway_level

    metadata_key -- the metadata key to filter observations against.
      e.g. 'KEGG_Pathways','COG_Category'

    search_only_pathway_level -- only check a given level of the pathway
    hierarchy.  NOTE: KEGG labels pathway levels from 1 not 0, and this
    function matches KEGG.  So specifying 2 will search KEGG level 2 or
    or python index 1 of the metadata.

    NOTES:
    metadata for metadata_key is expected to be a list of annotations,
    with each annotation a list of pathway levels

    COG Example:
    "metadata": {"COG_Description": "Uncharacterized conserved protein",\
      "COG_Category": [["POORLY CHARACTERIZED", "[S] Function unknown"]]}}

    KEGG Example:
    "metadata": {"KEGG_Description": ["cathepsin L [EC:3.4.22.15]"],\
      "KEGG_Pathways":\
      [["Human Diseases", "Immune System Diseases","Rheumatoid arthritis"],\
      ["Metabolism", "Enzyme Families", "Peptidases"],
      ["Organismal Systems", "Immune System",\
      "Antigen processing and presentation"],\
      ["Cellular Processes", "Transport and Catabolism", "Phagosome"],\
      ["Genetic Information Processing", "Folding, Sorting and Degradation",\
     "Chaperones and folding catalysts"],\
      ["Cellular Processes", "Transport and Catabolism", "Lysosome"]]}
    """
    pathway_index = None
    if search_only_pathway_level is not None:
        pathway_index = search_only_pathway_level - 1
        if pathway_index < 0:
            raise ValueError(\
              "make_pathway_filter_fn accepts only positive integer values of search_only_pathway_level")

    def filter_observation_by_pathway(obs_value,obs_id,obs_metadata):
        function_md = obs_metadata.get(metadata_key,None)
        if not function_md:
            #empty metadata
            err_str =\
              "Empty observation metadata.\
              Is the metadata key '%s' valid? Valid keys are: %s"\
              %(str(metadata_key),str(obs_metadata.keys()))
            raise ValueError(err_str)
        #Note that many pathway annotations per observation are allowed
        for annotation in function_md:
            for level,function in enumerate(annotation):
                #Skip this level if it isn't the desired level
                if pathway_index is not None:
                    if level != pathway_index:
                        continue
                if function in ok_values:
                    return True

        #If we've scanned all values and no match is found,
        #return False to discard the observation
        return False

    return filter_observation_by_pathway

def partition_metagenome_contributions(otu_table,genome_table, limit_to_functions=[],
        limit_to_functional_categories=[], metadata_key = 'KEGG_Pathways',remove_zero_rows=True,verbose=True):
    """Return a list of the contribution of each organism to each function, per sample
    (rewritten version using numpy)
    otu_table -- the BIOM Table object for the OTU table
    genome_table -- the BIOM Table object for the predicted genomes
    limit_to_functions -- a list of function ids to include.
      If empty, include all function ids

    limit_by_function_categories -- if provided limit by functional category.
      For example, this can be used to limit output by KEGG functional categories

    Output table as a list of lists with header
    Function\tOrganism\tSample\tCounts\tpercent_of_sample
    """

    if limit_to_functions:
        if verbose:
            print "Filtering the genome table to include only user-specified functions:",limit_to_functions
        ok_ids = frozenset(map(str,limit_to_functions))

        filter_by_set = lambda vals,gene_id,metadata: str(gene_id) in ok_ids
        genome_table = genome_table.filter(filter_by_set, axis='observation')

        if genome_table.is_empty():
            raise ValueError("User filtering by functions (%s) removed all results from the genome table"%(str(limit_to_functions)))

    if limit_to_functional_categories:
        fn_cat_filter = make_pathway_filter_fn(ok_values = frozenset(map(str,limit_to_functional_categories)),metadata_key=metadata_key)
        genome_table = genome_table.filter(fn_cat_filter, axis='observation', inplace=False)

        if genome_table.is_empty():
            raise ValueError("User filtering by functional categories (%s) removed all results from the genome table"%(str(limit_to_functional_categories)))

    otu_data,genome_data,overlapping_ids = extract_otu_and_genome_data(otu_table,genome_table)
    #We have a list of data with abundances and gene copy numbers
    lines=[]
    result = [["Gene","Sample","OTU","GeneCountPerGenome",\
            "OTUAbundanceInSample","CountContributedByOTU",\
            "ContributionPercentOfSample","ContributionPercentOfAllSamples",
                   "Kingdom","Phylum","Class","Order","Family","Genus","Species"]]


    #Zero-valued total counts will be set to epsilon
    epsilon = 1e-5

    for j,gene_id in enumerate(genome_table.ids(axis='observation')):
        all_gene_rows = []
        for k,sample_id in enumerate(otu_table.ids()):
            #Add raw counts for the gene in this sample to a list
            sample_gene_rows = []
            for i,otu_id in enumerate(overlapping_ids):
                otu_gene_count = genome_data[i][j]
                otu_abundance = otu_data[i][k]
                contribution =  otu_gene_count * otu_abundance
                if remove_zero_rows and contribution == 0.0:
                    #skip zero contributions
                    continue
                sample_gene_rows.append([gene_id,sample_id,otu_id,otu_gene_count,otu_abundance,contribution])
            #Now get the percentage of each genes contribution to the sample overall
            total_counts =max(epsilon,sum([float(row[-1]) for row in sample_gene_rows]))

            for row in sample_gene_rows:
                percent_of_sample = float(row[-1])/total_counts
                row.append(percent_of_sample)
            all_gene_rows.extend(sample_gene_rows)

        count_idx = -2 #Counts are now in the next to last position in each row
        total_counts =max(epsilon,sum([float(row[count_idx]) for row in all_gene_rows]))
        otu_index=2 #position of otu ids in the table

        o_md = otu_table.metadata(axis='observation')
        for row in all_gene_rows:
            percent_of_sample = float(row[count_idx])/total_counts
            row.append(percent_of_sample)

            #add taxonomy information for each OTU
            obs_index = otu_table.index(row[otu_index], 'observation')
            if o_md is not None and 'taxonomy' in o_md[obs_index]:
                row.extend(o_md[obs_index]['taxonomy'])

        lines.extend(all_gene_rows)
    result.extend(lines)

    return result
