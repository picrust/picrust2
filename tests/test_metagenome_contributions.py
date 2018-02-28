#!/usr/bin/env python
# File created on 22 Feb 2012
from __future__ import division

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2012, The PICRUST project"
__credits__ = ["Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "1.1.3"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"


from cogent.util.unit_test import TestCase, main
from biom.parse import parse_biom_table
from picrust.metagenome_contributions import partition_metagenome_contributions,\
  make_pathway_filter_fn
from picrust.predict_metagenomes import predict_metagenomes

class PartitionMetagenomeTests(TestCase):
    """ """

    def setUp(self):
        self.otu_table1 = parse_biom_table(otu_table1)
        self.otu_table_with_taxonomy = parse_biom_table(otu_table_with_taxonomy)
        self.genome_table1 = parse_biom_table(genome_table1)
        self.genome_table2 = parse_biom_table(genome_table2)
        self.predicted_metagenome_table1 = parse_biom_table(predicted_metagenome_table1)
        self.predicted_gene_partition_table = predicted_gene_partition_table
        self.predicted_gene_partition_table_with_taxonomy =\
          predicted_gene_partition_table_with_taxonomy

        #Examples of BIOM format value,id,metadata tuples
        #as returned when iterating over a table
        #metadata are defined at the bottom of this file.
        self.metadata_example = [(700.0,"Gene1",example_metadata1),\
          (250.0,"Gene2",example_metadata2),(0.0,"Gene3",example_metadata3)]


    def test_partition_metagenome_contributions_with_taxonomy(self):
        obs = partition_metagenome_contributions(self.otu_table_with_taxonomy,self.genome_table1)

        obs_text = "\n".join(["\t".join(map(str,i)) for i in obs])
        exp_text_list = [map(str,r.split()) for r in self.predicted_gene_partition_table_with_taxonomy.split('\n')]

        #BIOM adds spaces to metadata fields (not sure why), so add them here just for the taxonomy fields
        for row in exp_text_list[1:]:
            row[9]=' '+row[9]
            row[10]=' '+row[10]
            row[11]=' '+row[11]
            row[12]=' '+row[12]
            row[13]=' '+row[13]
            row[14]=' '+row[14]

        exp_text="\n".join(["\t".join(i) for i in exp_text_list])

        self.assertEqual(obs_text,exp_text)

    def test_make_pathway_filter_fn_KEGG(self):
        """make_pathway_filter_function functions with valid KEGG data"""
        filter_fn = make_pathway_filter_fn(["Phagosome"],\
         metadata_key = 'KEGG_Pathways')

        test_md = self.metadata_example
        exp_result = [True,True,False]
        obs_result = [filter_fn(*md) for md in test_md]
        self.assertEqual(obs_result,exp_result)

    def test_make_pathway_filter_fn_by_level(self):
        """make_pathway_filter_function functions specifying KEGG pathway cat level"""
        test_md = self.metadata_example
        #Example 1 has a top-level 'Unclassifed' annotation
        #Example 3 has a bottom-level 'Unclassified' annotation
        #Test 1.  When not specifying level, all levels should be
        #tested and examples 1 & 3 should both return True

        filter_fn = make_pathway_filter_fn(["Unclassified"],\
         metadata_key = 'KEGG_Pathways')
        exp_result = [True,False,True]
        obs_result = [filter_fn(*md) for md in test_md]
        self.assertEqual(obs_result,exp_result)

        #Test 2. When specifying level 1, only Example 1 is true
        filter_fn = make_pathway_filter_fn(["Unclassified"],\
         metadata_key = 'KEGG_Pathways',search_only_pathway_level = 1)
        exp_result = [True,False,False]
        obs_result = [filter_fn(*md) for md in test_md]
        self.assertEqual(obs_result,exp_result)

        #Test 3. When specifying level 3, only Example 3 is true
        filter_fn = make_pathway_filter_fn(["Unclassified"],\
         metadata_key = 'KEGG_Pathways',search_only_pathway_level = 3)
        exp_result = [False,False,True]
        obs_result = [filter_fn(*md) for md in test_md]
        self.assertEqual(obs_result,exp_result)

    def test_make_pathway_filter_fn_raises_on_bad_input(self):
        """make_pathway_filter_fn raises on invalid metadata level"""
        #No error on specifying level 1
        filter_fn = make_pathway_filter_fn(["Unclassified"],\
            metadata_key = 'KEGG_Pathways',search_only_pathway_level = 1)
        with self.assertRaises(ValueError):
            #Specifying level 0 should raise an informative error
            filter_fn = make_pathway_filter_fn(["Unclassified"],\
              metadata_key = 'KEGG_Pathways',search_only_pathway_level = 0)

    def test_pathway_filter_fn_raises_on_bad_md_key(self):
        """make_pathway_filter_fn child fn raises informatively on invalid metadata"""

        test_md = self.metadata_example

        #Assume user misspells 'KEGG_Pathways'
        filter_fn = make_pathway_filter_fn(["Unclassified"],\
          metadata_key = 'KEGG_Pathwayssssssss')

        with self.assertRaises(ValueError):
            obs_result = [filter_fn(*md) for md in test_md]

    def test_partition_metagenome_contributions(self):
        """partition_metagenome_contributions functions with valid input"""
        #For reference, the OTU table should look like this:
        ##OTU ID Sample1 Sample2 Sample3 Sample4
        #GG_OTU_1    1.0 2.0 3.0 5.0
        #GG_OTU_2    5.0 1.0 0.0 2.0
        #GG_OTU_3    0.0 0.0 1.0 4.0

        #...and the genome table will look like this:
        ##OTU ID GG_OTU_1    GG_OTU_3    GG_OTU_2
        #f1  1.0 2.0 3.0
        #f2  0.0 1.0 0.0
        #f3  0.0 0.0 1.0

        #For which predict metagenomes should produce a table like this:
        ##OTU ID    Sample1 Sample2 Sample3 Sample4
        #f1  16.0    5.0 5.0 19.0
        #f2  0.0 0.0 1.0 4.0
        #f3  5.0 1.0 0.0 2.0

        #First, sanity checks

        #We expect to see the contributions broken down by OTU
        metagenome_table = predict_metagenomes(self.otu_table1,self.genome_table1)
        obs = partition_metagenome_contributions(self.otu_table1,self.genome_table1)

        obs_text = "\n".join(["\t".join(map(str,i)) for i in obs])
        exp_text = "\n".join(["\t".join(map(str,r.split())) for r in \
          self.predicted_gene_partition_table.split('\n')])

        #Test that the percent of all samples is always smaller than
        #the percent of the current sample
        for l in obs[1:]:
            self.assertTrue(l[-1]<=l[-2])

        #Test that the summed contributions equal the metagenome table value
        sum_f1_sample1 = sum([i[5] for i in obs[1:] if (i[0]=="f1" and i[1]=="Sample1")])
        self.assertFloatEqual(sum_f1_sample1,16.0)

        sum_f2_sample1 = sum(\
          [i[5] for i in obs[1:] if (i[0]=="f2" and i[1]=="Sample1")])
        self.assertFloatEqual(sum_f2_sample1,0.0)

        sum_f3_sample1 = sum(\
          [i[5] for i in obs[1:] if (i[0]=="f3" and i[1]=="Sample1")])
        self.assertFloatEqual(sum_f3_sample1,5.0)

        for l in obs[1:]:
            gene,sample,OTU,gene_count_per_genome,otu_abundance_in_sample,count,percent,percent_all = l
            #Test that genomes without genes don't contribute
            #Only GG_OTU_3 has f2, so for all others the gene
            #contribution should be 0,0
            if gene == "f2" and OTU != "GG_OTU_3":
                self.assertFloatEqual(count,0.0)
                self.assertFloatEqual(percent,0.0)
                self.assertFloatEqual(percent_all,0.0)
            #Ditto for GG_OTU_2 and f3
            if gene == "f3" and OTU != "GG_OTU_2":
                self.assertFloatEqual(count,0.0)
                self.assertFloatEqual(percent,0.0)
                self.assertFloatEqual(percent_all,0.0)

            #Test that OTUs absent from a sample don't contribute
            if sample == "Sample1" and OTU == "GG_OTU_3":
                self.assertFloatEqual(count,0.0)
                self.assertFloatEqual(percent,0.0)
                self.assertFloatEqual(percent_all,0.0)

        #Having validated that this looks OK, just compare to
        #hand-checked result
        self.assertEqual(obs_text,exp_text)

        #Check if "limit to functions" works and retrieves the correct information
        obs_limited = partition_metagenome_contributions(self.otu_table1,self.genome_table1,limit_to_functions=["f2"])
        for l in obs_limited[1:]:
            gene,sample,OTU,gene_count_per_genome,otu_abundance_in_sample,count,percent,percent_all = l
            self.assertEqual(gene,"f2")

otu_table1 = """{"rows": [{"id": "GG_OTU_1", "metadata": null}, {"id": "GG_OTU_2", "metadata": null}, {"id": "GG_OTU_3", "metadata": null}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 1.0], [0, 1, 2.0], [0, 2, 3.0], [0, 3, 5.0], [1, 0, 5.0], [1, 1, 1.0], [1, 3, 2.0], [2, 2, 1.0], [2, 3, 4.0]], "columns": [{"id": "Sample1", "metadata": null}, {"id": "Sample2", "metadata": null}, {"id": "Sample3", "metadata": null}, {"id": "Sample4", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2753", "matrix_type": "sparse", "shape": [3, 4], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2012-02-22T20:50:05.024661", "type": "OTU table", "id": null, "matrix_element_type": "float"}"""

otu_table_with_taxonomy = """{"rows": [{"id": "GG_OTU_1", "metadata": {"taxonomy": ["k__1", " p__", " c__", " o__", " f__", " g__", " s__"]}}, {"id": "GG_OTU_2", "metadata": {"taxonomy": ["k__2", " p__", " c__", " o__", " f__", " g__", " s__"]}}, {"id": "GG_OTU_3", "metadata": {"taxonomy": ["k__3", " p__", " c__", " o__", " f__", " g__", " s__"]}}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 1.0], [0, 1, 2.0], [0, 2, 3.0], [0, 3, 5.0], [1, 0, 5.0], [1, 1, 1.0], [1, 3, 2.0], [2, 2, 1.0], [2, 3, 4.0]], "columns": [{"id": "Sample1", "metadata": null}, {"id": "Sample2", "metadata": null}, {"id": "Sample3", "metadata": null}, {"id": "Sample4", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2753", "matrix_type": "sparse", "shape": [3, 4], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2012-02-22T20:50:05.024661", "type": "OTU table", "id": null, "matrix_element_type": "float"}"""

genome_table1 = """{"rows": [{"id": "f1", "metadata": null}, {"id": "f2", "metadata": null}, {"id": "f3", "metadata": null}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 1.0], [0, 1, 2.0], [0, 2, 3.0], [1, 1, 1.0], [2, 2, 1.0]], "columns": [{"id": "GG_OTU_1", "metadata": null}, {"id": "GG_OTU_3", "metadata": null}, {"id": "GG_OTU_2", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2753", "matrix_type": "sparse", "shape": [3, 3], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2012-02-22T20:49:58.258296", "type": "OTU table", "id": null, "matrix_element_type": "float"}"""

genome_table2 = """{"rows": [{"id": "f1", "metadata": null}, {"id": "f2", "metadata": null}, {"id": "f3", "metadata": null}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 1.0], [0, 1, 2.0], [0, 2, 3.0], [1, 1, 1.0], [2, 2, 1.0]], "columns": [{"id": "GG_OTU_21", "metadata": null}, {"id": "GG_OTU_23", "metadata": null}, {"id": "GG_OTU_22", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2753", "matrix_type": "sparse", "shape": [3, 3], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2012-02-22T20:49:58.258296", "type": "OTU table", "id": null, "matrix_element_type": "float"}"""

predicted_metagenome_table1 = """{"rows": [{"id": "f1", "metadata": null}, {"id": "f2", "metadata": null}, {"id": "f3", "metadata": null}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 16.0], [0, 1, 5.0], [0, 2, 5.0], [0, 3, 19.0], [1, 2, 1.0], [1, 3, 4.0], [2, 0, 5.0], [2, 1, 1.0], [2, 3, 2.0]], "columns": [{"id": "Sample1", "metadata": null}, {"id": "Sample2", "metadata": null}, {"id": "Sample3", "metadata": null}, {"id": "Sample4", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2753", "matrix_type": "sparse", "shape": [3, 4], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2012-02-22T16:01:30.837052", "type": "OTU table", "id": null, "matrix_element_type": "float"}"""

predicted_gene_partition_table =\
 """Gene    Sample  OTU GeneCountPerGenome  OTUAbundanceInSample    CountContributedByOTU   ContributionPercentOfSample ContributionPercentOfAllSamples    Kingdom    Phylum    Class    Order    Family    Genus    Species
f1  Sample1 GG_OTU_1    1.0 1.0 1.0 0.0625  0.0222222222222
f1  Sample1 GG_OTU_2    3.0 5.0 15.0    0.9375  0.333333333333
f1  Sample2 GG_OTU_1    1.0 2.0 2.0 0.4 0.0444444444444
f1  Sample2 GG_OTU_2    3.0 1.0 3.0 0.6 0.0666666666667
f1  Sample3 GG_OTU_1    1.0 3.0 3.0 0.6 0.0666666666667
f1  Sample3 GG_OTU_3    2.0 1.0 2.0 0.4 0.0444444444444
f1  Sample4 GG_OTU_1    1.0 5.0 5.0 0.263157894737  0.111111111111
f1  Sample4 GG_OTU_2    3.0 2.0 6.0 0.315789473684  0.133333333333
f1  Sample4 GG_OTU_3    2.0 4.0 8.0 0.421052631579  0.177777777778
f2  Sample3 GG_OTU_3    1.0 1.0 1.0 1.0 0.2
f2  Sample4 GG_OTU_3    1.0 4.0 4.0 1.0 0.8
f3  Sample1 GG_OTU_2    1.0 5.0 5.0 1.0 0.625
f3  Sample2 GG_OTU_2    1.0 1.0 1.0 1.0 0.125
f3  Sample4 GG_OTU_2    1.0 2.0 2.0 1.0 0.25"""

predicted_gene_partition_table_with_taxonomy =\
    """Gene    Sample  OTU GeneCountPerGenome  OTUAbundanceInSample    CountContributedByOTU   ContributionPercentOfSample ContributionPercentOfAllSamples    Kingdom    Phylum    Class    Order    Family    Genus    Species
f1  Sample1 GG_OTU_1    1.0 1.0 1.0 0.0625 0.0222222222222 k__1 p__ c__ o__ f__ g__ s__
f1  Sample1 GG_OTU_2    3.0 5.0 15.0    0.9375  0.333333333333    k__2      p__     c__     o__     f__     g__     s__
f1  Sample2 GG_OTU_1    1.0 2.0 2.0 0.4 0.0444444444444    k__1      p__     c__     o__     f__     g__     s__
f1  Sample2 GG_OTU_2    3.0 1.0 3.0 0.6 0.0666666666667    k__2      p__     c__     o__     f__     g__     s__
f1  Sample3 GG_OTU_1    1.0 3.0 3.0 0.6 0.0666666666667    k__1      p__     c__     o__     f__     g__     s__
f1  Sample3 GG_OTU_3    2.0 1.0 2.0 0.4 0.0444444444444    k__3      p__     c__     o__     f__     g__     s__
f1  Sample4 GG_OTU_1    1.0 5.0 5.0 0.263157894737  0.111111111111    k__1      p__     c__     o__     f__     g__     s__
f1  Sample4 GG_OTU_2    3.0 2.0 6.0 0.315789473684  0.133333333333    k__2      p__     c__     o__     f__     g__     s__
f1  Sample4 GG_OTU_3    2.0 4.0 8.0 0.421052631579  0.177777777778    k__3      p__     c__     o__     f__     g__     s__
f2  Sample3 GG_OTU_3    1.0 1.0 1.0 1.0 0.2    k__3      p__     c__     o__     f__     g__     s__
f2  Sample4 GG_OTU_3    1.0 4.0 4.0 1.0 0.8    k__3      p__     c__     o__     f__     g__     s__
f3  Sample1 GG_OTU_2    1.0 5.0 5.0 1.0 0.625    k__2      p__     c__     o__     f__     g__     s__
f3  Sample2 GG_OTU_2    1.0 1.0 1.0 1.0 0.125    k__2      p__     c__     o__     f__     g__     s__
f3  Sample4 GG_OTU_2    1.0 2.0 2.0 1.0 0.25    k__2      p__     c__     o__     f__     g__     s__"""


#Simple gene with only 1 pathway annotation
example_metadata1 = {"KEGG_Description": ["cathepsin L [EC:3.4.22.15]"],\
 "KEGG_Pathways": [["Cellular Processes", "Transport and Catabolism",\
    "Phagosome"],["Unclassified"]]}
#Complex annoation with phagosome
example_metadata2 = {"KEGG_Description": ["cathepsin L [EC:3.4.22.15]"],\
  "KEGG_Pathways":\
    [["Human Diseases", "Immune System Diseases","Rheumatoid arthritis"],\
     ["Metabolism", "Enzyme Families", "Peptidases"],
     ["Organismal Systems", "Immune System",\
       "Antigen processing and presentation"],\
     ["Cellular Processes", "Transport and Catabolism", "Phagosome"],\
     ["Genetic Information Processing", "Folding, Sorting and Degradation",\
       "Chaperones and folding catalysts"],\
     ["Cellular Processes", "Transport and Catabolism", "Lysosome"]]}
#Complex annotation without phagosome
example_metadata3 = {"KEGG_Description": ["cathepsin L [EC:3.4.22.15]"],\
  "KEGG_Pathways":\
    [["Human Diseases", "Immune System Diseases","Rheumatoid arthritis"],\
     ["Metabolism", "Enzyme Families", "Peptidases"],\
     ["Organismal Systems", "Immune System",\
       "Antigen processing and presentation"],\
     ["Cellular Processes", "Transport and Catabolism","Unclassified"],\
     ["Genetic Information Processing", "Folding, Sorting and Degradation",\
       "Chaperones and folding catalysts"],\
     ["Cellular Processes", "Transport and Catabolism", "Lysosome"]]}


if __name__ == "__main__":
    main()
