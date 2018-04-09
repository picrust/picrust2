#!/usr/bin/env python

from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2018, The PICRUSt Project"
__credits__ = ["Greg Caporaso","Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "2-alpha.6"

from numpy import array
from cogent.util.unit_test import TestCase, main
from biom.parse import parse_biom_table, get_axis_indices,\
  direct_slice_data
from picrust2.predict_metagenomes import predict_metagenomes,\
  calc_nsti,get_overlapping_ids,\
  extract_otu_and_genome_data,transfer_sample_metadata,\
  transfer_observation_metadata,transfer_metadata,\
  load_subset_from_biom_str,yield_subset_biom_str,\
  predict_metagenome_variances,variance_of_sum,variance_of_product,\
  sum_rows_with_variance

class PredictMetagenomeTests(TestCase):
    """ """

    def setUp(self):
        #Datasets for metagenome prediction
        self.otu_table1 = parse_biom_table(otu_table1)
        self.otu_table1_with_metadata = parse_biom_table(otu_table1_with_metadata)
        self.genome_table1 = parse_biom_table(genome_table1)
        self.genome_table1_with_metadata = parse_biom_table(genome_table1_with_metadata)
        self.genome_table2 = parse_biom_table(genome_table2)
        self.predicted_metagenome_table1 = parse_biom_table(predicted_metagenome_table1)
        self.predicted_metagenome_table1_with_metadata = parse_biom_table(predicted_metagenome_table1_with_metadata)

        #Datasets for variance estimation during metagenome prediction
        self.zero_variance_table1 = parse_biom_table(zero_variance_table1)
        self.variance_table1_var_by_otu = parse_biom_table(variance_table1_var_by_otu)
        self.variance_table1_var_by_gene = parse_biom_table(variance_table1_var_by_gene)
        self.variance_table1_one_gene_one_otu = parse_biom_table(variance_table1_one_gene_one_otu)

        self.predicted_metagenome_table1_zero_variance = parse_biom_table(predicted_metagenome_table1_zero_variance)
        self.predicted_metagenome_variance_table1_one_gene_one_otu =\
            parse_biom_table(predicted_metagenome_variance_table1_one_gene_one_otu)
        self.predicted_metagenome_table1_one_gene_one = parse_biom_table(predicted_metagenome_table1)

        #Datasets for testing confidence intervals
        self.predicted_metagenome_table1_one_gene_one_otu_upper_CI =\
          parse_biom_table(predicted_metagenome_table1_one_gene_one_otu_upper_CI)
        self.predicted_metagenome_table1_one_gene_one_otu_lower_CI =\
          parse_biom_table(predicted_metagenome_table1_one_gene_one_otu_lower_CI)

    def test_predict_metagenomes(self):
        """ predict_metagenomes functions as expected with valid input """
        actual = predict_metagenomes(self.otu_table1,self.genome_table1)
        self.assertEqual(str(actual),
                         str(self.predicted_metagenome_table1))

    def test_predict_metagenomes_value_error(self):
        """ predict_metagenomes raises ValueError when no overlapping otu ids """
        self.assertRaises(ValueError,predict_metagenomes,self.otu_table1,self.genome_table2)

    def test_predict_metagenome_variances_returns_zero_variance_from_zero_variance(self):
        """ predict_metagenomes outputs correct results given zero variance input"""

        curr_otu_table = self.otu_table1
        curr_genome_table = self.genome_table1
        curr_variance_table = self.zero_variance_table1
        curr_exp_metagenome_table = self.predicted_metagenome_table1
        curr_exp_metagenome_variance_table = self.predicted_metagenome_table1_zero_variance

        obs_prediction,obs_variances,obs_lower_CI_95,obs_upper_CI_95 =\
          predict_metagenome_variances(curr_otu_table,curr_genome_table,gene_variances=curr_variance_table)

        #Test that the prediction itself is as expected
        self.assertEqual(str(obs_prediction),
                         str(curr_exp_metagenome_table))

        #Test that the variance prediction is all zeros, as expected
        self.assertEqual(str(obs_variances),
                         str(curr_exp_metagenome_variance_table))

        #Test that with zero variance, the upper and lower CIs are equal to the expected value (i.e. the prediction)
        self.assertEqual(str(obs_lower_CI_95),
                         str(curr_exp_metagenome_table))
        self.assertEqual(str(obs_upper_CI_95),
                         str(curr_exp_metagenome_table))

    def test_predict_metagenome_variances_propagates_variance_in_gene_categories(self):
        """ predict_metagenomes correctly propagates the rank order of gene family variance"""
        curr_otu_table = self.otu_table1
        curr_genome_table = self.genome_table1
        curr_variance_table = self.variance_table1_var_by_gene
        curr_exp_metagenome_table = self.predicted_metagenome_table1
        obs_prediction,obs_variances,obs_lower_CI_95,obs_upper_CI_95 =\
          predict_metagenome_variances(curr_otu_table,curr_genome_table,gene_variances=curr_variance_table)

        #Check that the metagenome prediction hasn't changed
        self.assertEqual(str(obs_prediction),
                         str(curr_exp_metagenome_table))

    def test_predict_metagenome_variances_propagates_variance(self):
        """ predict_metagenomes correctly propagates differences in gene family variance as expected in a simple example"""

        curr_otu_table = self.otu_table1
        curr_genome_table = self.genome_table1
        curr_variance_table = self.variance_table1_one_gene_one_otu
        curr_exp_metagenome_table = self.predicted_metagenome_table1
        curr_exp_metagenome_varaiance_table = self.predicted_metagenome_variance_table1_one_gene_one_otu
        curr_exp_upper_CI_95 = self.predicted_metagenome_table1_one_gene_one_otu_upper_CI
        curr_exp_lower_CI_95 = self.predicted_metagenome_table1_one_gene_one_otu_lower_CI

        obs_prediction,obs_variances,obs_lower_CI_95,obs_upper_CI_95 =\
          predict_metagenome_variances(curr_otu_table,curr_genome_table,gene_variances=curr_variance_table)

        self.assertEqual(str(obs_prediction),
                         str(curr_exp_metagenome_table))
        #Expect no variance in f1 or f2 in any sample, and no variance in OTU 1 or 3.
        #Otu 2 occurs in all samples except sample 3, so all samples except 3 should
        #have variance.   The exact values follow from variance of scaled random variables or
        #The sum of random variables
        self.assertEqual(obs_variances,self.predicted_metagenome_variance_table1_one_gene_one_otu)

        #Check CIs against hand calculated CIs
        self.assertEqual(str(obs_upper_CI_95),
                         str(curr_exp_upper_CI_95))
        self.assertEqual(str(obs_lower_CI_95),
                         str(curr_exp_lower_CI_95))

    def test_predict_metagenomes_keeps_observation_metadata(self):
        """predict_metagenomes preserves Observation metadata in genome and otu table"""

        actual = predict_metagenomes(self.otu_table1_with_metadata,self.genome_table1_with_metadata)
        exp = self.predicted_metagenome_table1_with_metadata

        #NOTE: the expected data is  mapped to dicts below because otherwise the memory
        #location of the lambda function associated with the defaultdict
        #causes (artifactual) inequality of results

        actual_md = map(dict,sorted([md for md in actual.metadata(axis='observation')]))
        exp_md = map(dict,sorted([md for md in exp.metadata(axis='observation')]))
        for i,md in enumerate(actual_md):
            self.assertEqualItems(md,exp_md[i])
        for i,md in enumerate(exp_md):
            self.assertEqualItems(md,actual_md[i])

    def test_predict_metagenomes_keeps_sample_metadata(self):
        """predict_metagenomes preserves Sample metadata in genome and otu table"""
        #NOTE: could be consolidated with "_keeps_observation_metadata above

        actual = predict_metagenomes(self.otu_table1_with_metadata,\
          self.genome_table1_with_metadata,verbose=False)
        exp = self.predicted_metagenome_table1_with_metadata

        #Need to map to dicts, otherwise the memory location of the lambda function
        #associated with the defaultdict causes (artifactual) inequality of results

        actual_md = map(dict,sorted([md for md in actual.metadata()]))
        exp_md = map(dict,sorted([md for md in exp.metadata()]))
        for i,md in enumerate(actual_md):
            self.assertEqualItems(md,exp_md[i])
        for i,md in enumerate(exp_md):
            self.assertEqualItems(md,actual_md[i])

    def test_transfer_metadata_moves_sample_metadata_between_biom_tables(self):
        """transfer_metadata moves sample metadata values between BIOM format tables"""
        t1 = self.otu_table1
        exp = self.otu_table1_with_metadata
        actual = transfer_metadata(self.otu_table1_with_metadata,self.otu_table1,\
          "sample","sample",verbose=False)

        actual_md = map(dict,sorted([md for md in actual.metadata()]))
        exp_md = map(dict,sorted([md for md in exp.metadata()]))
        for i,md in enumerate(actual_md):
            self.assertEqualItems(md,exp_md[i])
        for i,md in enumerate(exp_md):
            self.assertEqualItems(md,actual_md[i])

    def test_transfer_metadata_moves_observation_metadata_between_biom_tables(self):
        """transfer_metadata moves observation metadata values between BIOM format tables"""
        t1 = self.genome_table1
        exp = self.genome_table1_with_metadata
        actual = transfer_metadata(self.genome_table1_with_metadata,\
          self.genome_table1,"observation","observation",verbose=False)

        actual_md = map(dict,sorted([md for md in actual.metadata(axis='observation')]))
        exp_md = map(dict,sorted([md for md in exp.metadata(axis='observation')]))
        for i,md in enumerate(actual_md):
            self.assertEqualItems(md,exp_md[i])
        for i,md in enumerate(exp_md):
            self.assertEqualItems(md,actual_md[i])


    def test_transfer_sample_metadata_moves_sample_metadata_between_biom_tables(self):
        """transfer_sample_metadata moves sample metadata values between BIOM format tables"""
        t1 = self.otu_table1
        exp = self.otu_table1_with_metadata
        actual = transfer_sample_metadata(self.otu_table1_with_metadata,\
          self.otu_table1,"sample",verbose=False)
        actual_md = map(dict,sorted([md for md in actual.metadata()]))
        exp_md = map(dict,sorted([md for md in exp.metadata()]))
        for i,md in enumerate(actual_md):
            self.assertEqualItems(md,exp_md[i])
        for i,md in enumerate(exp_md):
            self.assertEqualItems(md,actual_md[i])

    def test_transfer_observation_metadata_moves_observation_metadata_between_biom_tables(self):
        """transfer_sample_metadata moves sample metadata values between BIOM format tables"""
        t1 = self.genome_table1
        exp = self.genome_table1_with_metadata
        actual = transfer_observation_metadata(self.genome_table1_with_metadata,\
          self.genome_table1,"observation",verbose=False)

        actual_md = map(dict,sorted([md for md in actual.metadata(axis='observation')]))
        exp_md = map(dict,sorted([md for md in exp.metadata(axis='observation')]))
        for i,md in enumerate(actual_md):
            self.assertEqualItems(md,exp_md[i])
        for i,md in enumerate(exp_md):
            self.assertEqualItems(md,actual_md[i])

    def test_load_subset_from_biom_str_loads_a_subset_of_observations(self):
        """load_subset_from_biom_str loads a subset of observations from a valid BIOM format JSON string"""
        biom_str = otu_table1_with_metadata
        ids_to_load = ['GG_OTU_1','GG_OTU_2']
        axis = 'observations'
        #NOTE: this will fail currently due to a known bug in the BIOM direct_parse_key
        #as soon as that is updated this should pass, however
        #two_taxon_table = load_subset_from_biom_str(biom_str,ids_to_load,axis)
        #self.assertEqualItems(two_taxon_table.ObservationIds,ids_to_load)

        #Test that loading all ids is identical to just loading the table
        #exp = parse_biom(biom_str)
        #obs = load_subset_from_biom_str(biom_str,ids_to_load=['GG_OTU_1','GG_OTU_2','GG_OTU_3'],axis=axis)
        #self.assertEqual(obs,exp)

    def test_variance_of_sum_functions_as_expected_with_valid_input(self):
        """variance_of_sum functions as expected given two variances"""
        #Example drawn from:http://onlinestatbook.com/chapter4/variance_sum_law2.html
        var1 = 10000
        var2 = 11000
        r=0.50
        # expected variance = 10,000 + 11,000 + 2(0.5)*sqrt(10,000)*sqrt(11,000)=31488
        expected_var1 = 31488.088481701518
        observed_var = variance_of_sum(var1,var2,r,sign_of_varB=1)
        self.assertFloatEqual(expected_var1,observed_var)
        expected_var2 = 10511.911518298484
        observed_var = variance_of_sum(var1,var2,r,sign_of_varB=-1)
        self.assertFloatEqual(expected_var2,observed_var)
        #Test that this works for vector input
        var1=array([10000,10000,0])
        var2=array([11000,11000,0])
        observed_var = variance_of_sum(var1,var2,r)
        expected_var = array([expected_var1,expected_var1,0.0])
        self.assertFloatEqual(observed_var,expected_var)

    def test_sum_rows_with_variance(self):
        """sum_rows_with_variance sums the rows of a numpy array while accounting for variance"""
        data_array = array([[0,0],[0,1.0]])
        variance_array = array([[1.0,0],[0,1000.0]])
        exp_data_array = array([0.0,1.0])
        exp_variance_array = array([1.0,1000.0])
        obs_data_array,obs_variance_array =\
          sum_rows_with_variance(data_array,variance_array)
        self.assertFloatEqual(obs_data_array,exp_data_array)
        self.assertFloatEqual(obs_variance_array,exp_variance_array)

    def test_variance_of_product_functions_as_expected_with_valid_input(self):
        """variance_of_product functions as expected given two values and two variances"""
        varA = 100.0
        A = 1.0
        varB = 1000.0
        B = 10.0
        r=0.5
        #Expected (calc'd by hand) = (100/1)**2 + (1000/10)**2 + 2*0.5*sqrt(100)*sqrt(1000)/10.0
        #Equivalently = 10000 + 10000 + sqrt(1000) = 20000 + sqrt(1000) = 20000 + 31.622776601683793
        expected = 20031.622776601683793
        observed = variance_of_product(A,B,varA,varB,r=0.5)
        self.assertFloatEqual(observed,expected)
        #Test taat this works for vector input
        Av = array([A]*10)
        Bv = array([B]*10)
        varAv=array([varA]*10)
        varBv=array([varB]*10)
        rv = array([r]*10)
        result = variance_of_product(Av,Bv,varAv,varBv,rv)
        self.assertFloatEqual(result,array([expected]*10))

    def test_yield_subset_biom_str_yields_string_pieces_from_valid_input(self):
        """yield_subset_biom_str yields components of a biom string containing only a subset of ids, given a valid biom str"""
        biom_str = otu_table1_with_metadata
        ids_to_load = ['GG_OTU_1','GG_OTU_2']
        axis = 'observation'

        idxs, new_axis_md = get_axis_indices(biom_str,ids_to_load, axis)
        new_data = direct_slice_data(biom_str,idxs, axis)

        #NOTE: this will fail currently due to a known bug in the BIOM direct_parse_key
        #as soon as that is updated this should pass, however
        obs = [part for part in yield_subset_biom_str(biom_str,new_data,new_axis_md,axis)]
        exp = ['{', '"id": "GG_OTU_1"', ',',\
                '"format": "Biological Observation Matrix v0.9"', ',',\
                '"format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html"', ',',\
                '"type": "Gene table"', ',',\
#                '"generated_by": "QIIME 1.4.0-dev, svn revision 2753', ',',\
                '"generated_by": "QIIME 1.4.0-dev', ',',\
                '"date": "2012-02-22T20:50:05.024661"', ',',\
                '"matrix_type": "sparse"', ',',\
                '"matrix_element_type": "float"', ',',\
                '"data": [[0,0,1.0],[0,1,2.0],[0,2,3.0],[0,3,5.0],[1,0,5.0],[1,1,1.0],[1,3,2.0]], "shape": [2, 4]',',',\
                '"rows": [{"id": "GG_OTU_1", "metadata": null}, {"id": "GG_OTU_2", "metadata": null}]', ',',\
                '"columns": [{"id": "Sample1", "metadata": {"pH":7.0}}, {"id": "Sample2", "metadata": {"pH":8.0}}, {"id": "Sample3", "metadata": {"pH":7.0}}, {"id": "Sample4", "metadata": null}]', '}']
        #For now be aware that commas in generated_by
        #strings won't parse correnctly
        for i,piece in enumerate(exp):
            self.assertEqual(obs[i],piece)
        #TODO: when BIOM direct_parse_key is fixed this should pass
        #self.assertEqual(obs,exp)



otu_table1 = """{"rows": [{"id": "GG_OTU_1", "metadata": null}, {"id": "GG_OTU_2", "metadata": null}, {"id": "GG_OTU_3", "metadata": null}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 1.0], [0, 1, 2.0], [0, 2, 3.0], [0, 3, 5.0], [1, 0, 5.0], [1, 1, 1.0], [1, 3, 2.0], [2, 2, 1.0], [2, 3, 4.0]], "columns": [{"id": "Sample1", "metadata": null}, {"id": "Sample2", "metadata": null}, {"id": "Sample3", "metadata": null}, {"id": "Sample4", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2753", "matrix_type": "sparse", "shape": [3, 4], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2012-02-22T20:50:05.024661", "type": "Gene table", "id": null, "matrix_element_type": "float"}"""

otu_table1_with_metadata = """{"rows": [{"id": "GG_OTU_1", "metadata": null}, {"id": "GG_OTU_2", "metadata": null}, {"id": "GG_OTU_3", "metadata": null}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 1.0], [0, 1, 2.0], [0, 2, 3.0], [0, 3, 5.0], [1, 0, 5.0], [1, 1, 1.0], [1, 3, 2.0], [2, 2, 1.0], [2, 3, 4.0]], "columns": [{"id": "Sample1", "metadata": {"pH":7.0}}, {"id": "Sample2", "metadata": {"pH":8.0}}, {"id": "Sample3", "metadata": {"pH":7.0}}, {"id": "Sample4", "metadata": null}],"generated_by": "QIIME 1.4.0-dev, svn revision 2753", "matrix_type": "sparse", "shape": [3, 4], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2012-02-22T20:50:05.024661", "type": "Gene table", "id": null, "matrix_element_type": "float"}"""

genome_table1 = """{"rows": [{"id": "f1", "metadata": null}, {"id": "f2", "metadata": null}, {"id": "f3", "metadata": null}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 1.0], [0, 1, 2.0], [0, 2, 3.0], [1, 1, 1.0], [2, 2, 1.0]], "columns": [{"id": "GG_OTU_1", "metadata": null}, {"id": "GG_OTU_3", "metadata": null}, {"id": "GG_OTU_2", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2753", "matrix_type": "sparse", "shape": [3, 3], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2012-02-22T20:49:58.258296", "type": "Gene table", "id": null, "matrix_element_type": "float"}"""

zero_variance_table1  = """{"rows": [{"id": "f1", "metadata": null}, {"id": "f2", "metadata": null}, {"id": "f3", "metadata": null}], "format":
"Biological Observation Matrix v0.9", "data": [[0, 0, 0.0], [0, 1, 0.0], [0, 2, 0.0], [1, 1, 0.0], [2, 2, 0.0]], "columns": [{"id": "GG_OTU_1", "metadata": null}, {"id": "GG_OTU_3", "metadata": null}, {"id": "GG_OTU_2", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2753", "matrix_type": "sparse", "shape": [3, 3], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2012-02-22T20:49:58.258296", "type": "Gene table", "id": null, "matrix_element_type": "float"}"""

variance_table1 = """{"rows": [{"id": "f1", "metadata": null}, {"id": "f2", "metadata": null}, {"id": "f3", "metadata": null}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 0.0], [0, 1, 0.0], [0, 2, 1.0], [1, 1, 10.0], [2, 2, 100.0]], "columns": [{"id": "GG_OTU_1", "metadata": null}, {"id": "GG_OTU_3", "metadata": null}, {"id": "GG_OTU_2", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2753", "matrix_type": "sparse", "shape": [3, 3], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2012-02-22T20:49:58.258296", "type": "Gene table", "id": null, "matrix_element_type": "float"}"""

variance_table1_var_by_gene = """{"rows": [{"id": "f1", "metadata": null}, {"id": "f2", "metadata": null}, {"id": "f3", "metadata": null}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 0.0], [0, 1, 0.0], [0, 2, 0.0], [1, 1, 1.0],[2, 1, 10.0],[1, 2, 1.0], [2, 2, 10.0]], "columns": [{"id": "GG_OTU_1", "metadata": null}, {"id": "GG_OTU_3", "metadata": null}, {"id": "GG_OTU_2", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2753", "matrix_type": "sparse", "shape": [3, 3], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2012-02-22T20:49:58.258296", "type": "Gene table", "id": null, "matrix_element_type": "float"}"""

variance_table1_var_by_otu = """{"rows": [{"id": "f1", "metadata": null}, {"id": "f2", "metadata": null}, {"id": "f3", "metadata": null}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 0.0], [0, 1, 10.0], [0, 2, 100.0], [1, 1, 10.0],[2, 1, 10.0],[1, 2, 100.0], [2, 2, 100.0]], "columns": [{"id": "GG_OTU_1", "metadata": null}, {"id": "GG_OTU_3", "metadata": null}, {"id": "GG_OTU_2", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2753", "matrix_type": "sparse", "shape": [3, 3], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2012-02-22T20:49:58.258296", "type": "Gene table", "id": null, "matrix_element_type": "float"}"""

variance_table1_one_gene_one_otu = """{"rows": [{"id": "f1", "metadata": null}, {"id": "f2", "metadata": null}, {"id": "f3", "metadata": null}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 0.0], [0, 1, 0.0], [0, 2, 0.0], [1, 1, 0.0],[2, 1, 0.0],[1, 2, 0.0], [2, 2, 1000.0]], "columns": [{"id": "GG_OTU_1", "metadata": null}, {"id": "GG_OTU_3", "metadata": null}, {"id": "GG_OTU_2", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2753", "matrix_type": "sparse", "shape": [3, 3], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2012-02-22T20:49:58.258296", "type": "Gene table", "id": null, "matrix_element_type": "float"}"""

genome_table1_with_metadata = """{"rows": [{"id": "f1", "metadata": {"KEGG_description":"ko00100    Steroid biosynthesis"}}, {"id": "f2", "metadata": {"KEGG_description":"ko00195   Photosynthesis"}}, {"id": "f3", "metadata": {"KEGG_description":"ko00232    Caffeine metabolism"}}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 1.0], [0, 1, 2.0], [0, 2, 3.0], [1, 1, 1.0], [2, 2, 1.0]], "columns": [{"id": "GG_OTU_1", "metadata": {"confidence": 0.665,"taxonomy": ["Root", "k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Lachnospiraceae"]}}, {"id": "GG_OTU_3", "metadata": {"confidence": 1.0,"taxonomy": ["Root", "k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Lachnospiraceae"]}}, {"id": "GG_OTU_2", "metadata":{"confidence": 0.98,"taxonomy": ["Root", "k__Bacteria"]}}], "generated_by": "QIIME 1.4.0-dev, svn revision 2753", "matrix_type": "sparse", "shape": [3, 3], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2012-02-22T20:49:58.258296", "type": "Gene table", "id": null, "matrix_element_type": "float"}"""

genome_table2 = """{"rows": [{"id": "f1", "metadata": null}, {"id": "f2", "metadata": null}, {"id": "f3", "metadata": null}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 1.0], [0, 1, 2.0], [0, 2, 3.0], [1, 1, 1.0], [2, 2, 1.0]], "columns": [{"id": "GG_OTU_21", "metadata": null}, {"id": "GG_OTU_23", "metadata": null}, {"id": "GG_OTU_22", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2753", "matrix_type": "sparse", "shape": [3, 3], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2012-02-22T20:49:58.258296", "type": "Gene table", "id": null, "matrix_element_type": "float"}"""

predicted_metagenome_table1 = """{"rows": [{"id": "f1", "metadata": null}, {"id": "f2", "metadata": null}, {"id": "f3", "metadata": null}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 16.0], [0, 1, 5.0], [0, 2, 5.0], [0, 3, 19.0], [1, 2, 1.0], [1, 3, 4.0], [2, 0, 5.0], [2, 1, 1.0], [2, 3, 2.0]], "columns": [{"id": "Sample1", "metadata": null}, {"id": "Sample2", "metadata": null}, {"id": "Sample3", "metadata": null}, {"id": "Sample4", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2753", "matrix_type": "sparse", "shape": [3, 4], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2012-02-22T16:01:30.837052", "type": "Gene table", "id": null, "matrix_element_type": "float"}"""

predicted_metagenome_table1_one_gene_one_otu_upper_CI = """{"rows": [{"id": "f1", "metadata": null}, {"id": "f2", "metadata": null}, {"id": "f3", "metadata": null}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 16.0], [0, 1, 5.0], [0, 2, 5.0], [0, 3, 19.0], [1, 2, 1.0], [1, 3, 4.0], [2, 0, 315.0], [2, 1, 63.0], [2, 3, 126.0]], "columns": [{"id": "Sample1", "metadata": null}, {"id": "Sample2", "metadata": null}, {"id": "Sample3", "metadata": null}, {"id": "Sample4", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2753", "matrix_type": "sparse", "shape": [3, 4], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2012-02-22T16:01:30.837052", "type": "Gene table", "id": null, "matrix_element_type": "float"}"""

predicted_metagenome_table1_one_gene_one_otu_lower_CI = """{"rows": [{"id": "f1", "metadata": null}, {"id": "f2", "metadata": null}, {"id": "f3", "metadata": null}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 16.0], [0, 1, 5.0], [0, 2, 5.0], [0, 3, 19.0], [1, 2, 1.0], [1, 3, 4.0], [2, 0, 0.0], [2, 1, 0.0], [2, 3, 0.0]], "columns": [{"id": "Sample1", "metadata": null}, {"id": "Sample2", "metadata": null}, {"id": "Sample3", "metadata": null}, {"id": "Sample4", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2753", "matrix_type": "sparse", "shape": [3, 4], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2012-02-22T16:01:30.837052", "type": "Gene table", "id": null, "matrix_element_type": "float"}"""

predicted_metagenome_variance_table1_one_gene_one_otu = """{"rows": [{"id": "f1", "metadata": null}, {"id": "f2", "metadata": null}, {"id": "f3", "metadata": null}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 0.0], [0, 1, 0.0], [0, 2, 0.0], [0, 3, 0.0], [1, 2, 0.0], [1, 3, 0.0], [2, 0, 25000.0], [2, 1, 1000.0], [2, 3, 4000.0]], "columns": [{"id": "Sample1", "metadata": null}, {"id": "Sample2", "metadata": null}, {"id": "Sample3", "metadata": null}, {"id": "Sample4", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2753", "matrix_type": "sparse", "shape": [3, 4], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2012-02-22T16:01:30.837052", "type": "Gene table", "id": null, "matrix_element_type": "float"}"""


predicted_metagenome_table1_zero_variance = """{"rows": [{"id": "f1", "metadata": null}, {"id": "f2", "metadata": null}, {"id": "f3", "metadata": null}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 0.0], [0, 1, 0.0], [0, 2, 0.0], [0, 3, 0.0], [1, 2, 0.0], [1, 3, 0.0], [2, 0, 0.0], [2, 1, 0.0], [2, 3, 0.0]], "columns": [{"id": "Sample1", "metadata": null}, {"id": "Sample2", "metadata": null}, {"id": "Sample3", "metadata": null}, {"id": "Sample4", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2753", "matrix_type": "sparse", "shape": [3, 4], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2012-02-22T16:01:30.837052", "type": "Gene table", "id": null, "matrix_element_type": "float"}"""

predicted_metagenome_table1_with_metadata = """{"rows": [{"id": "f1", "metadata": {"KEGG_description":"ko00100    Steroid biosynthesis"}}, {"id": "f2", "metadata": {"KEGG_description":"ko00195   Photosynthesis"}}, {"id": "f3", "metadata": {"KEGG_description":"ko00232    Caffeine metabolism"}}], "format": "Biological Observation Matrix v0.9","data": [[0, 0, 16.0], [0, 1, 5.0], [0, 2, 5.0], [0, 3, 19.0], [1, 2, 1.0], [1, 3, 4.0], [2, 0, 5.0], [2, 1, 1.0], [2, 3, 2.0]], "columns": [{"id": "Sample1", "metadata": {"pH":7.0}}, {"id": "Sample2", "metadata": {"pH":8.0}}, {"id": "Sample3", "metadata": {"pH":7.0}}, {"id": "Sample4", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2753", "matrix_type": "sparse", "shape": [3, 4], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2012-02-22T16:01:30.837052", "type": "Gene table", "id": null, "matrix_element_type": "float"}"""

if __name__ == "__main__":
    main()
