#!/usr/bin/env python
# File created on 22 Feb 2012
from __future__ import division

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
__credits__ = ["Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "2-alpha.1"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"


from cogent.util.unit_test import TestCase, main
from picrust.evaluate_test_datasets import convert_vals_to_spearman_ranks,\
spearman_correlation, calc_spearman_t,evaluate_test_dataset,calculate_accuracy_stats_from_confusion_matrix,\
confusion_matrix_from_data,confusion_matrix_results_by_index,\
roc_points,roc_auc,gini_coefficient,calculate_accuracy_stats_from_observations,\
trapezoidal_area, average_y_values

from random import shuffle
from biom.parse import parse_biom_table

class EvaluateTestDatasetTests(TestCase):
    """ """

    def setUp(self):
       self.genome_table1 = parse_biom_table(genome_table1)
       self.genome_table2 = parse_biom_table(genome_table2)

    def test_evaluate_test_datasets(self):
        """evalutate_test_datasets returns data points and correlations"""

        obs= evaluate_test_dataset(self.genome_table1,self.genome_table2)


    def test_convert_vals_to_spearman_ranks(self):
        """convert_vals_to_spearman_ranks converts a list of floats to their relative ranks."""

        #Example from Spearman Wikipedia page
        #http://en.wikipedia.org/wiki/Spearman's_rank_correlation_coefficient
        #TODO: add some examples from more formal sources


        unordered_vals = [1.2,0.8,1.2,18,2.3]
        exp_ranks = [3.5,5,3.5,1,2]

        obs = convert_vals_to_spearman_ranks(unordered_vals)
        self.assertFloatEqual(obs,exp_ranks)

    def test_spearman_correlation(self):
        """spearman_correlation calculates Spearman rank correlations"""

        #Test data taken from Wikipedia:
        #http://en.wikipedia.org/wiki/Spearman's_rank_correlation_coefficient

        x_vals = [106,86,100,101,99,103,97,113,112,110]
        y_vals = [7,0,27,50,28,29,20,12,6,17]

        r,prob = spearman_correlation(x_vals,y_vals,tails='high')
        exp_rho = -0.175757575
        exp_prob_t_high = 0.6864058
        self.assertFloatEqual(r,exp_rho)
        self.assertFloatEqual(prob,exp_prob_t_high)

    def test_calc_spearman_t(self):
        """calc_spearman_t should produce an adjusted t statistic"""
        r = -0.175757575
        n = 10
        exp_prob_t_high = 0.6864058
        obs_prob_t_high = calc_spearman_t(r,n,'high')

        self.assertFloatEqual(obs_prob_t_high,exp_prob_t_high)

    def test_calculate_accuracy_stats_from_observations(self):
        """calculate_accuracy_stats_from_observations should produce accurate values"""

        exp = [1,1,1,0,0,0,1,1,1,0,1,0]
        fp_only_obs = [1,1,1,1,0,0,1,1,1,1,1,1]




    def test_calculate_accuracy_stats_from_confusion_matrix(self):
        """calc_accuracy_stats_from_confusion_matrix should produce accurate values"""
        #This really only works well for qualitative data, since
        # most values will be somewhat off with quantitative data.
        #The method is generic, however, so either can be used.

        #test data from Wikipedia: http://en.wikipedia.org/wiki/Receiver_operating_characteristic

        #Example A
        tp,fp,fn,tn = 63,28,37,72
        result = calculate_accuracy_stats_from_confusion_matrix(tp,fp,fn,tn)
        self.assertFloatEqual(result['sensitivity'],0.63)
        self.assertFloatEqual(result['false_positive_rate'],0.28)
        self.assertFloatEqual(result['specificity'],0.72)
        self.assertFloatEqual(result['accuracy'],0.675)

        #Example B
        tp,fp,fn,tn = 77,77,23,23
        result = calculate_accuracy_stats_from_confusion_matrix(tp,fp,fn,tn)
        self.assertFloatEqual(result['sensitivity'],0.77)
        self.assertFloatEqual(result['false_positive_rate'],0.77)
        self.assertFloatEqual(result['specificity'],0.23)
        self.assertFloatEqual(result['accuracy'],0.50)

        #Example C
        tp,fp,fn,tn = 24,88,76,12
        result = calculate_accuracy_stats_from_confusion_matrix(tp,fp,fn,tn)
        self.assertFloatEqual(result['sensitivity'],0.24)
        self.assertFloatEqual(result['false_positive_rate'],0.88)
        self.assertFloatEqual(result['specificity'],0.12)
        self.assertFloatEqual(result['accuracy'],0.18)



    def test_confusion_matrix_from_data(self):
        """confusion_matrix_from_data should work for binary and quantitative data"""

        #binary test
        exp = [1,1,1,0,0,0,1,1,1,0,1,0,1]
        obs = [1,1,0,1,0,0,1,1,1,1,1,1,0]
        exp_tp = 6.0
        exp_fp = 3.0
        exp_fn = 2.0
        exp_tn = 2.0

        tp,fp,fn,tn = confusion_matrix_from_data(obs,exp)
        self.assertEqual([tp,fp,fn,tn],\
          [exp_tp,exp_fp,exp_fn,exp_tn])

        #quantitative test
        #results should be indentical to binary test
        #since only presence/absence is considered
        exp = [1,1,1,0,0,0,1,1,1,0,1,0,1]
        obs = [13.7,6.5,0,1,0,0,2.3,1,1.0,1.3,1.5,1,0]
        exp_tp = 6.0
        exp_fp = 3.0
        exp_fn = 2.0
        exp_tn = 2.0

        tp,fp,fn,tn = confusion_matrix_from_data(obs,exp)
        self.assertEqual([tp,fp,fn,tn],\
          [exp_tp,exp_fp,exp_fn,exp_tn])



    def test_confusion_matrix_results_by_index(self):
        """confusion_matrix_results_by_index should return indices for true positives(tp), fp,fn, and tn."""
        #binary test
        exp = [1,1,1,0,0,0,1,1,1,0,1,0,1]
        obs = [1,1,0,1,0,0,1,1,1,1,1,1,0]

        exp_tp = [0,1,6,7,8,10]
        exp_fp = [3,9,11]
        exp_fn = [2,12]
        exp_tn = [4,5]

        tp,fp,fn,tn = confusion_matrix_results_by_index(obs,exp)
        self.assertEqual([tp,fp,fn,tn],\
          [exp_tp,exp_fp,exp_fn,exp_tn])

        #quantitative test
        #results should be indentical to binary test
        #since only presence/absence is considered
        exp = [1,1,1,0,0,0,1,1,1,0,1,0,1]
        obs = [13.7,6.5,0,1,0,0,2.3,1,1.0,1.3,1.5,1,0]

        exp_tp = [0,1,6,7,8,10]
        exp_fp = [3,9,11]
        exp_fn = [2,12]
        exp_tn = [4,5]

        tp,fp,fn,tn = confusion_matrix_results_by_index(obs,exp)
        self.assertEqual([tp,fp,fn,tn],\
          [exp_tp,exp_fp,exp_fn,exp_tn])

    def test_roc_points(self):
        """roc_points should calculate the points for a Receiver Operating Characteristics curve
        """

        #The set up here is a bit elaborate since I generate the test datasets
        #based on the values we need in the confusion matrix.
        #I test the intermediate results though, so any errors should be due
        #to the actual function, not the test

        tn_obs = 0
        tn_exp = 0

        fp_obs = 1
        fp_exp = 0

        tp_obs = 1
        tp_exp = 1

        fn_obs = 0
        fn_exp = 1

         #point A
        obs = [tp_obs] * 63 + [fp_obs] *28 + [fn_obs] * 37 + [tn_obs]*72
        exp = [tp_exp] * 63 + [fp_exp] *28 + [fn_exp] * 37 + [tn_exp]*72
        trial_a_results = confusion_matrix_from_data(obs,exp)
        #Check that this is correct
        self.assertEqual(trial_a_results,(63,28,37,72))
        trial_a = (obs,exp)


        #point B
        obs = [tp_obs] * 77 + [fp_obs] *77 + [fn_obs] * 23 + [tn_obs]*23
        exp = [tp_exp] * 77 + [fp_exp] *77 + [fn_exp] * 23 + [tn_exp]*23
        trial_b_results = confusion_matrix_from_data(obs,exp)
        #Check that this is correct
        self.assertEqual(trial_b_results,(77,77,23,23))
        trial_b = (obs,exp)

        #point c
        obs = [tp_obs] * 24 + [fp_obs] *88 + [fn_obs] * 76 + [tn_obs]*12
        exp = [tp_exp] * 24 + [fp_exp] *88 + [fn_exp] * 76 + [tn_exp]*12
        trial_c_results = confusion_matrix_from_data(obs,exp)
        #Check that this is correct
        self.assertEqual(trial_c_results,(24,88,76,12))
        trial_c_results = calculate_accuracy_stats_from_observations(obs,exp)
        #Check that this is correct
        self.assertFloatEqual(trial_c_results["false_positive_rate"],0.88)

        trial_c = (obs,exp)

        trials = [trial_a, trial_b,trial_c]


        #Finally the actual test

        obs_points = roc_points(trials)
        exp_points = [(0.28,0.63),(0.77,0.77),(0.88,0.24)]
        self.assertFloatEqual(obs_points,exp_points)




    def test_roc_auc(self):
        """roc_auc should calculate the trapezoidal approximation to the area under a ROC curve
        """

        #Uninformative single-point case:
        points = [(0.5,0.5)]
        exp = 0.50
        self.assertFloatEqual(roc_auc(points,add_endpoints=True),exp)

        #Perfect single-point case:
        points = [(0.0,1.0)]
        exp = 1.00
        self.assertFloatEqual(roc_auc(points,add_endpoints=True),exp)

        #Basic test of two-point perfect prediction
        points = [(0.0,1.00),(0.75,1.00)]
        #exp = 0.75*0.75+(0.25*0.75*0.25)*2  #Calculate using one rectangle and two triangles
        exp = 1.0
        self.assertFloatEqual(roc_auc(points,add_endpoints=True),exp)


        # Per Hand and Till 2001, AUC = 2*Gini - 1

        #So I derive the expected value from the Gini example below:

        # From Catalano et al 2009, table 3, page 11
        # Available here:
        #http://scholarcommons.usf.edu/cgi/viewcontent.cgi?article=1032&context=numeracy


        proportion_of_population = [0.1*i for i in range (11)]
        cumulative_portion_of_consumption =\
          [0.000,\
           0.023,\
           0.060,\
           0.110,\
           0.175,\
           0.254,\
           0.345,\
           0.459,\
           0.588,\
           0.754,\
           1.000]
        points = zip(proportion_of_population,\
          cumulative_portion_of_consumption)

        gini_obs = gini_coefficient(points)
        gini_exp = 0.346
        self.assertFloatEqual(gini_obs,gini_exp,eps=1e-3)

        roc_obs = roc_auc(points, add_endpoints=True)
        exp = gini_exp
        #Convert ROC to equivalent Gini index
        obs = abs(2.0*roc_obs - 1.0)
        self.assertFloatEqual(obs,exp,eps=1e-3)


    def test_trapezoidal_area(self):
        """trapezoidal_area should calcualte the area of a trapezoid, given X and Y coords"""

        #simple square
        x1 = 0.0
        x2 = 2.0
        y1 = 2.0
        y2 = 2.0

        exp = 4.0
        obs = trapezoidal_area(x1,y1,x2,y2)
        self.assertFloatEqual(obs,exp)

        #triangle
        x1 = 1.0
        x2 = 3.0
        y1 = 0.0
        y2 = 2.0

        exp = 2.0
        obs = trapezoidal_area(x1,y1,x2,y2)
        self.assertFloatEqual(obs,exp)
    def test_average_y_values(self):
        """average y values should return points for a function by averaging y values for points with same x coordinate"""

        points = [[0,0],[1,3],[1,5],[3,3],[4,1]]
        obs = average_y_values(points)
        exp = [(0,0.0),(1,4.0),(3,3.0),(4,1.0)]
        self.assertEqual(obs,exp)


    def test_gini_coefficient(self):
        """gini_coefficient should calculate the trapezoidal approximation to the Gini coefficient"""
        # From Catalano et al 2009, table 3, page 11
        # Available here:
        #http://scholarcommons.usf.edu/cgi/viewcontent.cgi?article=1032&context=numeracy
        proportion_of_population = [0.1*i for i in range (11)]
        cumulative_portion_of_consumption =\
          [0.000,\
           0.023,\
           0.060,\
           0.110,\
           0.175,\
           0.254,\
           0.345,\
           0.459,\
           0.588,\
           0.754,\
           1.000]
        points = zip(proportion_of_population,\
          cumulative_portion_of_consumption)

        obs = gini_coefficient(points)
        exp = 0.346
        self.assertFloatEqual(obs,exp,eps=1e-3)

genome_table1 = """{"rows": [{"id": "f1", "metadata": null}, {"id": "f2", "metadata": null}, {"id": "f3", "metadata": null}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 1.0], [0, 1, 2.0], [0, 2, 3.0], [1, 1, 1.0], [2, 2, 1.0]], "columns": [{"id": "GG_OTU_1", "metadata": null}, {"id": "GG_OTU_3", "metadata": null}, {"id": "GG_OTU_2", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2753", "matrix_type": "sparse", "shape": [3, 3], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2012-02-22T20:49:58.258296", "type": "OTU table", "id": null, "matrix_element_type": "float"}"""

genome_table2 = """{"rows": [{"id": "f1", "metadata": null}, {"id": "f2", "metadata": null}, {"id": "f3", "metadata": null}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 1.0], [0, 1, 3.0], [0, 2, 4.0], [1, 1, 1.0], [2, 2, 2.0]], "columns": [{"id": "GG_OTU_21", "metadata": null}, {"id": "GG_OTU_23", "metadata": null}, {"id": "GG_OTU_22", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2753", "matrix_type": "sparse", "shape": [3, 3], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2012-02-22T20:49:58.258296", "type": "OTU table", "id": null, "matrix_element_type": "float"}"""

if __name__ == "__main__":
    main()
