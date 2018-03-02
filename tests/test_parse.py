#!/usr/bin/env python
# File created on 22 Feb 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "2-alpha.1"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"
 

from cogent.util.unit_test import TestCase, main
from picrust.parse import extract_ids_from_table, parse_asr_confidence_output

class ParseTests(TestCase):
    """ """
    
    def setUp(self):
        """ """
        self.counts_f1 = counts_f1.split('\n')
        self.counts_bad1 = counts_bad1.split('\n')
        self.counts_bad2 = counts_bad2.split('\n')
    
    def test_extract_ids_from_table_qiime_legacy_otu_table(self):
        """extract_ids_from_table extracts ids from a legacy QIIME 1.3 OTU table
        """
        exp = ['39','238','277','1113','1206','1351',\
          '1637','2043','3472','3757','4503','4546',\
          '5403','5479']
        obs = extract_ids_from_table(otu_table_lines)
        self.assertEqual(obs,exp)
        
    
    def test_parse_asr_confidence_output(self):
        """parse_asr_confidence_output parses PICRUST per node confidence output in ML format (includes sigma and loglik values)"""
        
        asr_lines = asr_confidence_output

        min_vals ={\
          'node1': -36.0545,\
          'node10': -15.3024,\
          "'f__Acaryochloridaceae__g__Acaryochloris'":-7.1233,\
          "'f__Synechococcaceae'": -9.2953
        }
        max_vals ={\
          'node1': 44.0219,\
          'node10': 23.2732,\
          "'f__Acaryochloridaceae__g__Acaryochloris'":14.682,\
          "'f__Synechococcaceae'": 17.0651
        }
        sigma = 5986.9627
        loglik = -12353.4949
        
        obs_min_vals,obs_max_vals, params,column_mapping =\
         parse_asr_confidence_output(asr_confidence_output)
        obs_sigma = params['sigma'][0]
        obs_loglik = params['loglik'][0]
        self.assertFloatEqual(obs_sigma,sigma)
        self.assertFloatEqual(obs_loglik, loglik)
        
        #Test all min values
        for k in min_vals.keys():
            self.assertFloatEqual(obs_min_vals[k],min_vals[k])

        #Test all min values
        for k in max_vals.keys():
            self.assertFloatEqual(obs_max_vals[k],max_vals[k])


    def test_parse_asr_confidence_output_pic_format(self):
        """parse_asr_confidence_output parses PICRUST per node confidence output in pic format (*no* sigma or loglik values)"""
        
        asr_lines = asr_confidence_output_pic_format

        min_vals ={\
          'node1': 0.7184,\
          'node10': 1.371,\
          "0.237":1.6493,\
          "'p_Euryarchaeota'": 0.9584
        }
        max_vals ={\
          'node1': 3.7092,\
          'node10': 3.1288,\
          "0.237":3.5431,\
          "'p_Euryarchaeota'": 2.3388
        }
        # No sigma or loglik values, so params should be empty 
        obs_min_vals,obs_max_vals, params,column_mapping =\
         parse_asr_confidence_output(asr_lines)
        #obs_sigma = params['sigma'][0]
        #obs_loglik = params['loglik'][0]
        #self.assertFloatEqual(obs_sigma,sigma)
        #self.assertFloatEqual(obs_loglik, loglik)
        self.assertEqual(params,{})
        
        #Test all min values
        for k in min_vals.keys():
            self.assertFloatEqual(obs_min_vals[k],min_vals[k])

        #Test all min values
        for k in max_vals.keys():
            self.assertFloatEqual(obs_max_vals[k],max_vals[k])



asr_confidence_output=[\
"nodes\t16S_rRNA_Count",\
"node1\t-36.0545|44.0219",\
"node10\t-15.3024|23.2732",\
"'f__Acaryochloridaceae__g__Acaryochloris'\t-7.1233|14.682",\
"'f__Synechococcaceae'\t-9.2953|17.0651",\
"sigma\t5986.9627|NaN",\
"loglik\t-12353.4949"]

asr_confidence_output_pic_format=[\
"nodes\t16S_rRNA_Count",\
"node1\t0.7184|3.7092",\
"node10\t1.371|3.1288",\
"'p_Euryarchaeota'\t0.9584|2.3388",\
"0.237\t1.6493|3.5431"]

otu_table_lines =[\
'# QIIME v1.3.0 OTU table',\
'#OTU ID\tS1\tS2\tS3\tS4\tS5\tS6\tS7\tS8\tS9\tConsensus Lineage',\
'39\t0\t0\t1\t0\t0\t5\t3\t0\t1\tArchaea;Euryarchaeota;Halobacteria;Halobacteriales;Halobacteriaceae',\
 '238\t0\t1\t5\t6\t0\t1\t5\t2\t3\tArchaea;Euryarchaeota',\
 '277\t0\t0\t2\t0\t0\t0\t0\t5\t1\tArchaea;Euryarchaeota',\
 '1113\t0\t0\t0\t0\t1\t5\t0\t0\t0\tArchaea;Euryarchaeota;Halobacteria;Halobacteriales;Halobacteriaceae',\
 '1206\t0\t0\t3\t2\t0\t1\t0\t0\t0\tArchaea;Euryarchaeota',\
 '1351\t0\t0\t0\t0\t0\t0\t1\t11\t5\tArchaea;Euryarchaeota',\
 '1637\t0\t3\t0\t2\t0\t0\t0\t3\t0\tArchaea;Euryarchaeota',\
 '2043\t0\t0\t0\t1\t0\t0\t1\t0\t0\tArchaea;Euryarchaeota;Halobacteria;Halobacteriales;Halobacteriaceae;Halorubrum',\
 '3472\t0\t0\t2\t0\t0\t0\t0\t2\t1\tArchaea;Euryarchaeota',\
 '3757\t0\t0\t1\t2\t0\t3\t0\t0\t0\tArchaea',\
 '4503\t0\t0\t0\t0\t0\t1\t1\t0\t0\tArchaea;Euryarchaeota;Halobacteria;Halobacteriales;Halobacteriaceae',\
 '4546\t0\t3\t0\t3\t0\t0\t9\t6\t2\tArchaea;Euryarchaeota',\
 '5403\t1\t2\t1\t1\t0\t1\t0\t1\t3\tArchaea;Euryarchaeota',\
 '5479\t0\t0\t1\t1\t0\t0\t0\t0\t0\tArchaea']

counts_f1 = """GG_OTU_1	2
GG_OTU_2	4
GG_OTU_3	1"""

counts_bad1 = """GG_OTU_1	hello
GG_OTU_2	4
GG_OTU_3	1"""

counts_bad2 = """GG_OTU_1	42
GG_OTU_2	4
GG_OTU_3	0"""

if __name__ == "__main__":
    main()
