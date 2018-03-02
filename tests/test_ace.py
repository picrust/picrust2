#!/usr/bin/env python
# File created on 27 Feb 2012
from __future__ import division

__author__ = "Morgan Langille"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
__credits__ = ["Morgan Langille"]
__license__ = "GPL"
__version__ = "2-alpha.1"
__maintainer__ = "Morgan Langille"
__email__ = "morgan.g.i.langille@gmail.com"
__status__ = "Development"
 

from cogent.util.unit_test import TestCase, main
from picrust.wrap_asr import ape_ace_wrapper
from cogent.app.util import get_tmp_filename
from cogent.util.misc import remove_files
from cogent import LoadTable
from cogent.util.table import Table

class AceTests(TestCase):
    """ """
    
    def setUp(self):
        """ """
        #create a tmp tree file
        self.in_tree1_fp = get_tmp_filename(prefix='AceTests',suffix='.nwk')
        self.in_tree1_file = open(self.in_tree1_fp,'w')
        self.in_tree1_file.write(in_tree1)
        self.in_tree1_file.close()

        #create a tmp tree file (with underscores in tip names)
        self.in_tree2_fp = get_tmp_filename(prefix='AceTests',suffix='.nwk')
        self.in_tree2_file = open(self.in_tree2_fp,'w')
        self.in_tree2_file.write(in_tree2)
        self.in_tree2_file.close()

        #create a tmp trait file
        self.in_trait1_fp = get_tmp_filename(prefix='AceTests',suffix='.tsv')
        self.in_trait1_file=open(self.in_trait1_fp,'w')
        self.in_trait1_file.write(in_trait1)
        self.in_trait1_file.close()

        #create another tmp trait file (need to test table with only single column seperately)
        self.in_trait2_fp = get_tmp_filename(prefix='AceTests',suffix='.tsv')
        self.in_trait2_file=open(self.in_trait2_fp,'w')
        self.in_trait2_file.write(in_trait2)
        self.in_trait2_file.close()

        #create a tmp trait file (with underscores in tip names)
        self.in_trait3_fp = get_tmp_filename(prefix='AceTests',suffix='.tsv')
        self.in_trait3_file=open(self.in_trait3_fp,'w')
        self.in_trait3_file.write(in_trait3)
        self.in_trait3_file.close()

        self.files_to_remove = [self.in_tree1_fp,self.in_trait1_fp,self.in_trait2_fp, self.in_trait3_fp, self.in_tree2_fp]

    def tearDown(self):
        remove_files(self.files_to_remove)
                       

    def test_ape_ace_wrapper_ml(self):
        """ test_ape_ace_wrapper with method 'ML' functions as expected with valid input
        """
        actual,actual_ci= ape_ace_wrapper(self.in_tree1_fp, self.in_trait1_fp, method="ML")
        expected=Table(['nodes','trait1','trait2'],[['14','2.9737','2.5436'],['12','2.3701','2.7056'],['11','0.8370','2.9706'],['10','4.4826','2.1388']])
        self.assertEqual(actual.tostring(),expected.tostring())
        expected_ci=Table(['nodes','trait1','trait2'],\
                              [['14','1.4467|4.5007','2.1979|2.8894'],\
                               ['12','0.9729|3.7674','2.3892|3.0219'],\
                               ['11','0.147|1.527','2.8143|3.1268'],\
                               ['10','3.4227|5.5426','1.8988|2.3788'],\
                               ['sigma','1.9742|0.6981','0.1012|0.0359'],\
                               ['loglik','-6.7207','5.1623'],\
                               ])
        self.assertEqual(actual_ci.tostring(),expected_ci.tostring())

   

    def test_ape_ace_wrapper_pic(self):
        """ test_ape_ace_wrapper with method 'pic' functions as expected with valid input
        """
        actual,actual_ci= ape_ace_wrapper(self.in_tree1_fp,self.in_trait1_fp, method="pic")
        expected=Table(['nodes','trait1','trait2'],[['14','2.9737','2.5436'],['12','1.2727','3'],['11','0.6667','3'],['10','5','2']])
        self.assertEqual(actual.tostring(),expected.tostring())
        expected_ci=Table(['nodes','trait1','trait2'],\
                              [['14','0.7955|5.1519','0.3655|4.7218'],\
                               ['12','-1.1009|3.6464','0.6264|5.3736'],\
                               ['11','-0.4068|1.7402','1.9265|4.0735'],\
                               ['10','3.3602|6.6398','0.3602|3.6398'],\
                               ])
        self.assertEqual(actual_ci.tostring(),expected_ci.tostring())


    def test_ape_ace_wrapper_pic_single_trait(self):
        """ test_ape_ace_wrapper with method 'pic' functions as expected with single column trait table
        """
        actual,ci= ape_ace_wrapper(self.in_tree1_fp,self.in_trait2_fp, method="pic")
        expected=Table(['nodes','trait1'],[['14','2.9737'],['12','1.2727'],['11','0.6667'],['10','5']])
        self.assertEqual(actual.tostring(),expected.tostring())

    def test_ape_ace_wrapper_pic_with_funky_tip_labels(self):
        """ test_ape_ace_wrapper for a tree with underscores in tip labels
        """
        actual,ci= ape_ace_wrapper(self.in_tree2_fp,self.in_trait3_fp, method="pic")
        expected=Table(['nodes','trait1','trait2'],[['14','2.9737','2.5436'],['12','1.2727','3'],['11','0.6667','3'],['10','5','2']])
        self.assertEqual(actual.tostring(),expected.tostring())
        
                                                       
in_tree1="""(((1:0.1,2:0.2)11:0.6,3:0.8)12:0.2,(4:0.3,D:0.4)10:0.5)14;"""

in_tree2="""((('abc_123':0.1,2:0.2)11:0.6,3:0.8)12:0.2,('NC_2345':0.3,D:0.4)10:0.5)14;"""


in_trait1="""tips	trait1	trait2
1	1	3
2	0	3
3	2	3
4	5	2
D	5	2"""

in_trait2="""tips	trait1
1	1
2	0
3	2
4	5
D	5"""

in_trait3="""tips	trait1	trait2
abc_123	1	3
2	0	3
3	2	3
NC_2345	5	2
D	5	2"""


if __name__ == "__main__":
    main()
