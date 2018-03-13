#!/usr/bin/env python
# File created on 23 Nov 2011
from __future__ import division


__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2015, The PICRUSt Project"
__credits__ = ["Greg Caporaso", "Morgan Langille", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "2-alpha.2"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"


from biom.parse import parse_biom_table
from cogent.parse.tree import DndParser
from cogent.util.unit_test import main, TestCase
import os
from os.path import abspath, dirname, exists
from picrust.util import (
    atomic_write,
    biom_meta_to_string,
    convert_biom_to_precalc,
    convert_precalc_to_biom,
    get_picrust_project_dir,
    PicrustNode,
    transpose_trait_table_fields,
)
import StringIO
from tempfile import mkstemp


class PicrustNodeTests(TestCase):
    def setUp(self):
        pass

    def test_getsubtree(self):
        """testing getting a subtree
        """
        otu_names = ['NineBande', 'Mouse', 'HowlerMon', 'DogFaced']
        newick = '(((Human,HowlerMon),Mouse),NineBande,DogFaced);'
        newick_reduced = '((Mouse,HowlerMon),NineBande,DogFaced);'
        tree = DndParser(newick, constructor = PicrustNode)

        subtree = tree.getSubTree(otu_names)
        new_tree = DndParser(newick_reduced, constructor = PicrustNode)
        # check we get the same names
        self.assertEqual(*[len(t.Children) for t in (subtree,new_tree)])
        self.assertEqual(subtree.getNewick(), new_tree.getNewick())


#    def test_multifurcating(self):
#        """Coerces nodes to have <= n children"""
#        t_str = "((a:1,b:2,c:3)d:4,(e:5,f:6,g:7)h:8,(i:9,j:10,k:11)l:12)m:14;"
#        t = DndParser(t_str)

        # can't break up easily... sorry 80char
#        exp_str = "((a:1.0,(b:2.0,c:3.0):0.0)d:4.0,((e:5.0,(f:6.0,g:7.0):0.0)h:8.0,(i:9.0,(j:10.0,k:11.0):0.0)l:12.0):0.0)m:14.0;"
#        obs = t.multifurcating(2)
#        self.assertEqual(obs.getNewick(with_distances=True), exp_str)
#        self.assertNotEqual(t.getNewick(with_distances=True),
#                            obs.getNewick(with_distances=True))

#        obs = t.multifurcating(2, 0.5)
#        exp_str = "((a:1.0,(b:2.0,c:3.0):0.5)d:4.0,((e:5.0,(f:6.0,g:7.0):0.5)h:8.0,(i:9.0,(j:10.0,k:11.0):0.5)l:12.0):0.5)m:14.0;"
#        self.assertEqual(obs.getNewick(with_distances=True), exp_str)

#        t_str = "((a,b,c)d,(e,f,g)h,(i,j,k)l)m;"
#        exp_str = "((a,(b,c))d,((e,(f,g))h,(i,(j,k))l))m;"
#        t = DndParser(t_str, constructor=PicrustNode)
#        obs = t.multifurcating(2)
#        self.assertEqual(obs.getNewick(with_distances=True), exp_str)
#        obs = t.multifurcating(2, eps=10) # no effect on TreeNode type
#        self.assertEqual(obs.getNewick(with_distances=True), exp_str)

#        self.assertRaises(TreeError, t.multifurcating, 1)

    def test_bifurcating(self):
        """Coerces nodes to have <= 2 children"""
        t_str = "((a:1,b:2,c:3)d:4,(e:5,f:6,g:7)h:8,(i:9,j:10,k:11)l:12)m:14;"
        t = DndParser(t_str)

        # can't break up easily... sorry 80char
        exp_str = "((a:1.0,(b:2.0,c:3.0):0.0)d:4.0,((e:5.0,(f:6.0,g:7.0):0.0)h:8.0,(i:9.0,(j:10.0,k:11.0):0.0)l:12.0):0.0)m:14.0;"
        obs = t.bifurcating()



class UtilTests(TestCase):
    """ Tests of the picrust/util.py module """

    def setUp(self):
        """ Initialize variables: run before each test """
        self.precalc_in_biom = parse_biom_table(precalc_in_biom)


    def tearDown(self):
        """ Clean up: run after each test """
        pass

    def test_convert_precalc_to_biom(self):
        """ convert_precalc_to_biom as expected with valid input """
        # test if table conversion works (and the ablity to handle a string as input)
        result_table = convert_precalc_to_biom(precalc_in_tab)
        self.assertEqual(result_table,self.precalc_in_biom)

        # test partial loading (and the ability to handle a file handle as input)
        ids_to_load = ['OTU_1','OTU_2']
        two_taxon_table = convert_precalc_to_biom(StringIO.StringIO(precalc_in_tab),ids_to_load)
        self.assertEqualItems(two_taxon_table.ids(), ids_to_load)

    def test_convert_precalc_to_biom_value_error(self):
        """ convert_precalc_to_biom raises ValueError when no overlapping otu ids or additional ids """
        self.assertRaises(ValueError,convert_precalc_to_biom,precalc_in_tab,['bogus_id1','bogus_id2'])
        self.assertRaises(ValueError,convert_precalc_to_biom,precalc_in_tab,['OTU_1','bogus_id2'])


    def test_convert_biom_to_precalc(self):
        """ convert_biom_to_precalc as expected with valid input """

        result = convert_biom_to_precalc(parse_biom_table(precalc_in_biom))
        self.assertEqual(result,precalc_in_tab)

    def test_biom_meta_to_string(self):
        """ biom_meta_to_string functions as expected """

        meta_simple="foo;bar"
        expect_simple="foo:bar"
        result_simple=biom_meta_to_string(meta_simple)
        self.assertEqual(result_simple,expect_simple)

        meta_list=['foo;','bar']
        expect_list="foo:;bar"
        result_list=biom_meta_to_string(meta_list)
        self.assertEqual(result_list,expect_list)

        meta_list_of_list=[['foo;','bar'],['hello','world|']]
        expect_list_of_list="foo:;bar|hello;world:"
        result_list_of_list=biom_meta_to_string(meta_list_of_list)
        self.assertEqual(result_list_of_list,expect_list_of_list)


    def test_get_picrust_project_dir(self):
        """get_picrust_project_dir functions as expected"""

        # Do an explicit check on whether the file system containing
        # the current file is case insensitive.
        case_insensitive_filesystem = \
         exists(__file__.upper()) and exists(__file__.lower())

        actual = get_picrust_project_dir()
        # I base the expected here off the imported location of
        # picrust/util.py here, to handle cases where either the user has
        # PICRUST in their PYTHONPATH, or when they've installed it with
        # setup.py.
        # If util.py moves this test will fail -- that
        # is what we want in this case, as the get_picrust_project_dir()
        # function would need to be modified.
        import picrust.util
        util_py_filepath = abspath(abspath(picrust.util.__file__))
        expected = dirname(dirname(util_py_filepath))

        if case_insensitive_filesystem:
            # make both lowercase if the file system is case insensitive
            actual = actual.lower()
            expected = expected.lower()
        self.assertEqual(actual,expected)

    def test_transpose_trait_table_fields(self):
        """transpose_table_fields should correctly transpose table fields"""
        # The purpose is to make it easy to transpose tab-delimited tables
        # before analysis starts (and to have an internal method where needed

        # Get header
        test_header = "genomes\t123456\t123457\t123458\t123459\n"
        test_data_fields =[\
          ["gene1",0,3,5,0],\
          ["gene2",1,2,1,0],\
          ["gene3",0,0,0,1]]

        exp_header  = "genomes\tgene1\tgene2\tgene3\n"

        exp_data_fields =[\
          ['123456',0,1,0],\
          ['123457',3,2,0],\
          ['123458',5,1,0],\
          ['123459',0,0,1]]

        new_header,new_data_fields = transpose_trait_table_fields(test_data_fields,\
          header=test_header,id_row_idx=0)

        self.assertEqual(new_header,exp_header)
        self.assertEqual(new_data_fields,exp_data_fields)

    def test_atomic_write_deletes_on_fail(self):
        _, tempfile_path = mkstemp()
        os.remove(tempfile_path)  # We reuse the path, not the file itself

        try:
            with atomic_write(tempfile_path) as f:
                f.write('fail')
                raise Exception
        except Exception:
            pass

        self.assertFalse(exists(tempfile_path))

    def test_atomic_write_transparent_on_success(self):
        _, tempfile_path = mkstemp()
        os.remove(tempfile_path)  # We reuse the path, not the file itself

        with atomic_write(tempfile_path) as f:
            f.write('success')

        with open(tempfile_path, 'rb') as f:
            self.assertEqual('success', f.read())

        os.remove(tempfile_path)


precalc_in_tab="""#OTU_IDs	f1	f2	f3	metadata_NSTI
metadata_simple	f1_desc	f2_desc	f3_desc
metadata_list	f1;l1;l2	f2;l1;l2	f3;l2;l3
metadata_list_of_lists	f1;l1;l2	f2;l1;l2|f2;l1a;l2a	f3;l3;l2
OTU_1	1.0	2.0	3.0	1.2
OTU_2	0.0	0.0	0.0	2.3
OTU_3	4.0	4.0	4.0	0.5"""

precalc_in_biom="""{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "Gene table","generated_by": "test","date": "2013-07-20T23:44:49.970541","matrix_type": "dense","matrix_element_type": "float","shape": [3, 3],"data": [[1.0,0.0,4.0],[2.0,0.0,4.0],[3.0,0.0,4.0]],"rows": [{"id": "f1", "metadata": {"simple": "f1_desc", "list": ["f1", "l1", "l2"], "list_of_lists": [["f1", "l1", "l2"]]}},{"id": "f2", "metadata": {"simple": "f2_desc", "list": ["f2", "l1", "l2"], "list_of_lists": [["f2", "l1", "l2"], ["f2", "l1a", "l2a"]]}},{"id": "f3", "metadata": {"simple": "f3_desc", "list": ["f3", "l2", "l3"], "list_of_lists": [["f3", "l3", "l2"]]}}],"columns": [{"id": "OTU_1", "metadata": {"NSTI": "1.2"}},{"id": "OTU_2", "metadata": {"NSTI": "2.3"}},{"id": "OTU_3", "metadata": {"NSTI": "0.5"}}]}"""

if __name__ == "__main__":
    main()
