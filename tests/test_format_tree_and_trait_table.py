#!/usr/bin/env python

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
__credits__ = ["Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "1.1.3"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"


from cogent.util.unit_test import main, TestCase
from cogent import LoadTree
from cogent.parse.tree import DndParser
from picrust.format_tree_and_trait_table import reformat_tree_and_trait_table,\
  nexus_lines_from_tree,add_branch_length_to_root,\
  set_min_branch_length,make_nexus_trees_block,\
  filter_table_by_presence_in_tree,convert_trait_table_entries,\
  yield_trait_table_fields,ensure_root_is_bifurcating,\
  filter_tree_tips_by_presence_in_table,print_node_summary_table,\
  add_to_filename,make_id_mapping_dict,make_translate_conversion_fn,\
  make_char_translation_fn,remove_spaces, format_tree_node_names


"""
Tests for format_tree_and_trait_tables.py
"""

class TestFormat(TestCase):
    """Tests of format_tree_and_trait_table.py"""

    def setUp(self):
        self.SimpleTree = \
          DndParser("((A:0.02,B:0.01)E:0.05,(C:0.01,D:0.01)F:0.05)root;")
    
        self.SimpleTreeWithSpaces = \
          DndParser("((E coli:0.02,S typhimurium :0.01)Gamma proteobacteria:0.05,(C\t:0.01,D:0.01)F:0.05)root;")
        
        self.SimplePolytomyTree = \
          DndParser("((A:0.02,B:0.01,B_prime:0.03)E:0.05,(C:0.01,D:0.01)F:0.05)root;")
        
        self.SimpleUnlabelledTree = \
          DndParser("((A:0.02,B:0.01):0.05,(C:0.01,D:0.01):0.05)root;")
        
        #First number is GG id, the second is IMG
        self.GreengenesToIMG = \
          [('469810','645058788'),\
          ('457471','645058789'),\
          ('266998','641736109')]

    def reformat_tree_and_trait_table(self):
        """Test the main function under various conditions"""
        pass
   
    def test_reformat_tree_node_names(self):
        """Reformat_tree_node_names reformats tree nodes using supplied fns"""
        
        tree = self.SimpleTreeWithSpaces
        
        obs_tree = format_tree_node_names(tree,[remove_spaces])
        obs = obs_tree.getNewick(with_distances=True)

        exp =\
          "((E_coli:0.02,S_typhimurium :0.01)Gamma_proteobacteria:0.05,(C:0.01,D:0.01)F:0.05)root;"



    def test_nexus_lines_from_tree(self):
        """Nexus lines from tree should return NEXUS formatted lines..."""
        obs =  nexus_lines_from_tree(self.SimpleTree)
        exp = ['#NEXUS', 'begin trees;', '\ttranslate', '\t\t0 A,', '\t\t1 B,',\
                '\t\t2 C,', '\t\t3 D;',\
                '\t\ttree PyCogent_tree = ((0:0.02,1:0.01):0.05,(2:0.01,3:0.01):0.05)root;', 'end;']
        self.assertEqual(obs,exp)
        
    def test_add_branch_length_to_root(self):
        """Add branch length to root should add epsilon branch lengths"""
        pass

    def test_set_min_branch_length(self):
        """set_min_branch_length should set a minimum branch length"""
        tree = self.SimpleTree
        
        obs = set_min_branch_length(tree,min_length = 0.04)
        exp = "((A:0.04,B:0.04)E:0.05,(C:0.04,D:0.04)F:0.05)root;"
        self.assertEqual(obs.getNewick(with_distances=True),exp)

    def test_make_nexus_trees_block(self):
        """make_nexus_trees_block should output a trees block in NEXUS format"""
        pass

    def test_filter_table_by_presence_in_tree(self):
        """filter_trait_table_by_presence_in_tree should filter trait table"""
        pass

    def test_convert_trait_table_entries(self):
        """convert_trait_entries should convert labels,values using conversion fns"""
        lines =[\
          ['organism 1','0','0.3','15','1','6'],\
          ['organism 2','1','1','13','-1',-1],\
          ['organism 3','2','0','12','0.9','5']]
        
        val_conv_fns = []
        
        #Handle replacement of ambiguous values
        replacement_dict ={'-':0,'-1':0,-1:0,'NULL':0,None:0}
        replace_ambig_fn = make_translate_conversion_fn(replacement_dict)
        val_conv_fns.append(replace_ambig_fn)
        label_conv_fns = [remove_spaces]
       
        #Convert values to integers (in string form)
        val_conv_fns.append(lambda x: str(int(float(x))))
       
        obs = [l for l in convert_trait_table_entries(lines,\
                value_conversion_fns = val_conv_fns,\
                label_conversion_fns=label_conv_fns)]
        exp =[\
          ['organism_1','0','0','15','1','6'],\
          ['organism_2','1','1','13','0','0'],\
          ['organism_3','2','0','12','0','5']]

        self.assertEqual(obs,exp)

    def test_make_translate_conversion_fn(self):
        """make_translate_conversion_fn should make a conversion_fn given dicts """
        
        translation_dict =\
          {'-1':0,-1:0,'Unknown':0,'?':0}

        conversion_fn =\
          make_translate_conversion_fn(translation_dict)
        
        values = [0,-1,15,'10','-1','Unknown','?']
        exp = ['0','0','15','10','0','0','0']
        obs = map(conversion_fn,values)
        self.assertEqual(obs,exp)

    def test_make_char_translation_conversion_fn(self):
        """make_translate_conversion_fn should make a conversion_fn given dicts """
        #Replace semicolons and colons with underscores
        translation_dict =\
          {';':'_',':':'_'}

        conversion_fn =\
          make_char_translation_fn(translation_dict)
        
        values = [24687,23458,'Proteobacteria;Gammaproteobacteria',\
                  'test;test:test']
        exp = ['24687','23458','Proteobacteria_Gammaproteobacteria',\
                  'test_test_test']
        obs = map(conversion_fn,values)
        self.assertEqual(obs,exp)


    def test_remove_spaces(self):
        """spaces_to_underscores should replaces whitespace with underscores"""

        text = "E.\tcoli"
        obs = remove_spaces(text)
        exp = "E._coli"
        self.assertEqual(obs,exp)

        text = "\t\nE.\tcoli\n\t\t"
        obs = remove_spaces(text)
        exp = "E._coli"
        self.assertEqual(obs,exp)
    
    def test_yield_trait_table_fields(self):
        """yield_trait_table_fields should successively yield trait table fields"""
        pass

    def test_ensure_root_is_bifurcating(self):
        """ensure_root_is_bifurcating ensures that the root of the tree is bifuracting"""
        pass

    def test_filter_tree_tips_by_presence_in_table(self):
        """filter_tree_tips_by_presence_in_table filters appropriately"""
        tree = self.SimpleTree
        table_fields =[['A',1,0,0,1],['B',1,0,1,1]]
        obs_tree = filter_tree_tips_by_presence_in_table(tree,table_fields,\
          verbose = False)
        exp_tree = "(A:0.02,B:0.01)root;"
        #print "valid tree tips:", ",".join([t.Name for t in tree.tips()])
        #print "obs_tree:",obs_tree
        self.assertEqual(obs_tree.getNewick(with_distances=True),exp_tree)



    def test_print_node_summary_table(self):
        """print_node_summary_table prints a summary of tree nodes"""
        
        # Make sure it works for the simple tree
        tree = self.SimpleTree
        obs = [ l for l in print_node_summary_table(tree)]
        
        # Fields should be name, number of children, length, parent
        exp = ['A\t0\t0.02\tE', 'B\t0\t0.01\tE', 'E\t2\t0.05\troot',\
                'C\t0\t0.01\tF', 'D\t0\t0.01\tF', 'F\t2\t0.05\troot',\
                'root\t2\tNone\tNone']
        self.assertEqual(obs,exp) 

        # Now try it out with a tree containing a polytomy
        tree = self.SimplePolytomyTree
        obs = [ l for l in print_node_summary_table(tree)]
        
        # Fields should be name, number of children, length, parent
        exp = ['A\t0\t0.02\tE', 'B\t0\t0.01\tE','B_prime\t0\t0.03\tE',\
                'E\t3\t0.05\troot','C\t0\t0.01\tF', 'D\t0\t0.01\tF',\
                'F\t2\t0.05\troot', 'root\t2\tNone\tNone']
        self.assertEqual(obs,exp) 



    def test_add_to_filename(self):
        """add_to_filename adds text between name and file extension"""
        
        #Basic Test
        filename = "my_tree.newick"
        new_suffix = "formatted"
        delimiter = "_"

        obs = add_to_filename(filename,new_suffix,delimiter)
        exp = "my_tree_formatted.newick"
        
        self.assertEqual(obs,exp)

        #Make sure it works with multiple extensions
        
        # NOTE: this fn is based on splitext from os.path
        # therefore only the rightmost extension is considered.
        
        
        filename = "my_tree.this.is.a.newick"
        new_suffix = "formatted"
        delimiter = "."

        obs = add_to_filename(filename,new_suffix,delimiter)
        exp = "my_tree.this.is.a.formatted.newick"
        
        self.assertEqual(obs,exp)

    def test_make_id_mapping_dict(self):
        """Make id mapping dict should generate two mapping dicts"""

        mappings = self.GreengenesToIMG
        trait_to_tree_mappings =\
        make_id_mapping_dict(mappings)
        
        obs = trait_to_tree_mappings
        exp  = {
          '645058788':'469810',\
          '645058789':'457471',\
          '641736109':'266998'}
        self.assertEqualItems(obs,exp)
        
        #obs = tree_to_trait_mappings
        #exp = {
        ##  '469810':'645058788',\
        #  '457471':'645058789',\
        #  '266998':'641736109'}
        #
        #self.assertEqualItems(obs,exp)

    def remap_trait_table_organims(self):
        """Remap trait table organisms should remap ids"""
     
        trait_to_tree_mappings = {
          '645058788':'469810',\
          '645058789':'457471',\
          '641736109':'266998'}
         
        img_lines = [\
         "645058788\t0\t0\t0\n",\
         "645058789\t0\t0\t0\n",\
         "641736109\t3\t1\t1\n"]

        exp_remapped_lines = [ 
         "469810\t0\t0\t0\n",\
         "457471\t0\t0\t0\n",\
         "266998\t3\t1\t1\n"]

        obs_remapped_lines = \
          remap_trait_table_organisms(\
          trait_table_lines = img_lines,\
          trait_table_to_tree_mapping_dict = trait_to_tree_mappings,\
          input_delimiter = "\t", output_delimiter = "\t")

        self.assertEqual(exp_remapped_lines, obs_remapped_lines)
        
        # Now test data with non-mapping entries, and
        # specifying space-delimited output

        img_lines = [\
         "645058788\t0\t0\t0\n",\
         "645058789\t0\t0\t0\n",\
         "armadillo\t0\t0\t0\n",\
         "641736109\t3\t1\t1\n"]

        exp_remapped_lines = [ 
         "469810 0 0 0\n",\
         "457471 0 0 0\n",\
         "None 0 0 0\n",\
         "266998 3 1 1\n"]

        obs_remapped_lines = \
          remap_trait_table_organisms(\
          trait_table_lines = img_lines,\
          trait_table_to_tree_mapping_dict = trait_to_tree_mappings,\
          input_delimiter = "\t", output_delimiter = " ",\
          default_entry="None")

        self.assertEqual(exp_remapped_lines, obs_remapped_lines)


if __name__ == "__main__":
    main()
