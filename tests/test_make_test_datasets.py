#!/usr/bin/env python

from unittest import main, TestCase 
from cogent import LoadTree
from cogent.parse.tree import DndParser
from picrust.make_test_datasets import exclude_tip, yield_test_trees,\
  make_distance_based_exclusion_fn, make_distance_based_tip_label_randomizer
from collections import defaultdict
"""
Tests for make_test_datasets.py
"""

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
__credits__ = ["Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "2-alpha.1"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"


class TestMakeTestTrees(TestCase):
    """Tests of make_test_trees.py"""

    def setUp(self):
        self.SimpleTree = \
          DndParser("((A:0.02,B:0.01)E:0.05,(C:0.01,D:0.01)F:0.05)root;")
    
        self.SimplePolytomyTree = \
                DndParser("((A:0.02,B:0.01,B_prime:0.03)E:0.05,(C:0.01,D:0.01)F:0.05)root;")
    
    def test_exclude_tip(self):
        """exclude_tip should yield a holdout tree"""

        #Try excluding tip 'B'
        test_tree = self.SimpleTree.deepcopy()
        
        obs = exclude_tip(test_tree.getNodeMatchingName('B'),test_tree)
        obs_newick = obs.getNewick(with_distances=True)
        exp_newick = "((C:0.01,D:0.01)F:0.05,A:0.07)root;" 
        alt_newick = "(A:0.07,(C:0.01,D:0.01)F:0.05)root;"
        #exp_newick and alt_newick represent 
        #two ways of expressing the same tree
        self.assertTrue(obs_newick in [exp_newick,alt_newick])
        
        #Make sure the holdout works with a polytomy
        test_tree = self.SimplePolytomyTree.deepcopy()
        
        obs = exclude_tip(test_tree.getNodeMatchingName('B'),test_tree)
        obs_newick = obs.getNewick(with_distances=True)
        exp_newick = "((A:0.02,'B_prime':0.03)E:0.05,(C:0.01,D:0.01)F:0.05)root;" 
        self.assertEqual(obs_newick,exp_newick)

        #Make sure we raise if the tip is invalid        
         
        test_tree = self.SimpleTree.deepcopy()
        
        self.assertRaises(ValueError,exclude_tip,\
          test_tree.getNodeMatchingName('E'),test_tree)

    def test_make_distance_based_exclusion_fn(self):
        """make_distance_based_exclusion_fn should return a working function"""

        exclude_similar_strains =\
            make_distance_based_exclusion_fn(0.03)
        
        #Test that new function is documented
        exp_doc = 'Exclude neighbors of tip within 0.030000 branch length units'
        self.assertEqual(exp_doc,exclude_similar_strains.__doc__)
        
        #Test that the function works
        
        test_tree = self.SimpleTree.deepcopy()
        #print test_tree.getNewick(with_distances=True)
        tip = test_tree.getNodeMatchingName('C')
        obs = exclude_similar_strains(tip,test_tree).getNewick(with_distances=True) 
        exp = "(A:0.02,B:0.01)root;"
        self.assertEqual(obs,exp)
        

        #Test on a tree where a single node will remain
        test_tree = \
          DndParser("((A:0.02,B:0.01)E:0.05,(C:0.06,D:0.01)F:0.05)root;")
        #print test_tree.getNewick(with_distances=True)
        tip = test_tree.getNodeMatchingName('D')
        obs = exclude_similar_strains(tip,test_tree).getNewick(with_distances=True) 
        exp = "((A:0.02,B:0.01)E:0.05,C:0.11)root;"
        self.assertEqual(obs,exp)
        
        
        #Test that we raise if distance is too large
        test_tree = self.SimpleTree.deepcopy()
        test_fn = make_distance_based_exclusion_fn(300.0)
        tip = test_tree.getNodeMatchingName('C')
        
        self.assertRaises(ValueError,test_fn,tip,test_tree)


    def test_make_distance_based_randomizer(self):
        """make_distance_based_randomizer should randomize tip labels wthin d.  NOTE: results are tested by cumulative binomial distribution, so ~1/500 test runs may fail stochasitically."""

        tip_randomizer =\
            make_distance_based_tip_label_randomizer(0.50)

        #Test that the function works
        test_replicates = 5000
        results = defaultdict(int)
        total = 0
        for i in range(test_replicates):
            test_tree = self.SimpleTree.deepcopy()
            #print test_tree.asciiArt()
            tip = test_tree.getNodeMatchingName('C')
            obs = tip_randomizer(tip,test_tree)
            obs_newick = obs.getNewick(with_distances=True) 
            #print obs.asciiArt()     
        
            results[obs_newick]+= 1
            total += 1
        n_unique_trees = len(results.keys())

        #Since only the 4 tips are scrambled
        #the number of unique trees is just
        #the number of permutations 4! = 4*3*2 = 24
        #(but note that half of these are equivalent
        #trees with a,b swapped with b,a)
        self.assertEqual(n_unique_trees,24)
         
        for tree in sorted(results.keys()):
            n_obs = results[tree]
            # Treat the n trials as Bernoulli trials
            # each has success change 1/24
            # so we want the inverse cumulative
            # binomial for 5000 trials with p ~ .041666666
            # n_successes > 253 has p 0.001 (166 gives a p of 0.001
            #on the other tail)
               
            self.assertTrue(n_obs >= 166)
            self.assertTrue(n_obs < 253)   
    
    def test_yield_test_trees(self):
        """yield_test_trees should yield modified test trees"""
        
        start_tree = self.SimpleTree.deepcopy()


        # First, test with simple tree and exclude_tip
        # (not a typical use, but easy to predict
        # and still important to test)

        
        test_fn = exclude_tip
        obs = [o for o in yield_test_trees(start_tree, test_fn)]
        test_tips =  [tip for tip,tree in obs]
        test_trees = [tree for tip,tree in obs]
        #print test_tips
        #print test_trees 
        #Test that each tree  excludes the correct tip
        for i,exp_tip in enumerate(start_tree.tips()):
            #print exp_tip
            node_names = [obs_tip.Name for obs_tip in test_trees[i].tips()]
            #print node_names
            self.assertTrue(exp_tip.Name not in node_names)
            
        #Test that the topology is correct for an example tree 
        self.assertTrue(test_trees[1].getNewick(with_distances=True) in\
                ["((C:0.01,D:0.01)F:0.05,A:0.07)root;","(A:0.07,(C:0.01,D:0.01)F:0.05)root;"])
        
        #Second, let's test a more realistic test function
        test_fn = make_distance_based_exclusion_fn(0.03)
        
        obs = [o for o in yield_test_trees(start_tree, test_fn)]
        test_tips =  [tip for tip,tree in obs]
        test_trees = [tree for tip,tree in obs]
        
        #Test that each tree excludes the tip of interest 
        #NOTE: this behavior changed some time ago
        #We want to exclude the test tip.
        for i,exp_tip in enumerate(start_tree.tips()):
            #print exp_tip
            node_names = [obs_tip.Name for obs_tip in test_trees[i].tips()]
            #print node_names
            self.assertFalse(exp_tip.Name in node_names)
            
        #Test that the topology is correct for all tips
        
        
        #Tip 'A'
        self.assertEqual(test_trees[0].getNewick(with_distances=True),\
            "(C:0.01,D:0.01)root;")
        
        #Tip 'B'
        self.assertEqual(test_trees[1].getNewick(with_distances=True),\
            "(C:0.01,D:0.01)root;")
        
        #Tip 'C'
        self.assertEqual(test_trees[2].getNewick(with_distances=True),\
            "(A:0.02,B:0.01)root;")

        #Tip 'D'
        self.assertEqual(test_trees[3].getNewick(with_distances=True),\
            "(A:0.02,B:0.01)root;")











if __name__ == "__main__":
    main()
