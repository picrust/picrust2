#!/usr/bin/env python

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
__credits__ = ["Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "2-alpha.1"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"

from math import e,sqrt
from cogent.util.unit_test import main,TestCase
from numpy import array,arange,array_equal,around
from cogent import LoadTree
from cogent.parse.tree import DndParser
from cogent.app.util import get_tmp_filename
from cogent.util.misc import remove_files
from cogent.maths.stats.special import ndtri
from warnings import catch_warnings
from picrust.predict_traits  import assign_traits_to_tree,\
  predict_traits_from_ancestors, get_most_recent_ancestral_states,\
  fill_unknown_traits, equal_weight,linear_weight,\
  inverse_variance_weight, make_neg_exponential_weight_fn,\
  weighted_average_tip_prediction, get_interval_z_prob,\
  thresholded_brownian_probability, update_trait_dict_from_file,\
  biom_table_from_predictions, get_nearest_annotated_neighbor,\
  predict_nearest_neighbor, predict_random_neighbor,\
  calc_nearest_sequenced_taxon_index,\
  variance_of_weighted_mean,fit_normal_to_confidence_interval,\
  get_most_recent_reconstructed_ancestor,\
  normal_product_monte_carlo, get_bounds_from_histogram,\
  get_nn_by_tree_descent,get_brownian_motion_param_from_confidence_intervals


"""
Tests for predict_traits.py
"""


in_trait1="""nodes	trait2	trait1
3	3	1
A	5	2.5
D	5	2"""

in_trait2="""tips	trait1	trait2	trait3
1	1	3	1
2	0	3	2
3	2	3	3"""

in_bad_trait="""tips	not_trait1	trait2	trait3
1	1	3	1
2	0	3	2
3	2	3	3"""



class TestPredictTraits(TestCase):
    """Tests of predict_traits.py"""

    def setUp(self):
        self.SimpleTree = \
          DndParser("((A:0.02,B:0.01)E:0.05,(C:0.01,D:0.01)F:0.05)root;")
        
        
        #Set up a tree with obvious differences in the rate of gene content
        #evolution to test confidence interval estimation
        #Features:  
        # --trait 1 is has ~ 10 fold higher confidence intervals than trait 0. 
        # Trait 2 is 10 fold higher than trait 1
        
        # -- of predicted nodes B and D, D has a ~10 fold longer branch

        self.SimpleUnequalVarianceTree =\
          DndParser("((A:0.01,B:0.01)E:0.05,(C:0.01,D:0.10)F:0.05)root;")
        traits = {"A":[1.0,1.0,1.0],"C":[1.0,1.0,1.0],"E":[1.0,1.0,1.0],"F":[1.0,1.0,1.0]}
        self.SimpleUnequalVarianceTree = assign_traits_to_tree(traits,\
          self.SimpleUnequalVarianceTree,trait_label="Reconstruction")
        self.SimpleUnequalVarianceTree.getNodeMatchingName('E').upper_bound = [2.0,20.0,200.0]
        self.SimpleUnequalVarianceTree.getNodeMatchingName('E').lower_bound = [-1.0,-19.0,-199.0]
        self.SimpleUnequalVarianceTree.getNodeMatchingName('F').upper_bound = [2.0,20.0,200.0]
        self.SimpleUnequalVarianceTree.getNodeMatchingName('F').lower_bound = [-1.0,-19.0,-199.0]
        
        #Set up a tree with a three-way polytomy
        self.SimplePolytomyTree = \
          DndParser("((A:0.02,B:0.01,B_prime:0.03)E:0.05,(C:0.01,D:0.01)F:0.05)root;")
    
        self.SimpleTreeTraits =\
            {"A":[1.0,1.0],"E":[1.0,1.0],"F":[0.0,1.0],"D":[0.0,0.0]}
        
        self.PartialReconstructionTree =\
                DndParser("((((B:0.01,C:0.01)I3:0.01,A:0.01)I2:0.01,D:0.01)I1:0.01)root;")

        self.CloseToI3Tree =\
                DndParser("((((B:0.01,C:0.95)I3:0.01,A:0.01)I2:0.95,D:0.05)I1:0.95)root;")
        
        self.CloseToI1Tree =\
                DndParser("((((B:0.95,C:0.95)I3:0.95,A:0.01)I2:0.02,D:0.05)I1:0.05)root;")

        self.BetweenI3AndI1Tree=\
                DndParser("((((B:0.01,C:0.1)I3:0.02,A:0.01)I2:0.02,D:0.05)I1:0.02)root;")


        self.PartialReconstructionTraits =\
                {"B":[1.0,1.0],"C":[1.0,1.0],"I3":[1.0,1.0],"I1":[0.0,1.0],"D":[0.0,1.0]}

        self.GeneCountTraits =\
                {"B":[1.0,1.0],"C":[1.0,2.0],"I3":[1.0,1.0],"I1":[0.0,3.0],"D":[0.0,5.0]}

        #create a tmp trait file
        self.in_trait1_fp = get_tmp_filename(prefix='Predict_Traits_Tests',suffix='.tsv')
        self.in_trait1_file=open(self.in_trait1_fp,'w')
        self.in_trait1_file.write(in_trait1)
        self.in_trait1_file.close()

        #create another tmp trait file (with columns in different order)
        self.in_trait2_fp = get_tmp_filename(prefix='Predict_Traits_Tests',suffix='.tsv')
        self.in_trait2_file=open(self.in_trait2_fp,'w')
        self.in_trait2_file.write(in_trait2)
        self.in_trait2_file.close()


        #create a tmp trait file with a incorrect trait name
        self.in_bad_trait_fp = get_tmp_filename(prefix='Predict_Traits_Tests',suffix='.tsv')
        self.in_bad_trait_file=open(self.in_bad_trait_fp,'w')
        self.in_bad_trait_file.write(in_bad_trait)
        self.in_bad_trait_file.close()

        self.files_to_remove = [self.in_trait1_fp,self.in_trait2_fp,self.in_bad_trait_fp]

    def tearDown(self):
        remove_files(self.files_to_remove)
    
    def test_nearest_neighbor_prediction(self):
        """nearest_neighbor_prediction predicts nearest neighbor's traits"""
        traits = self.SimpleTreeTraits
        tree = self.SimpleTree
        result_tree = assign_traits_to_tree(traits,tree,trait_label="Reconstruction")
        
        #Test with default options
        results = predict_nearest_neighbor(tree, nodes_to_predict =["B","C"])
        self.assertEqual(results["B"],array([1.0,1.0]))
        self.assertEqual(results["C"],array([0.0,0.0]))
        
        #Test allowing ancestral NNs
        results = predict_nearest_neighbor(tree, nodes_to_predict =["B","C"],\
         tips_only = False)
        self.assertEqual(results["C"],array([0.0,1.0]))

        #Test allowing self to be NN AND Ancestral NNs
        results = predict_nearest_neighbor(tree, nodes_to_predict =["A","B","C","D"],\
         tips_only = False,use_self_in_prediction=True)

        self.assertEqual(results["A"],array([1.0,1.0]))
        self.assertEqual(results["B"],array([1.0,1.0]))
        self.assertEqual(results["C"],array([0.0,1.0]))
        self.assertEqual(results["D"],array([0.0,0.0]))

 
    def test_calc_nearest_sequenced_taxon_index(self):
        """calc_nearest_sequenced_taxon_index calculates the NSTI measure"""
        traits = self.SimpleTreeTraits
        tree = self.SimpleTree
        result_tree = assign_traits_to_tree(traits,tree,trait_label="Reconstruction")
        #Expected distances:
        # A --> A 0.0
        # B --> A 0.03
        # C --> D 0.02
        # D --> D 0.0
        # = 0.05/4.0 = 0.0125
        exp = 0.0125
        verbose = False
        #Test with default options
        obs_nsti,obs_distances = calc_nearest_sequenced_taxon_index(tree,verbose=verbose)
        self.assertFloatEqual(obs_nsti,exp)
        self.assertFloatEqual(obs_distances["A"],0.0)
        self.assertFloatEqual(obs_distances["B"],0.03)
        self.assertFloatEqual(obs_distances["C"],0.02)
        self.assertFloatEqual(obs_distances["D"],0.00)

        #Test calcing the index while 
        #limiting prediction to B and C
        
        # B --> A 0.03
        # C --> D 0.02
        
        exp = 0.025
        obs_nsti,obs_distances = calc_nearest_sequenced_taxon_index(tree,\
          limit_to_tips = ["B","C"],verbose=False)
        self.assertFloatEqual(obs_nsti,exp)
        self.assertFloatEqual(obs_distances["B"],0.03)
        self.assertFloatEqual(obs_distances["C"],0.02)
    
    def test_get_nn_by_tree_descent(self):
        """calc_nearest_sequenced_taxon_index calculates the NSTI measure"""
        traits = self.SimpleTreeTraits
        tree = self.SimpleTree
        result_tree = assign_traits_to_tree(traits,tree,trait_label="Reconstruction")
        #Expected distances:
        # A --> A 0.0
        # B --> A 0.03
        # C --> D 0.02
        # D --> D 0.0
        # = 0.05/4.0 = 0.0125
        exp = 0.0125
        #Test with default options
        nn,distance = get_nn_by_tree_descent(tree,"B",verbose=True)
        self.assertEqual(nn.Name,"A")
        self.assertFloatEqual(distance,0.03)
        
        nn,distance = get_nn_by_tree_descent(tree,"A",verbose=True)
        self.assertEqual(nn.Name,"A")
        self.assertFloatEqual(distance,0.00)
        
        nn,distance = get_nn_by_tree_descent(tree,"A",filter_by_property=False,verbose=True)
        self.assertEqual(nn.Name,"B")
        self.assertFloatEqual(distance,0.03)
        
        nn,distance = get_nn_by_tree_descent(tree,"C",verbose=True)
        self.assertEqual(nn.Name,"D")
        self.assertFloatEqual(distance,0.02)
        #self.assertFloatEqual(obs_distances["A"],0.0)
        #self.assertFloatEqual(obs_distances["B"],0.03)
        #self.assertFloatEqual(obs_distances["C"],0.02)
        #self.assertFloatEqual(obs_distances["D"],0.00)

        #Test calcing the index while 
        #limiting prediction to B and C
        
        # B --> A 0.03
        # C --> D 0.02
        
        exp = 0.025
        obs_nsti,obs_distances = calc_nearest_sequenced_taxon_index(tree,\
          limit_to_tips = ["B","C"],verbose=False)
        self.assertFloatEqual(obs_nsti,exp)
        self.assertFloatEqual(obs_distances["B"],0.03)
        self.assertFloatEqual(obs_distances["C"],0.02)


    def test_predict_random_neighbor(self):
        """predict_random_neighbor predicts randomly"""
        traits = self.SimpleTreeTraits
        tree = self.SimpleTree
        result_tree = assign_traits_to_tree(traits,tree)
        
        #If there is only one other valid result, this
        #should always be predicted
        
        #self.SimpleTreeTraits =\
        #            {"A":[1.0,1.0],"E":[1.0,1.0],"F":[0.0,1.0],"D":[0.0,0.0]}
        
        #If self predictions are disallowed, then the prediction for A should
        #always come from node D, and be 0,0.   

        results = predict_random_neighbor(tree,['A'],\
          trait_label = "Reconstruction",\
          use_self_in_prediction=False)

        self.assertEqual(results['A'],[0.0,0.0])

        #If use_self is True, ~50% of predictions should be [1.0,1.0] and
        # half should be [0.0,0.0]

        #Pick repeatedly and make sure frequencies are
        #reasonable.  The technique is fast, so 
        #many iterations are reasonable.
        
        iterations = 100000
        a_predictions = 0
        d_predictions = 0
        for i in range(iterations):
            results = predict_random_neighbor(tree,['A'],\
              trait_label = "Reconstruction",\
              use_self_in_prediction=True)
            #print results
            if results['A'] == [1.0,1.0]:
                #print "A pred"
                a_predictions += 1
            elif results['A'] == [0.0,0.0]:
                #print "D pred"
                d_predictions +=1
            else:
                raise RuntimeError(\
                  "Bad prediction result: Neither node A nor node D traits used in prediction")
        #print "All a predictions:",a_predictions
        #print "All d predictions:",d_predictions
        ratio = float(a_predictions)/float(iterations)
        #print "Ratio:", ratio
        self.assertFloatEqual(ratio,0.5,eps=1e-2)




    def test_get_nearest_annotated_neightbor(self):
        """get_nearest_annotated_neighbor finds nearest relative with traits"""
        traits = self.SimpleTreeTraits
        tree = self.SimpleTree
        result_tree = assign_traits_to_tree(traits,tree)
 

       
        #Test ancestral NN matching
        nn =  get_nearest_annotated_neighbor(tree,'A',\
              tips_only=False, include_self=False)
        
        self.assertEqual(nn.Name,'E')
        
        nn =  get_nearest_annotated_neighbor(tree,'B',\
              tips_only=False, include_self=False)
        
        self.assertEqual(nn.Name,'E')
        
 
        nn =  get_nearest_annotated_neighbor(tree,'C',\
              tips_only=False, include_self=False)
        
        self.assertEqual(nn.Name,'F')
        
  
        nn =  get_nearest_annotated_neighbor(tree,'D',\
              tips_only=False, include_self=False)
        
        self.assertEqual(nn.Name,'F')
        
       
        #Test tip only, non-self matching
        nn =  get_nearest_annotated_neighbor(tree,'A',\
              tips_only=True, include_self=False)
        
        self.assertEqual(nn.Name,'D')
        
        nn =  get_nearest_annotated_neighbor(tree,'B',\
              tips_only=True, include_self=False)
        
        self.assertEqual(nn.Name,'A')

 
        nn =  get_nearest_annotated_neighbor(tree,'C',\
              tips_only=True, include_self=False)
        
        self.assertEqual(nn.Name,'D')

        nn =  get_nearest_annotated_neighbor(tree,'D',\
              tips_only=True, include_self=False)
        
        self.assertEqual(nn.Name,'A')

    def test_biom_table_from_predictions(self):
        """format predictions into biom format"""
        traits = self.SimpleTreeTraits
        tree = self.SimpleTree
        
        #print "Starting tree:",tree.asciiArt()
        # Test on simple tree
        result_tree = assign_traits_to_tree(traits,tree)
        nodes_to_predict = [n.Name for n in result_tree.tips()]
        #print "Predicting nodes:", nodes_to_predict
        predictions = predict_traits_from_ancestors(result_tree,\
          nodes_to_predict)

        biom_table=biom_table_from_predictions(predictions,["trait1","trait2"])
        
    def test_equal_weight(self):
        """constant_weight weights by a constant"""
        w = 1.0
        d = 0.1
        for i in range(100):
            obs = equal_weight(i)
            exp = w
            self.assertFloatEqual(obs,exp)
    
    def test_make_neg_exponential_weight_fn(self):
        """make_neg_exponential_weight_fn returns the specified fn"""
        
        exp_base = 10
        weight_fn = make_neg_exponential_weight_fn(exp_base)
        
        d = 10.0
        obs = weight_fn(d)
        exp = 10.0**-10.0
        self.assertFloatEqual(obs,exp)

        #Test for base two
        exp_base = 2
        weight_fn = make_neg_exponential_weight_fn(exp_base)
        
        d = 16.0
        obs = weight_fn(d)
        exp = 2.0**-16.0
        self.assertFloatEqual(obs,exp)


    def test_linear_weight(self):
        """linear_weight weights linearly"""
        
        max_d = 1.0
        d = 0.90
        obs = linear_weight(d,max_d)
        exp = 0.10
        self.assertFloatEqual(obs, exp)

        d = 0.0
        obs = linear_weight(d,max_d)
        exp = 1.0
        self.assertFloatEqual(obs, exp)

        max_d = 3.0
        d = 1.5
        obs = linear_weight(d,max_d)
        exp = 0.50
        self.assertFloatEqual(obs, exp)
    
    def test_inverse_variance_weight(self):
        """inverse_variance_weight"""
        #TODO: test this works with arrays of variances 
        var = 1000.0
        for d in range(1,10):
            d = float(d)
            obs = inverse_variance_weight(d,var)
            exp = 1.0/1000.0
            self.assertFloatEqual(obs,exp)

        #Now test the special case of zero variance
        var = 0.0
        for d in range(1,10):
            d = float(d)
            obs = inverse_variance_weight(d,var)
            exp = 1.0/1e-10
            self.assertFloatEqual(obs,exp)


    def test_assign_traits_to_tree(self):
        """assign_traits_to_tree should map reconstructed traits to tree nodes"""
        
        # Test that the function assigns traits from a dict to a tree node
        traits = self.SimpleTreeTraits
        tree = self.SimpleTree
        
        # Test on simple tree
        result_tree = assign_traits_to_tree(traits,tree)
        
        # Test that each node is assigned correctly
        for node in result_tree.preorder():
            obs = node.Reconstruction 
            exp = traits.get(node.Name, None)
            self.assertEqual(obs,exp)
        
        # Test on polytomy tree
        
        tree = self.SimplePolytomyTree
        result_tree = assign_traits_to_tree(traits,tree)
        
        # Test that each node is assigned correctly
        for node in result_tree.preorder():
            obs = node.Reconstruction 
            exp = traits.get(node.Name, None)
            self.assertEqual(obs,exp)
    
    def test_assign_traits_to_tree_quoted_node_name(self):
        """Assign_traits_to_tree should remove quotes from node names"""
        # Test that the function assigns traits from a dict to a tree node
        traits = self.SimpleTreeTraits
        tree = self.SimpleTree
        #Make one node quoted
        tree.getNodeMatchingName('A').Name="'A'"
        tree.getNodeMatchingName('B').Name='"B"'

        # Test on simple tree
        result_tree = assign_traits_to_tree(traits,tree,fix_bad_labels=True)
        #Setting fix_bad_labels to false produces NoneType predictions when
        #labels are quoted
        
        # Test that each node is assigned correctly
        for node in result_tree.preorder():
            obs = node.Reconstruction 
            exp = traits.get(node.Name.strip("'").strip('"'), None)
            self.assertEqual(obs,exp)
        
        # Test on polytomy tree
        
        tree = self.SimplePolytomyTree
        result_tree = assign_traits_to_tree(traits,tree)
        
        # Test that each node is assigned correctly
        for node in result_tree.preorder():
            obs = node.Reconstruction 
            exp = traits.get(node.Name, None)
            self.assertEqual(obs,exp)

    def test_update_trait_dict_from_file(self):
        """update_trait_dict_from_file should parse input trait tables (asr and genome) and match traits between them"""
        header,traits=update_trait_dict_from_file(self.in_trait1_fp)
        self.assertEqual(header,["trait2","trait1"])
        self.assertEqual(traits,{3:[3,1],'A':[5,2.5],'D':[5,2]})

        #test that we get a warning when header from other trait table doesn't match perfectly.
        with catch_warnings(record=True) as w:
            header2,traits2=update_trait_dict_from_file(self.in_trait2_fp,header)
            self.assertEqual(header2,["trait2","trait1"])
            self.assertEqual(traits2,{1:[3,1], 2:[3,0], 3:[3,2]})
            assert len(w) == 1
            assert issubclass(w[-1].category, UserWarning)
            assert "Missing" in str(w[-1].message)
                    

        #try giving a trait table with a trait that doesn't match our header
        self.assertRaises(RuntimeError,update_trait_dict_from_file,self.in_bad_trait_fp,header)

    def test_predict_traits_from_ancestors(self):
        """predict_traits_from_ancestors should propagate ancestral states"""
        # Testing the point predictions first (since these are easiest) 
        # When the node is very close to I3, prediction should be approx. I3

        traits = self.PartialReconstructionTraits
        tree = assign_traits_to_tree(traits,self.CloseToI3Tree)
        
        nodes_to_predict = ['A'] 
        prediction = predict_traits_from_ancestors(tree=tree,\
          nodes_to_predict=nodes_to_predict) 
        
        exp = traits["I3"]
        #print "PREDICTION:",prediction 
        for node in nodes_to_predict:
            self.assertFloatEqual(around(prediction[node]),exp)

        #TODO: need to add test case where a very hard to predict
        # single value is present in a sequenced genome.  Then
        # test that use_self_in_prediction controls whether this is used
        

    def test_predict_traits_from_ancestors_correctly_predicts_variance(self):
        """predict_traits_from_ancestors should correctly report variance due to branch lengths and rates of gene copy number evolution """
        tree = self.SimpleUnequalVarianceTree
        #All values are 1, but variance in the prediction should vary
        #due to vary unequal branch lengths (between taxa) and brownian
        #motion parameters (between traits)
        nodes_to_predict = ['B','D']
        bm_fixed_10_fold = [1.0,10.0,100.0]
        prediction,variances,confidence_intervals = predict_traits_from_ancestors(tree=tree,\
          nodes_to_predict=nodes_to_predict,calc_confidence_intervals=True,\
          lower_bound_trait_label='lower_bound',upper_bound_trait_label='upper_bound',
          brownian_motion_parameter = bm_fixed_10_fold,trait_label="Reconstruction")
        
        #All traits are 1, so all predictions should be 1
        exp_predictions = {'B':[1.0,1.0,1.0],'D':[1.0,1.0,1.0]}
        self.assertEqualItems(prediction,exp_predictions)
        #We don't expect variances to be exactly 10 fold increasing
        #but do expect they should be in rank order
        for tip in ['B','D']:
            tip_vars = variances[tip]['variance']
            self.assertTrue(tip_vars[0]<tip_vars[1]) 
            self.assertTrue(tip_vars[1]<tip_vars[2])
        
        #Also note that trait D is on a much longer branch, so we expect
        #it to have higher variance
        self.assertTrue((array(variances['B']['variance'])<array(variances['D']['variance'])).all())
            
    
    
    def test_fill_unknown_traits(self):
        """fill_unknown_traits should propagate only known characters"""


        # Only the missing values in to_update should be 
        # filled in with appropriate values from new
        to_update = array([1,0,1,None,1,0])
        new = array([None,None,1,1,1,1])
    
        obs = fill_unknown_traits(to_update,new)
        exp = array([1,0,1,1,1,0])

        self.assertTrue(array_equal(obs,exp))

        #Try the reverse update

        obs = fill_unknown_traits(new,to_update)
        exp = array([1,0,1,1,1,1])
        self.assertTrue(array_equal(obs,exp))

        # Ensure that if to_update is None, the value of new is returned
        obs = fill_unknown_traits(None, new)
        #print "Obs:",obs
        exp = new
        self.assertTrue(array_equal(obs,exp))

    def test_weighted_average_tip_prediction(self):
        """Weighted average node prediction should predict node values"""
        
        
        # When the node is very close to I3, prediction should be approx. I3

        traits = self.PartialReconstructionTraits
        tree = assign_traits_to_tree(traits,self.CloseToI3Tree)
        
        node_to_predict = "A"
        node = tree.getNodeMatchingName(node_to_predict)
        most_recent_reconstructed_ancestor =\
          get_most_recent_reconstructed_ancestor(node)
        
        prediction = weighted_average_tip_prediction(tree=tree,\
          node=node,\
          most_recent_reconstructed_ancestor=\
          most_recent_reconstructed_ancestor)
            
        
        exp = traits["I3"]
        
        self.assertFloatEqual(around(prediction),exp)


        # When the node is very close to I1, prediction should be approx. I1


        traits = self.PartialReconstructionTraits
        tree = assign_traits_to_tree(traits,self.CloseToI1Tree)
        node_to_predict = "A"
        #print "tree:",tree.asciiArt()
        node = tree.getNodeMatchingName(node_to_predict)
        most_recent_reconstructed_ancestor =\
          get_most_recent_reconstructed_ancestor(node)
        prediction = weighted_average_tip_prediction(tree=tree,\
          node=node,\
          most_recent_reconstructed_ancestor=\
          most_recent_reconstructed_ancestor)
        exp = traits["I1"]
        #print "prediction:",prediction
        #print "exp:",exp
        a_node = tree.getNodeMatchingName('A')
        #for node in tree.preorder():
        #    print node.Name,node.distance(a_node),node.Reconstruction
        self.assertFloatEqual(around(prediction),exp)

        # Try out the B case with exponential weighting
        
        traits = self.PartialReconstructionTraits
        tree = assign_traits_to_tree(traits,self.CloseToI3Tree)
        weight_fn = make_neg_exponential_weight_fn(exp_base=e)
        
        
        node_to_predict = "A"
        node = tree.getNodeMatchingName(node_to_predict)
        most_recent_reconstructed_ancestor =\
          get_most_recent_reconstructed_ancestor(node)
        prediction = weighted_average_tip_prediction(tree=tree,\
          node=node,\
          most_recent_reconstructed_ancestor=\
          most_recent_reconstructed_ancestor)

        #prediction = weighted_average_tip_prediction(tree=tree,\
        #  node_to_predict=node_to_predict,weight_fn=weight_fn) 
        exp = traits["B"]
        self.assertFloatEqual(around(prediction),exp)

        # Try out the I1 case with exponential weighting
        
        traits = self.PartialReconstructionTraits
        tree = assign_traits_to_tree(traits,self.CloseToI1Tree)
        weight_fn = make_neg_exponential_weight_fn(exp_base=e)
        #weight_fn = linear_weight
        
        node_to_predict = "A"
        node = tree.getNodeMatchingName(node_to_predict)
        most_recent_reconstructed_ancestor =\
          get_most_recent_reconstructed_ancestor(node)
        prediction = weighted_average_tip_prediction(tree=tree,\
          node=node,\
          most_recent_reconstructed_ancestor=\
          most_recent_reconstructed_ancestor)

        exp = traits["I1"]
        self.assertFloatEqual(around(prediction),exp)

        # Try out the balanced case where children and ancestors 
        # should be weighted a equally with exponential weighting
        
        # We'll  try this with full gene count data to ensure 
        # that case is tested

        traits = self.GeneCountTraits
        tree = assign_traits_to_tree(traits,self.BetweenI3AndI1Tree)
        weight_fn = make_neg_exponential_weight_fn(exp_base=e)
        
        node_to_predict = "A"
        
        node = tree.getNodeMatchingName(node_to_predict)
        most_recent_reconstructed_ancestor =\
          get_most_recent_reconstructed_ancestor(node)
        prediction = weighted_average_tip_prediction(tree=tree,\
          node=node,\
          most_recent_reconstructed_ancestor=\
          most_recent_reconstructed_ancestor)


        
        
        
        #prediction = weighted_average_tip_prediction(tree=tree,\
        #  node_to_predict=node_to_predict,weight_fn=weight_fn) 
        
        exp = (array(traits["I1"]) + array(traits["I3"]))/2.0
        self.assertFloatEqual(prediction,exp)
        
        #TODO: test the case with partial missing data (Nones)

        #TODO: test the case with fully missing data for either
        # the ancestor or the children. 

        #TODO: Test with polytomy trees

        # These *should* work, but until they're tested we don't know

    def test_get_interval_z_prob(self):
        """get_interval_z_prob should get the probability of a Z-score on an interval"""

        #Approximate expected values were calculated from
        #the table of z-values found in:
        
        #Larson, Ron; Farber, Elizabeth (2004). 
        #Elementary Statistics: Picturing the World. P. 214, 
        #As recorded here: http://en.wikipedia.org/wiki/Standard_normal_table

        #-- Test 1 --
        #Expected values for 0 - 0.01

        obs = get_interval_z_prob(0.0,0.01)
        #Larson & Farber reported values are:
        #For z of 0.00, p= 0.5000
        #For z of 0.01, p= 0.5040
        
        exp = 0.0040
        self.assertFloatEqual(obs,exp,eps=0.01)
        #Error is around 1e-5 from estimate
        
        #-- Test 2 --
        # 0.75 - 0.80
        obs = get_interval_z_prob(0.75,0.80)
        #Larson & Farber reported values are:
        #For z of 0.75, p= 0.7734
        #For z of 0.80, p= 0.7881

        exp = 0.7881 - 0.7734
        self.assertFloatEqual(obs,exp,eps=0.01)

    def test_thresholded_brownian_probability(self):
        """Brownian prob should return dict of state probabilities"""
        #x = thresholded_brownian_probability(2.2755, 1001**0.5, 0.03, min_val = 0.0,increment = 1.00,trait_prob_cutoff = 1e-4)
        #lines =  ["\t".join(map(str,[k,x[k]]))+"\n" for k in sorted(x.keys())]
        #for line in lines:
        #    print line
        
        #print "Total prob:", sum(x.values())
        start_state = 3.0
        var = 30.00
        d = 0.03
        min_val = 0.0
        increment = 1.0
        trait_prob_cutoff =  1e-200

        obs = thresholded_brownian_probability(start_state,d,var,min_val,\
          increment,trait_prob_cutoff)
        #TODO: Need to calculate exact values for this minimal case 
        #with the Larson & Farber Z tables, by hand.
        
        #For now test for sanity
        
        #Test that no probabilities below threshold are included
        self.assertTrue(min(obs.values()) > trait_prob_cutoff)
        #Test that start values +1 or -1 are equal
        self.assertEqual(obs[2.0],obs[4.0])
        #Test that the start state is the highest prob value
        self.assertEqual(max(obs.values()),obs[start_state])
        

    def test_fit_normal_to_confidence_interval(self):
        """fit_normal_to_confidence_interval should return a mean and variance given CI"""

        #Lets use a normal distribution to generate test values
        normal_95 = ndtri(0.95)
        mean = 0
        upper = mean + normal_95
        lower = mean - normal_95
        obs_mean,obs_var =\
          fit_normal_to_confidence_interval(upper,lower,confidence=0.95)
        exp_mean = mean
        exp_var = 1.0
        self.assertFloatEqual(obs_mean,exp_mean)
        self.assertFloatEqual(obs_var,exp_var)
        
        #An alternative normal:
        normal_99 = ndtri(0.99)
        mean = 5.0
        upper = mean + normal_99
        lower = mean - normal_99
        obs_mean,obs_var =\
          fit_normal_to_confidence_interval(upper,lower,confidence=0.99)
        exp_mean = mean
        exp_var = 1.0
        self.assertFloatEqual(obs_mean,exp_mean)
        self.assertFloatEqual(obs_var,exp_var)
    
    def test_variance_of_weighted_mean(self):
        """variance_of_weighted_mean calculates the variance of a weighted mean"""
        
        #Just a hand calculated example using the formula from here:
        #http://en.wikipedia.org/wiki/Weighted_mean
       

        #TODO: test if this works for arrays of variances

        #If all weights and standard deviations are equal, then
        #variance = stdev/sqrt(n)
        weights = array([0.5,0.5])
        sample_stdevs = array([4.0,4.0])
        variances = sample_stdevs**2
        exp = 4.0/sqrt(2.0)
        obs = variance_of_weighted_mean(weights,variances)
        self.assertFloatEqual(obs,exp)

        #If standard deviations are equal, but weights are not, the result
        #is equal to stdev*sqrt(sum(squared_weights))

        weights = array([0.1,0.9])
        sample_stdevs = array([4.0,4.0])
        variances = sample_stdevs**2
        exp_unbalanced = 4.0*sqrt(sum(weights**2))
        obs = variance_of_weighted_mean(weights,variances)
        self.assertEqual(obs,exp_unbalanced)

        #If all standard deviations are equal:
        #The minimal value for the variance is when all weights are equal
        #the maximal value is when one weight is 1.0 and another is 0.0

        sample_variances = array([3.0,3.0,3.0,3.0])
        
        balanced_weights = array([0.25,0.25,0.25,0.25])
        two_weights = array([0.0,0.50,0.50,0.0])
        unbalanced_weights = array([0.0,1.0,0.0,0.0])

        balanced_variance = variance_of_weighted_mean(balanced_weights,sample_variances)
        two_weight_variance = variance_of_weighted_mean(two_weights,sample_variances)
        unbalanced_variance = variance_of_weighted_mean(unbalanced_weights,sample_variances)
        
        #We expect balanced_variance < two-weight_variance < unbalanced_variance
        self.assertTrue(balanced_variance < two_weight_variance)
        self.assertTrue(balanced_variance < unbalanced_variance)
        self.assertTrue(two_weight_variance < unbalanced_variance)


        #Check that doing this for two 1D arrays is equal to using a single 2d array
        weights1 = array([0.1,0.9])
        weights2 = array([0.5,0.5])
        vars1 = array([4.0,4.0])
        vars2 = array([1000.0,1000.0])
        obs1 = variance_of_weighted_mean(weights1,vars1)
        obs2 = variance_of_weighted_mean(weights2,vars2)
        
        #Expect that calculating the result as a single 2D array
        #gives identical results to calculating as two 1D arrays
        exp = array([obs1,obs2])
        
        combined_weights = array([[0.1,0.9],[0.5,0.5]])
        combined_vars = array([[4.0,4.0],[1000.0,1000.0]])
        combined_obs = variance_of_weighted_mean(combined_weights,combined_vars)

        self.assertFloatEqual(combined_obs,exp)

        
        
    def test_normal_product_monte_carlo(self):
        """normal_product_monte_carlo calculates the confidence limits of two normal distributions empirically"""
        
        # Need good test data here.  
        #The APPL statistical language apparently has an analytical
        # solution to the product normal that could be used

        #Result for product of two standard normal distributions
        lower,upper = normal_product_monte_carlo(0.0,1.0,0.0,1.0)
        #print "95% confidence limit for product of two standard normal distributions:",lower,upper
       # 1.60 corresponds to the value for the 0.10 (10%) confidence limit
       #when using a two-tailed test.
       #Therefore for the one tailed upper limit, I believe we expect 1.60 to 
       #correspond to a type I error rate of 0.05

        #self.assertFloatEqual(lower,-1.60,eps=.1)
        #self.assertFloatEqual(upper,1.60,eps=.1)

        #result = normal_product_monte_carlo(1.0/3.0,1.0,2.0,1.0)
        #print result
        mean1 = 0.4
        mean2 = 1.2
        v1 = 1.0
        v2 = 1.0
        lower,upper = normal_product_monte_carlo(mean1,v1,mean2,\
          v2,confidence=0.95)
        #print "confidence limit for product of two normal distributions:",\
        #    lower,upper

        lower_estimate = mean1*mean2 + lower
        upper_estimate = mean1*mean2 + upper
        #self.assertFloatEqual(lower_estimate,-1.8801,eps=.1)
        #self.assertFloatEqual(upper_estimate,2.3774,eps=.1)


    def test_get_bounds_from_histogram(self):
        """Get bounds from histogram finds upper and lower tails of distribution at specified confidence levels"""
        
        #Test a simple array

        test_hist = array([0.01,0.98,0.01])
        test_bin_edges = arange(3)
        obs_lower,obs_upper = get_bounds_from_histogram(test_hist,test_bin_edges,confidence=0.90)
        #Upper and lower bounds should be conservative, and therefore exclude the center
        exp_lower = 1
        exp_upper = 2
        self.assertFloatEqual(obs_lower,exp_lower)
        self.assertFloatEqual(obs_upper,exp_upper)
        
        # Confirm that summing the histogram over given indices
        # gives <= confidence % of the mass

        obs_sum_lower = sum(test_hist[:obs_lower])
        self.assertTrue(obs_sum_lower <= 0.05*sum(test_hist))
        obs_sum_upper = sum(test_hist[obs_upper:])
        self.assertTrue(obs_sum_upper <= 0.05*sum(test_hist))

        #Repeat for a more complex test case

        test_hist =array([1.0,2.0,0.0,5.0,25.0,2.0,50.0,10.0,5.0,1.0])
        test_bin_edges = array(arange(len(test_hist)+1))
        obs_lower,obs_upper = get_bounds_from_histogram(test_hist,test_bin_edges,confidence=0.90)
        
        exp_lower = 3
        exp_upper = 9
        self.assertFloatEqual(obs_lower,exp_lower)
        self.assertFloatEqual(obs_upper,exp_upper)

        obs_sum_lower = sum(test_hist[:obs_lower])
        self.assertTrue(obs_sum_lower <= 0.05*sum(test_hist))
        obs_sum_upper = sum(test_hist[obs_upper:])
        self.assertTrue(obs_sum_upper <= 0.05*sum(test_hist))


    
    def test_get_brownian_motion_param_from_confidence_intervals(self):
        """Get brownian motion parameters from confidence intervals"""
        #TODO: Ensure this works with arrays of brownian motions

        tree = self.SimpleTree
        
        #Test one-trait case
        traits = {"A":[1.0],"C":[2.0],"E":[1.0],"F":[1.0]}
        tree = assign_traits_to_tree(traits,tree,trait_label="Reconstruction") 
        tree.getNodeMatchingName('E').upper_bound = [2.0]  
        tree.getNodeMatchingName('F').upper_bound = [1.0]
        tree.getNodeMatchingName('E').lower_bound = [0.0]  
        tree.getNodeMatchingName('F').lower_bound = [1.0]
        
        brownian_motion_parameter =\
          get_brownian_motion_param_from_confidence_intervals(tree,\
          upper_bound_trait_label="upper_bound",\
          lower_bound_trait_label="lower_bound",\
          trait_label="Reconstruction",\
          confidence=0.95)


        #self.assertFloatEqual(brownian_motion_parameter,[1.0])    
        self.assertEqual(len(brownian_motion_parameter),1) 
        
        #Test two-trait case
        
        traits = self.SimpleTreeTraits
        tree = self.SimpleTree
        result_tree = assign_traits_to_tree(traits,tree,trait_label="Reconstruction") 
        
        true_brownian_motion_param = 5.0
        
        #E_histogram = thresholded_brownian_probability(1.0,\
        #     true_brownian_motion_param,d=0.01)
        #E_true_lower,E_true_upper = get_bounds_from_histogram(E_histogram,test_bin_edges,confidence=0.95)
         
        #set up tree with confidence intervals
        #{"A":[1.0,1.0],"E":[1.0,1.0],"F":[0.0,1.0],"D":[0.0,0.0]}
        #DndParser("((A:0.02,B:0.01)E:0.05,(C:0.01,D:0.01)F:0.05)root;")
        
        tree.getNodeMatchingName('E').upper_bound = [1.0,1.0]  
        tree.getNodeMatchingName('F').upper_bound = [1.0,2.0]
        tree.getNodeMatchingName('E').lower_bound = [-2.0,-2.0]  
        tree.getNodeMatchingName('F').lower_bound = [-1.0,0.0]
        
        brownian_motion_parameter =\
          get_brownian_motion_param_from_confidence_intervals(tree,\
          upper_bound_trait_label="upper_bound",\
          lower_bound_trait_label="lower_bound",\
          trait_label="Reconstruction",\
          confidence=0.95)


        #self.assertFloatEqual(brownian_motion_parameter,[1.0,1.0])    
        self.assertEqual(len(brownian_motion_parameter),2) 
    

if __name__ == "__main__":
    main()


