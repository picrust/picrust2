#!/usr/bin/env python
# File created on Jan 26 2012
from __future__ import division

__author__ = "Jesse RR Zaneveld"
__copyright__ = "Copyright 2015, The PICRUSt Project"
__credits__ = ["Jesse Zaneveld", "Morgan Langille"]
__license__ = "GPL"
__version__ = "2-alpha.1"
__maintainer__ = "Jesse RR Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"

from collections import defaultdict
from math import e
from copy import copy
from random import choice
from cogent.util.option_parsing import parse_command_line_parameters, make_option
from numpy.ma import masked_object
from numpy.ma import array as masked_array
from numpy import apply_along_axis,array,around,mean,maximum as numpy_max, minimum as numpy_min,\
  sqrt,sum,amax,amin,where, logical_not, argmin, histogram, add
from numpy.random import normal
from cogent.maths.stats.distribution import z_high
from cogent.maths.stats.special import ndtri
from cogent import LoadTable
from warnings import warn
from biom.table import Table


def biom_table_from_predictions(predictions, trait_ids,
                                observation_metadata={}, sample_metadata={},
                                convert_to_int=True):

    # Convert Nones to empty dicts to avoid problems with BIOM's parser
    if observation_metadata is None:
        obs_md = {}
    else:
        obs_md = observation_metadata

    if sample_metadata is None:
        sample_md = {}
    else:
        sample_md = sample_metadata

    organism_ids = predictions.keys()
    data = array(predictions.values()).T

    observation_metadata = [obs_md.get(obs_id, {}) for obs_id in trait_ids]
    sample_metadata = [sample_md.get(sample_id, {})
                       for sample_id in organism_ids]

    biom_table = Table(data, trait_ids, organism_ids,
                       observation_metadata=observation_metadata,
                       sample_metadata=sample_metadata)

    return biom_table


def assign_traits_to_tree(traits, tree, fix_bad_labels=True,
                          trait_label="Reconstruction"):
    """Assign a dict of traits to a PyCogent tree

    traits -- a dict of traits, keyed by node names
    tree -- a PyCogent phylonode object
    trait_label -- a string defining the attribute in which
    traits will be recorded.  For example, if this is set to 'Reconstruction',
    the trait will be attached to node.Reconstruction
    """
    if fix_bad_labels:
         fixed_traits = {}
         for t in traits:
             new_t = str(t)
             fixed_traits[new_t.strip('"').strip("'")] = traits[t]
         traits = fixed_traits

    for node in tree.preorder():
        node_name = node.Name.strip()
        if fix_bad_labels:
            node_name = node_name.strip().strip("'").strip('"')
        value_to_assign = traits.get(node_name, None)
        setattr(node, trait_label, value_to_assign)

    if 'root' in traits:
        setattr(tree.root(), trait_label, traits['root'])

    return tree

#Note that I needed to add an empty variance parameter to each fn,
#and a distance to the variance weighting to support a common interface

def linear_weight(d, max_d=1.0, variance=None):
    return (max_d - d) / max_d


def make_neg_exponential_weight_fn(exp_base=2.0):
    """Return function that exponentially weights exp_base by -d"""
    def neg_exponential_weight(d,variance=None):
        return exp_base**(-1*d)

    return neg_exponential_weight


def equal_weight(d, constant_weight=1.0, variance=None):
    return constant_weight


def inverse_variance_weight(d, variance, min_variance=1e-10):
    """weight by the inverse of variance
    NOTE: if variance is zero, weight by the inverse of min_variance"""
    adj_var = max(variance, min_variance)
    return 1.0 / adj_var



def fit_normal_to_confidence_interval(upper, lower, mean=None,
                                      confidence=0.95):
    """Return the mean and variance for  a normal distribution given confidence intervals
    upper -- upper bound

    lower -- lower bound

    mean -- mean, if known.  If not specified, will be set to be
    halfway between upper & lower bounds

    confidence -- the confidence for this interval, e.g.
    0.95 for the 95% confidence intervals
    """
    #Start by calculating Z-scores using the inverse normal distribution,
    #starting with the given confidence value
    z = ndtri(confidence)

    #Now we need to fit a normal distribution given the confidence
    #limits
    if mean is None:
        mean = (upper + lower) / 2.0

    #Z-score = deviation from mean/ stdev (equivalently x-mu/sigma)
    # Therefore with a little algebra we get stdev = (upper_CI - mean)/Z
    estimated_stdev = abs(upper - mean) / z
    variance = estimated_stdev**2  #i.e. sigma^2

    return mean, variance


def variance_of_weighted_mean(weights, variances, per_sample_axis=1):
    """Calculate the standard deviation of a weighted average

    weights -- a numpy array of floats representing the weights
    in the weighted average. Can be 1D or 2D (see per_sample_axis)

    variances -- a numpy array of variances for each observation (in
    the same order as variance

    axis_for_samples -- which axis to sum down to stay within samples.
    this value is ignored unless the array is 2d.

    This can be a little confusing, so consider this example:
    If this is a 1D array, it is just handled asa single sample.

    If it is a 2D array, then multiple weighted averages will be calculated
    for each column.   Lets say we set axis_for_samples equal to 1.
    (axis = axis_for_samples).

    Now, if we wanted to simultaneously
    calculate the variance of the weighted average for:

    case1: weights = [0.5,0.5], variances = [1.0,2.0]
    case2: weights = [0.1,0.9], variances = [1000.0,2000.0]

    We would supply the following numpy arrays:
    weights = [[0.5,0.5],[1.0,2.0]]
    variances = [[1.0,2.0],[1000.0,2000.0]]

    Notes:
    Formula from http://en.wikipedia.org/wiki/Weighted_mean

    The variance of the weighted mean is:
    sqrt(sum(stdev_squared * weight_squared))

    Variance is standard deviation squared, so we use that instead
    """
    if variances.shape != weights.shape:
        raise ValueError("Shape of variances (%s) doesn't match shape of weights (%s)" %(str(variances.shape),str(weights.shape)))
    if len(variances.shape) > 2:
        raise ValueError("Must supply a 1D or 2D array for variance.  Shape of supplied array was: %s" % str(variancs.shape))
    if len(variances.shape) == 1:
        per_sample_axis = None
    #Ensure weights are normalized

    if per_sample_axis is not None:
        normalized_weights = apply_along_axis(normalize_axis, per_sample_axis, weights)
    else:
        normalized_weights = normalize_axis(weights)
    squared_weights = weights**2

    return  sqrt(sum(squared_weights*variances, axis=per_sample_axis))

def normalize_axis(a):
    """Normalize a 1d numpy array by dividing each element by the total for the array"""
    return a /sum(a)

def thresholded_brownian_probability(start_state, var, d, min_val=0.0,
                                     increment=1.0, trait_prob_cutoff=0.01):
    """Calculates the probability of a given character value given a continuous, Brownian motion process

    start_state -- the starting, quantitiative value for the state (e.g. 0.36)
    var -- the variation (this is the Brownian motion parameter
      estimated during Ancestral State Reconstruction)
    d -- the branch length along which the change will occur
    min_val -- the minimum value for the character (0 for gene copy number for example)
    increment -- amount to increase gene copy number (almost always 1.0)
    trait_prob_cutoff -- the value below which you no longer care about rare possibilities

    Assumes evolution by Brownian motion (directionless, suitable for neutral evolution)

    See Felsenstein, J. "Inferring Phylogenies" p. 430
    for further discussion of threshold models of quantitative
    traits.  The usage here is as a predictor, rather than
    for phylogenetic inference.

    This approach might model things like gradual loss of a gene
    by degradation better than a sudden gain by HGT, because
    intermediate values between integer copy numbers are allowed

    Algorithm:
        Variance along a branch is calculated from the supplied variance
        For each possible state under consideration from min_val --> inf
        at increments of size 'increment', I calculate the probability of getting
        that discrete trait.

        This is calculated as the probability of having a continuous value
        in the range between that trait and the next one. So if min_val == 0 and
        increment = 1.0 (defaults) then any continuous trait value between -0.5 and
        1.5 is converted to a zero value, any value between 1.5 and 2.5 is converted to
        2, etc.

        The probability for this interval is calculated using z_high, subtracting
        the probability of the lower bound from the probability of the upper bound
        to get the prob. of being somewhere on the interval.

        The algorithm continues until the probability falls below the
        trait_prob_cutoff

    Once we have the probability of each possible change in copy number,
    we just add or subtract these from the start state to get a final value

    """
    #First we calculate the probabilities of gaining or losing
    #some number of gene copies, then apply those to the actual copy numbers
    trait_initial_value = start_state
    trait_variance =  var*d
    trait_std_dev = sqrt(trait_variance)
    mean_exp_change = 0.0 #By defn of brownian motion
    # Now we calculate the probability
    # of drawing each value of interest using the normal distribution

    #Start with a range centered around the start state

    half_width = increment / 2.0
    curr_trait_offset = 0
    i=  curr_trait_offset
    j = curr_trait_offset + half_width
    z_i = i/trait_std_dev
    z_j = j/trait_std_dev
    p = get_interval_z_prob(z_i,z_j)*2

    #The above is multiplied by two
    #because we only calculate one half of the interval
    #centered on the mean.
    result = defaultdict(float)
    result[start_state] += p

    while p > trait_prob_cutoff:
        #strategy:  calculate the high_values
        #then mirror to get the values below start state
        curr_trait_offset += increment

        i= curr_trait_offset
        j= curr_trait_offset + increment

        z_i = i/trait_std_dev
        z_j = j/trait_std_dev
        p = get_interval_z_prob(z_i,z_j)
        if p < trait_prob_cutoff:
            break

        #p = prob of deviating from the mean
        # an amount in the range between i,j
        high_value = start_state + curr_trait_offset

        result[high_value] += p
        #Now we set the low value to the same since the distribution
        #is symmetrical.
        #Note that because there is a lower threshold, multiple values
        #could accumulate, resulting in a higher P for the lower threshold
        low_value = max(min_val,start_state - curr_trait_offset)
        result[low_value] += p

    return result


def get_interval_z_prob(low_z,high_z):
    """Get the probability of a range of z scores
    """
    if low_z > high_z:
        raise ValueError(\
          "lower z value low_z must be lower than upper bound z value high_z.  Were parameters reversed?")
    #The low z represents a lower deviation from the mean,
    # and so will produce the higher probability
    high_prob = z_high(low_z)
    low_prob = z_high(high_z)
    interval_z_prob = high_prob - low_prob
    return interval_z_prob


def brownian_motion_var(d, brownian_motion_parameter):
    """Return the increase in variance due to brownian motion between two nodes
    d -- the distance between the two nodes on the tree
    brownian_motion_parameter -- variance for brownian_motion (must be calculated
    on the same tree as d).  Typically this comes from an ancestral
    state reconstruction method such as the ace function in the ape R package"""

    result =  d * array(brownian_motion_parameter)
    return result

def weighted_average_variance_prediction(tree, node,
    most_recent_reconstructed_ancestor=None, ancestral_variance=None,
    brownian_motion_parameter=None,
    trait_label="Reconstruction",
    weight_fn=linear_weight, verbose=True):
    """Predict the variance of the estimate of node traits

    tree -- a PyCogent PhyloNode tree object.   In this case,
      the tree nodes should be decorated with trait arrays stored in the
      attribute specified in trait_label.   By default this is
      node.Reconstruction.

    node_to_predict -- the Name of the node for which a reconstruction
      is desired.  Will be matched to the node.Name property on the tree.

    trait_label -- this is the node attribute where reconstructions/trait
      values are stored

    weight_fn -- this is a function that takes a distance
    on the phylogenetic tree and returns a weight.  Examples include
    linear_weight (equals distance), equal_weight (a fixed value that
    disregards distance), or neg_exponential_weight (neg. exponential
    weighting by branch length)

    brownian_motion_parameter -- the brownian motion parameter inferred
    from the ancestral state reconstruction.  This is distinct from the
    variance in ancestral state estimates due to error in reconstruction,
    and is instead more analagous to a rate of evolution.
    """
    parent_node =  node.Parent
    most_rec_recon_anc = most_recent_reconstructed_ancestor

    #Preparation
    # To handle empty (None) values, we fill unknown values
    # in the ancestor with values in the tips, and vice versa

    #STEP 1:  Infer traits, weight for most recently
    #reconstructed ancestral node
    ancestor_distance =  parent_node.distance(most_rec_recon_anc)
    ancestor_weight = weight_fn(ancestor_distance)

    #STEP 2:  Infer Parent node traits

    #Note ancestral_variance is an array with length
    #equal to number of traits

    #Variance due to error in ancestor, and evolution between last
    #reconstructed ancestor and the parent of the node to predict

    n_traits = len(ancestral_variance)

    all_weights = []
    all_variances = []

    ancestor_to_parent_variance =\
        brownian_motion_var(ancestor_distance, brownian_motion_parameter)

    all_variances.append(array(ancestral_variance) + ancestor_to_parent_variance)

    #Weight is per ancestor, but we need a 2D array, with the other
    #axis equal to the number of traits, so that we can weight each
    # trait x organism combination
    all_weights.append(array([ancestor_weight]*n_traits))

    total_weights = array([ancestor_weight]*len(ancestral_variance))

    #Need to end up with 2D arrays of weights and variances
    #where axis 1 is within-sample

    for child in parent_node.Children:
        child_traits = getattr(child,trait_label,None)
        if child_traits is None:
            #These nodes are skipped, so they
            #don't contribute to variance in the estimate
            continue

        if len(child_traits) != n_traits:
            raise ValueError("length of ancestral traits [%i] is not equal to the length of traits [%] in node %s" %(n_traits,len(child_traits),child.Name))
        distance_to_parent = parent_node.distance(child)

        organism_weight = weight_fn(distance_to_parent)
        #We've calculated a weight for the *organism*
        #Now we need to convert that to an array with length
        #equal to the number of traits
        weights_per_trait = array([organism_weight]*n_traits)

        #Generate an array of variances for the distance from the
        #sibling node to the parent node, of length equal to the number
        #of traits
        child_to_parent_var = brownian_motion_var(distance_to_parent,\
          brownian_motion_parameter)

        all_variances.append(child_to_parent_var)
        all_weights.append(weights_per_trait)

        #Need to apply the variance effects of weighting using the set of
        #weights for each trait

    # STEP 3: Predict target node given parent

    if len(all_weights[0]) == 1:
            all_weights = [w[0] for w in all_weights]
    if len(all_variances[0]) == 1:
            all_variances = [v[0] for v in all_variances]
    all_weights = array(all_weights)
    all_variances = array(all_variances)

    parent_variance_all_traits =\
      variance_of_weighted_mean(all_weights,all_variances,per_sample_axis=0)

    #This is the variance added due to evolution between the parent and the
    #predicted node
    d_node_to_parent = node.distance(parent_node)

    parent_to_node_variance = brownian_motion_var(d_node_to_parent,brownian_motion_parameter)
    #We assume variance from the parent to the node is independent of
    #variance from the variance in the parent itself.

    result = parent_variance_all_traits + parent_to_node_variance
    return result

def weighted_average_tip_prediction(tree, node,\
  most_recent_reconstructed_ancestor=None, trait_label="Reconstruction",\
  weight_fn=linear_weight, verbose=False):
    """Predict node traits, combining reconstructions with tip nodes

    tree -- a PyCogent PhyloNode tree object.   In this case,
      the tree nodes should be decorated with trait arrays stored in the
      attribute specified in trait_label.   By default this is
      node.Reconstruction.

    node_to_predict -- the Name of the node for which a reconstruction
      is desired.  Will be matched to the node.Name property on the tree.

    trait_label -- this is the node attribute where reconstructions/trait
      values are stored

    weight_fn -- this is a function that takes a distance
    on the phylogenetic tree and returns a weight.  Examples include
    linear_weight (equals distance), equal_weight (a fixed value that
    disregards distance), or neg_exponential_weight (neg. exponential
    weighting by branch length)
    """
    parent_node =  node.Parent

    most_rec_recon_anc = most_recent_reconstructed_ancestor

    #Preparation
    # To handle empty (None) values, we fill unknown values
    # in the ancestor with values in the tips, and vice versa

    #STEP 1:  Infer traits, weight for most recently
    #reconstructed ancestral node

    if most_rec_recon_anc is not None:
        anc_traits = getattr(most_rec_recon_anc,trait_label,None)
        ancestor_distance =  parent_node.distance(most_rec_recon_anc)
        ancestor_weight = weight_fn(ancestor_distance)
    else:
        anc_traits = ancestor_distance = ancestor_weight = None

    #STEP 2:  Infer Parent node traits

    if anc_traits is not None:
        prediction = array(anc_traits)*ancestor_weight
        total_weights = array([ancestor_weight]*len(prediction))
    else:
        prediction = None
        total_weights = None

    for child in parent_node.Children:

        #TODO: abstract to its own fn: get_child_trait_contribution
        #this will also allow recycling the function for weighted_only
        #prediction

        child_traits = getattr(child,trait_label,None)
        #print child.Name,":",child_traits
        if child_traits is None:
            continue

        distance_to_parent = parent_node.distance(child)
        weight = weight_fn(distance_to_parent)

        if prediction is None and total_weights is None:
            # No ancestral states available
            # Therefore, initialize weights and prediction
            # to this sibling node
            weights = array([weight]*len(child_traits))
            total_weights = weight
            prediction = array(child_traits)*weights
            continue

        prediction += array(child_traits)*weight
        total_weights += weight

    # STEP 3: Predict target node given parent

    #Without probabilites, we're left just predicting
    # the parent

    if prediction is  None:
        return None

    prediction = prediction/total_weights
    return prediction

def predict_random_neighbor(tree,nodes_to_predict,\
        trait_label="Reconstruction",use_self_in_prediction=True,\
        verbose=False):
    """Predict traits by selecting a random, annotated tip
    tree-- PhyloNode tree object, decorated with traits (see trait_label)

    trait_label -- the attribute where arrays of reconstructed traits
    are stored.  That is, if the label is 'Reconstruction', then
    node.Reconstruction or getattr(node,Reconstruction) should return
    an array of traits.

    use_self_in_prediction -- if True, allow for random prediction
    of self as one possible outcome.

    verbose -- print verbose output.
    """

    results = {}
    n_traits = None
    annotated_nodes = \
      [t for t in tree.tips() if \
       getattr(t,trait_label,None) is not None]

    for node_label in nodes_to_predict:

        if verbose:
            print "Predicting traits for node:",node_label

        node_to_predict = tree.getNodeMatchingName(node_label)
        if not use_self_in_prediction:
            possible_nodes = [t for t in annotated_nodes if \
                t.Name != node_label]
        else:
            possible_nodes = annotated_nodes

        node_to_predict = choice(possible_nodes)

        if verbose:
            print "Predicting using node:", node_to_predict.Name

        prediction = getattr(node_to_predict,trait_label)
        results[node_label] = prediction
    return results

def predict_nearest_neighbor(tree,nodes_to_predict,\
  trait_label="Reconstruction",use_self_in_prediction=True,\
  tips_only = True, verbose = False):
    """Predict node traits given labeled ancestral states

    tree -- a PyCogent phylonode object, with each node decorated with the
    attribute defined in trait label (e.g. node.Reconstruction = [0,1,1,0])

    nodes_to_predict -- a list of node names for which a trait
    prediction should be generated

    trait_label -- a string defining the attribute in which the
    trait to be reconstructed is stored.  This attribute should
    contain a numpy array of trait values (which can be set to None if
    not known)

    use_self_in_prediction -- if set to True, nodes that already
    have a trait value will be predicted using that trait value.
    If set to False, each node will ignore it's own traits when performing
    predictions, which can be useful for validation (otherwise a tree
    would need to be generated in which known nodes have their data removed)

    verbose -- output verbose debugging info

    """
    closest_annotated_node = None
    results = {}
    n_traits = None
    for node_label in nodes_to_predict:
        if verbose:
            print "Predicting traits for node:",node_label
        node_to_predict = tree.getNodeMatchingName(node_label)

        traits = getattr(node_to_predict,trait_label)

        # Do a little checking to make sure trait values either look
        # like valid numpy arrays of equal length, are not specified
        # or are set to None.

        if traits is not None:
            if n_traits is None:
                #Check that trait length is consistent
                try:
                    n_traits = len(traits)
                except TypeError:
                    raise TypeError("Node trait values must be arrays!  Couldn't call len() on %s" % traits)

            if traits and len(traits) != n_traits:
                raise ValueError(\
                  "The number of traits in the array for node %s (%i) does not match other nodes (%i)" %(\
                   node_to_predict,len(traits),n_traits))


        if not use_self_in_prediction:
            # ignore knowledge about self without modifying tree
            traits = None

        nearest_annotated_neighbor =\
          get_nearest_annotated_neighbor(tree,node_label,\
          trait_label=trait_label, tips_only = tips_only,\
          include_self = use_self_in_prediction)
        #print "NAN:", nearest_annotated_neighbor
        if nearest_annotated_neighbor is None:
            raise ValueError("Couldn't find an annotated nearest neighbor for node %s on tree" % node_label)

        results[node_label] = getattr(nearest_annotated_neighbor,trait_label)
    return results

def get_nearest_annotated_neighbor(tree,node_name,\
    trait_label="Reconstruction",tips_only= True, include_self=True):
    """Return the nearest annotated node, and its distance

    tree -- PhyloNode object, decorated with traits in the
    attribute specified in trait_label
    node -- name of the node of interest
    trait_label -- attribute where traits are stored where
    available
    tips_only -- if True, consider only extant, tip nodes
    as neighbors.  if False, allow the nearest neighbor to be
    ancestral.
    """
    n1 = tree.getNodeMatchingName(node_name)
    min_dist = 99999999999999999.0
    curr_best_match = None
    if tips_only:
        neighbors = tree.tips()
    else:
        neighbors = tree.preorder()
    #print neighbors
    for n2 in neighbors:
        traits = getattr(n2,trait_label,None)
        #print n2.Name, traits
        if not traits:
            continue
        if not include_self and n1.Name == n2.Name:
            continue
        dist = n1.distance(n2)
        if dist < min_dist:
            curr_best_match = n2
            min_dist = dist
    return curr_best_match


def calc_nearest_sequenced_taxon_index(tree,limit_to_tips = [],\
        trait_label="Reconstruction",include_self=True, verbose = True):
    """Calculate an index of the average distance to the nearest sequenced taxon on the tree"""

    distances = []
    if verbose:
        print "Finding all tree tips (may take a moment for large trees)"
    tree_tips = tree.tips()
    #TODO: there may be a way to use itertips her instead
    #if tips to examine is specified.  This would be much
    #more memory-efficient

    if limit_to_tips:
        # limit to specified tips if this value is passed
        # this is both faster and allows customized metrics for each OTU table
        # rather than just generically for all sequenced genomes + greengenes
        limit_to_tips = set(limit_to_tips)
        tips_to_examine = [t for t in tree_tips if t.Name in limit_to_tips]
    else:
        # If no set is specficied, calculate for all tips
        tips_to_examine = tree_tips

    # tips_to_examine is used as a lookup downstream
    tips_to_examine = set(tips_to_examine)

    if verbose:
        print "Calculating Nearest Sequenced Taxon Index (NTSI):"
        print "NOTE: this can be slow.  Run without -a to disable"

    #Next, build up a list of tips that are annotated with traits
    if verbose:
        print "Building a list of annotated tree tips"
    annot_tree_tip_names =\
        [t for t in tree_tips if getattr(t,trait_label,None) is not None]

    #Build up the set of tips we actually need info about:
    #The sequenced/annotated tips plus those we want to predict.
    annot_plus_to_predict = annot_tree_tip_names
    for tip in tips_to_examine:
        if tip not in annot_plus_to_predict:
            annot_plus_to_predict.append(tip)

    if verbose:
        print "Getting dists, tips for %i tips" % len(annot_plus_to_predict)


    #Armed with the set of all annotated tips, plus all tips
    #that we care to calculate NSTI for, we can restrict
    #our analysis to just these tips.

    #Restricting the space is critical for working with
    #large trees....

    dists,tree_tips =\
      tree.tipToTipDistances(annot_plus_to_predict)

    # caching here saves time with some options by preventing recalculation

    if verbose:
        print "Checking which tips are annotated with trait data..."

    not_annot_tip_indices =\
      [i for i,t in enumerate(tree_tips) if getattr(t,trait_label,None) is None]

    if verbose:
        print "Not annotated tip indices:",not_annot_tip_indices
    big_number = 1e250
    min_distances = {}

    if verbose:
        print "Finding min dists for each node..."

    if verbose and include_self:
        print "(annotated nodes of interest use themselves as the nearest neighbor)"

    for i,t in enumerate(tree_tips):
        if t not in tips_to_examine:
            continue
        #Copy to avoid messing up source tree when resetting dists
        #NOTE: need to check if this is strictly necessary?
        dist_to_curr_tip = copy(dists[i])
        #Mask non-annotated indices by setting to
        #a non-minimal value
        dist_to_curr_tip[not_annot_tip_indices] = big_number
        #if verbose:
        #    print "distances:",dist_to_curr_tip

        if not include_self:
            #Mask out self
            dist_to_curr_tip[i] = big_number

        min_dist = dist_to_curr_tip.min()
        if verbose:
            print t.Name," d(NN):",min_dist
        min_distances[t.Name]=min_dist

    # Average the nearest sequenced neighbor in each case to get a composite score
    nsti =  sum(min_distances.values())/float(len(min_distances))
    if verbose:
        print "NSTI:",nsti
    return nsti,min_distances

def get_nn_by_tree_descent(tree,node_of_interest,filter_by_property = "Reconstruction",verbose=False):
    """An alternative method for getting the NN of a node using tree descent

    The idea with this method is that if we descend the tree from a node of interest,
    then as soon as we hit a node whose length is greater than the distance from the node of interest
    to its sibling nodes, we know we have the set of nodes where the NN must reside.
    """

    #Start at a node of interest
    start_node = tree.getNodeMatchingName(node_of_interest)
    #if verbose:
    #    print "Found node:", start_node.Name

    found_NN_clade = False
    curr_base_node = start_node.Parent

    while not found_NN_clade:
        #print "curr_base_node:", curr_base_node.Name
        #print dir(curr_base_node)
        if filter_by_property:
            possible_NNs = [n for n in curr_base_node.tipChildren() if getattr(n,"Reconstruction",None) is not None]
        else:
            #consider all non-self possibilities
            possible_NNs = [n for n in curr_base_node.tipChildren() if n !=start_node]

        dists = [start_node.distance(n) for n in possible_NNs]

        curr_min = min(dists)
        if curr_base_node.Length + curr_base_node.distance(start_node) >= curr_min:
            found_NN_clade = True
            break
        curr_base_node=curr_base_node.Parent

    distance = curr_min
    nearest_neighbor = possible_NNs[argmin(dists)]

    #Find siblings
    return nearest_neighbor,distance


def get_brownian_motion_param_from_confidence_intervals(tree,\
  upper_bound_trait_label,lower_bound_trait_label,trait_label="Reconstruction",confidence=0.95):
    """Extract a Brownian motion parameter for a trait from confidence interval output
    tree -- PhyloNode tree object, decorated with confidence interval data
    upper_bound_trait_label -- the PhyloNode Property in which the upper bound for 95% CIs are stored
    upper_bound_trait_label -- the PhyloNode Property in which the upper lower bound for 95% CIs are stored

    The function infers a brownian motion parameter (sigma) for pic ASR data where only
    95% CIs are available using an approximate, quick and dirty method applicable to large trees with
    many examples (e.g. all sequenced genomes).

    Essentially, the decay in confidence from each tip with only one sequenced genome
    to its parent is used to calibrate the decrease in confidence over branch length.

    Procedure:
    - Iterate over all the tips in the tree
    - If it has a trait value, it must be a sequenced genome/characterized organism
        - Get its parent
        - If the parent has multiple sequenced children, ignore it
        - If the parent has only one sequenced child, calculate the brownian motion parameter

    """
    #variances = []  # holds arrays of variances across traits for each sequenced tip
    #distances = []  # holds a single distnace for each sequenced tip to its parent
    variances = None
    tips_with_traits = 0
    for tip in tree.iterTips():
        if getattr(tip,trait_label,None) is not None:
            tips_with_traits += 1
            # In a characterized tip
            tip_parent = tip.Parent
            more_than_one_annotated_child = False
            for c in tip_parent.Children:
                #TODO: Should this just be the weighted average variance of
                #children?
                if c.Name == tip.Name:
                    continue
                else:
                   if getattr(c,trait_label,None):
                       more_than_one_annotated_child = True
                       break
            if more_than_one_annotated_child:
                #print "Skipping node...multiple annotated children"
                continue

            dist = tip.distance(tip_parent)
            #parent_variance = tip_parent.distance(tip)
            upper_bounds = getattr(tip_parent,upper_bound_trait_label)
            lower_bounds = getattr(tip_parent,lower_bound_trait_label)
            traits = getattr(tip_parent,trait_label)
            if traits is None:
                #print "Skipping node...no ASR value for parent"
                continue
                #raise ValueError("No reconstructed trait values available for internal node.  Cannot infer brownian motion parameter from confidence intervals")
            per_trait_variances = []
            #print "LEN TRAITS:",len(traits)
            #print "upper_bounds:",upper_bounds
            #print "lower_bounds:",lower_bounds
            for i in range(len(traits)):
                #print "Finding variance for trait %i of tip %s" %(i,tip.Name)
                means, var = fit_normal_to_confidence_interval(upper_bounds[i],lower_bounds[i],\
                    mean=traits[i], confidence = 0.95)
                per_trait_variances.append(var)
            #print tip.Name,"\tPER TRAIT VARIANCES:\t",per_trait_variances
            #Empirically these are invariant by tip, so just take one example
            variances=array(per_trait_variances)
            distances=dist
            break
    if tips_with_traits == 0:
        raise ValueError("No tips have trait values annotated under label:"+trait_label)

    if variances is None:
	raise ValueError("Example variances could not be calculated for CIs. There may not be any cases of tip pairs where one has known and the other unknown trait values.")

    #now just average variances/d for all examples to get the brownian motion param
    #print "variances:", variances
    #print "distances:", distances
    #Length should be # of traits
    #brownian_motion_params = array([v/(distances[i]+eps) for i,v in enumerate(variances)])
    eps = 1e-10 # to avoid divide by zero errors
    brownian_motion_params = array(variances/(dist+eps))
    #print "brownian motion parameters:",brownian_motion_params
    #Length of each per tip brownian_motion_params entry should equal # of traits
    # Now average the per-tip arrays to get the array for the traits overall

    #average_brownian_motion_params = mean(brownian_motion_params,axis=1,keepdims=True)
    #print "average axis0:",average_brownian_motion_params
    #print "average axis1:",mean(brownian_motion_params,axis=1)
    #return average_brownian_motion_params
    return brownian_motion_params



def predict_traits_from_ancestors(tree,
                                  nodes_to_predict,\
                                  trait_label="Reconstruction",\
                                  weight_fn=linear_weight, 
                                  verbose = False,\
                                  calc_confidence_intervals=False,
                                  brownian_motion_parameter=None,\
                                  upper_bound_trait_label=None,
                                  lower_bound_trait_label=None,
                                  round_predictions=True):
    """Predict node traits given labeled ancestral states

    tree -- a PyCogent phylonode object, with each node decorated with the
    attribute defined in trait label (e.g. node.Reconstruction = [0,1,1,0])

    nodes_to_predict -- a list of node names for which a trait
    prediction should be generated

    trait_label -- a string defining the attribute in which the
    trait to be reconstructed is stored.  This attribute should
    contain a numpy array of trait values (which can be set to None if
    not known)

    verbose -- output verbose debugging info

    calc_confidence_invervals -- if true calculate confidence intervals for
    the trait prediction.  If set to true, the 95% confidence limits of
    the ASR must be provided for each internal node in the attributes
    specified by 'upper_bound_trait_label' and 'lower_bound_trait_label',
    and the Brownian motion parameter must be specified.

    brownian_motion_parameter -- a parameter describing the rate of evolution
    by a Brownian motion process

    upper_bound_trait_label -- the node attribute where 95% upper confidence
    limits for the ancestral state reconstruction can be found.  For example,
    if this is set to 'upper_limit' then node.upper_limit must return an array of floats
    (with length equal to the number of traits) for all internal nodes

    lower_bound_trait_label -- as upper_bound_trait_label, but for the lower
    confidence limit

    Output depends on whether calculate_confidence_intervals is True.
    If False:
        Returns an array of predicted trait values
    If True:
        Returns predicted trait values, a dict  of  variances for each trait
        and a dict with confidence interval information.
    """
    #if we're calculating confidence intervals, make sure we have the relevant information
    if calc_confidence_intervals:
        if upper_bound_trait_label is None or lower_bound_trait_label is None \
          or brownian_motion_parameter is None:
            err_text = "predict_traits_from_ancestors: you must specify upper_bound_trait_label, lower_bound_trait_label, and brownian_motion_parameter in order to calculate confidence intervals fro the prediction"
            raise ValueError(err_text)


    #result_tree = tree.deepcopy()
    results = {}
    variance_result = {}
    confidence_interval_results = defaultdict(dict)
    n_traits = None
    #Set up a dict to hold alredy sequenced genomes/tips with known character values
    tips_with_prior_info = {}

    #Calculate 1% progress so verbose output isn't too overwhelming.
    one_percent_progress = max(1,int((len(nodes_to_predict))*0.01))
    if verbose:
        print "Every %ith node prediction will be printed" % one_percent_progress

    #Interate through nodes, calculating trait predictions,
    #and (if requested) variances and confidence intervals

    # cache nodes to avoid tree traversals
    nodes_to_predict = set(nodes_to_predict)
    node_lookup = dict([(n.Name,n) for n in tree.tips() \
                         if n.Name in nodes_to_predict])

    print_this_node = False
    for i,node_label in enumerate(nodes_to_predict):
        if verbose:
            #Only prent every 1/100 tips predicted
            if i%one_percent_progress==0:
                print_this_node = True
            else:
                print_this_node = False

        if print_this_node:
            print "Predicting traits for node %i/%i:%s" %(i,len(nodes_to_predict),node_label)
        node_to_predict = node_lookup[node_label]

        traits = getattr(node_to_predict,trait_label)

        # Do a little checking to make sure trait values either look
        # like valid numpy arrays of equal length, are not specified
        # or are set to None.

        if traits is not None:
            if n_traits is None:
                #Check that trait length is consistent
                try:
                    n_traits = len(traits)
                except TypeError:
                    raise TypeError("Node trait values must be arrays!  Couldn't call len() on %s" % traits)

            if traits and len(traits) != n_traits:
                raise ValueError(\
                  "The number of traits in the array for node %s (%i) does not match other nodes (%i)" %(\
                   node_to_predict,len(traits),n_traits))

        #If we already know the answer (because a genome is sequenced), skip
        #the prediction step, and use the known answer.
        if traits is not None:
            # if we already know the traits (e.g. from a sequenced genome)
            # just predict those traits
            tips_with_prior_info[node_label]=traits
            #These will overwrite predictions downstream
            #predictions are still performed to roughly estimate variance
            #which will still be non-zero due to within-OTU effects

        #Find ancestral states, possibly with variance values
        if calc_confidence_intervals:
            ancestral_states,ancestral_variance =\
              get_most_recent_ancestral_states(node_to_predict,trait_label,\
              upper_bound_trait_label=upper_bound_trait_label,\
              lower_bound_trait_label=lower_bound_trait_label)

        #Find most recent ancestral node with ASR values
        most_recent_reconstructed_ancestor =\
            get_most_recent_reconstructed_ancestor(node_to_predict,trait_label)
        #print "Calc_confidence_intervals:",calc_confidence_intervals
        #print "most_recent_reconstructed_ancestor",most_recent_reconstructed_ancestor
        #Perform point estimate of trait values using weighted-average
        #prediction from last reconstructed ancestor
        prediction =\
              weighted_average_tip_prediction(tree,node_to_predict,\
              most_recent_reconstructed_ancestor =\
              most_recent_reconstructed_ancestor,\
              weight_fn = weight_fn)
        #round all predictions to whole numbers

        if round_predictions:
            prediction=around(prediction)

        results[node_label] = prediction

        #Now calculate variance of the estimate if requested
        if calc_confidence_intervals:
            variances =\
              weighted_average_variance_prediction(tree,node_to_predict,\
              weight_fn = weight_fn,\
              most_recent_reconstructed_ancestor =\
              most_recent_reconstructed_ancestor,\
              ancestral_variance=ancestral_variance,\
              brownian_motion_parameter=brownian_motion_parameter)

            variance_result[node_label] = {"variance":variances.tolist()}
           
            lower_95_CI,upper_95_CI = calc_confidence_interval_95(prediction,
                                                                  variances,
                                                                  round_CI=round_predictions)
            confidence_interval_results[node_label]['lower_CI']=lower_95_CI
            confidence_interval_results[node_label]['upper_CI']=upper_95_CI


        if print_this_node:
            n_traits_to_print = min(len(prediction),50)
            print "First %i trait predictions:%s" %(n_traits_to_print,\
              ','.join(map(str,list(prediction[:n_traits_to_print]))))

            if calc_confidence_intervals:
                print "First %i trait variances:%s" %(n_traits_to_print,\
                  ','.join(map(str,list(variances[:n_traits_to_print]))))
                print "Lower 95% confidence interval:",lower_95_CI[:n_traits_to_print]
                print "Upper 95% confidence interval:",upper_95_CI[:n_traits_to_print]

    #Overwrite known results from the dict of known results
    results.update(tips_with_prior_info)

    if calc_confidence_intervals:
        return results,variance_result, confidence_interval_results
    else:
        return results

def calc_confidence_interval_95(predictions,variances,round_CI=True,\
        min_val=0.0,max_val=None):
    """Calc the 95% confidence interval given predictions and variances"""
    stdev = sqrt(variances)
    pred = predictions
    CI_95 =  1.96*stdev
    lower_95_CI = pred - CI_95
    upper_95_CI = pred + CI_95
    if round_CI:
        lower_95_CI = around(lower_95_CI)
        upper_95_CI = around(upper_95_CI)
    if min_val is not None:
        lower_95_CI = numpy_max(min_val,lower_95_CI)
        upper_95_CI = numpy_max(min_val,upper_95_CI)

    if max_val is not None:
        lower_95_CI = numpy_min(max_val,lower_95_CI)
        upper_95_CI = numpy_min(max_val,upper_95_CI)

    return lower_95_CI,upper_95_CI




def fill_unknown_traits(traits_to_update, new_traits):
    """Returns  traits_to_update, with all None values replaced with the corresponding index in new_traits

    traits_to_update -- a numpy array, with None values representing missing data
    new_traits -- a numpy array, with None values representing missing data

    returns an array that is traits_to_update, with all None values replaced with
    the corresponding values in new_traits.

    This is useful when some, but not all traits of an organism are known.
    For example, PCR of specific genes may have fixed some traits, but you would
    like to predict the rest.

    Note -- in the special case where traits_to_update is the value None,
    new_traits will be returned.  This is done so that we don't have to guess
    or specify beforehand how many traits there will be when filling in tip nodes
    """
    if traits_to_update is None:
        return array(new_traits)
    masked_traits_to_update = masked_object(traits_to_update,None)

    result = where(masked_traits_to_update.mask ==  True,\
      new_traits,masked_traits_to_update)

    return array(result)

def get_most_recent_ancestral_states(node,trait_label,\
    upper_bound_trait_label=None, lower_bound_trait_label=None):
    """Traverse ancestors of node until a reconstructed value is found

    node -- a PhyloNode object
    trait_label -- the trait attribute corresponding to


    """
    for ancestor in node.ancestors():
        trait = getattr(ancestor,trait_label)
        if trait is not None:
            if not upper_bound_trait_label and not lower_bound_trait_label:
                return trait
            else:
                upper_bound = list(getattr(ancestor,upper_bound_trait_label,None))
                lower_bound = list(getattr(ancestor,lower_bound_trait_label,None))
                ancestral_variances = []
                for i in range(len(upper_bound)):
                    mu, var = \
                      fit_normal_to_confidence_interval(upper_bound[i],\
                      lower_bound[i],mean=trait[i], confidence = 0.95)
                    ancestral_variances.append(var)
                return trait, array(ancestral_variances)


    # If we get through all ancestors, and no traits are found,
    # then there are no most recent ancestral states
    return None


def get_most_recent_reconstructed_ancestor(node,trait_label="Reconstruction"):
    """Traverse ancestors of node until a reconstructed value is found

    node -- a PhyloNode object
    trait_label -- the trait attribute corresponding to


    """
    for ancestor in node.ancestors():
        trait = getattr(ancestor,trait_label)
        if trait is not None:
            return ancestor
    # If we get through all ancestors, and no traits are found,
    # then there are no most recent reconstructed ancestors
    return None

def update_trait_dict_from_file(table_file, header = [],input_sep="\t"):
    """Update a trait dictionary from a table file

    table_file --  File name of a trait table.

    The first line should be a header line, with column headers equal to trait
    (e.g. gene family) names, while the row headers should be organism
    ids that match the tree.

    trait_dict -- a dictionary of traits, keyed by organism.
    Items in trait dict will be overwritten if present.
    """
    #First line should be headers
    table=LoadTable(filename=table_file,header=True,sep=input_sep)

    #do some extra stuff to match columns if a header is provided
    if header:
        #error checking to make sure traits in ASR table are a subset of traits in genome table
        if set(header) != set(table.Header[1:]):
            if set(header).issubset(set(table.Header[1:])):
                diff_traits = set(table.Header[1:]).difference(set(header))
                warn("Missing traits in given ASR table with labels:{0}. Predictions will not be produced for these traits.".format(list(diff_traits)))
            else:
                raise RuntimeError("Given ASR trait table contains one or more traits that do not exist in given genome trait table. Predictions can not be made.")

        #Note: keep the first column heading at the beginning not sorted (this is the name for the row ids
        sorted_header=[table.Header[0]]
        sorted_header.extend(header)
        table = table.getColumns(sorted_header)

    traits = {}
    for fields in table:
        try:
            traits[fields[0]] = map(float,fields[1:])
        except ValueError:
            err_str =\
                    "Could not convert trait table fields:'%s' to float" %(fields[1:])
            raise ValueError(err_str)

    return table.Header[1:],traits

def normal_product_monte_carlo(mean1,variance1,mean2,variance2,confidence =0.95, n_trials = 5000):
    """Estimate the lower & upper confidence limits for the product of two normal distributions

    mean1 -- mean for the first normal distribution
    variance1 -- variance for the first normal distribution
    mean2 -- mean for the second normal distribution
    variance2 -- variance for the second normal distribution

    confidence -- the desired confidence interval
    n_trials -- number of monte carlo trials to use to simulate distribution
    """

    #numpy.random.normal takes stdev rather than variance, so convert variance to stdev
    stdev1 = variance1**0.5
    stdev2 = variance2**0.5

    sim_dist1 = normal(mean1,stdev1,n_trials)
    sim_dist2 = normal(mean2,stdev2,n_trials)

    #Note that multiplying arrays in numpy is element-wise not matrix multiplication
    sim_combined_dist = sim_dist1 * sim_dist2
    hist,bin_edges = histogram(sim_combined_dist,bins=n_trials)
    lower,upper = get_bounds_from_histogram(hist,bin_edges,confidence)
    return lower,upper

def get_bounds_from_histogram(hist,bin_edges,confidence=0.95):
    """Return the bins corresponding to upper and lower confidence limits

    hist  -- an array of histogram values
    bin_edges -- an array of bin edges (length equal to hist + 1)
    confidence -- the percent confidence required.  For example 0.95
    will result in the 5% and 95% confidence limits being returned

    NOTE:  since data are binned, the 5% or 95% confidence intervals
    will not be exact.  We choose to return the more conservative
    bin, so actual confidence will be <= 5% or >= 95%

    """

    total = sum(hist)
    normed_hist = hist/total
    #print "normed_hist:",normed_hist

    #The frequency of draws allowed to fall on *each* side
    #of the confidence interval is half the total
    confidence_each_tail = (1.0 - confidence)/2.0

    cum_sum = add.accumulate(normed_hist)
    #print "cum_sum:",cum_sum
    upper_indices = where(cum_sum > 1.0-confidence_each_tail)
    lower_indices = where(cum_sum < confidence_each_tail)
    #print "lower_indices:",lower_indices
    #print "upper_indices:",upper_indices
    lower_index = amax(lower_indices)
    upper_index = amin(upper_indices)

    lower_bin_index = lower_index + 1
    upper_bin_index = upper_index + 1

    upper_bound = bin_edges[upper_bin_index]
    lower_bound = bin_edges[lower_bin_index]

    return lower_bound,upper_bound
