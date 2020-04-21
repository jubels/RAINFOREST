#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 09:34:35 2019

@author: tschafers

Utility class with several helper functions (Splitting, Random sampling).

"""

from sklearn.utils import resample
from collections import Counter  
import pandas as pd


def test_split(index, value, data):
    """
    Divides a Pandas DF into feature match/mismatch subsets.
    Returns partition indexes
    """
    pos = data.index[data[index] == value].tolist()
    neg = data.index[data[index] != value].tolist()
    return pos, neg

def get_split(data, targets, tree, criterion, n_features, initial_score=0, replace=False):
    """
    Calculate the split score for a given dataset for a given node.
        - Randomly samples n features (SNP's)
        - Finds best splitting SNP and genotype to divide dataset into left/right partitions.
    Returns result dictionary
    """
    #Defaults
    features = resample(data.columns.tolist(),replace=replace, n_samples=n_features) 
    b_index = None
    b_value = 0
    b_score = initial_score
    b_groups = None
    ### ###################
    for f in features:
        #Check if feature exists within current tree.
        if tree.get_node(f) is None:
            #Retrieve all attributes for current feature.                
            attributes = data[f].unique()
            for a in data[f].unique() :
                # Ignore Heterozygous Genotypes
                if a != 1 and len(attributes) == 3:
                    idx_pos, idx_neg = test_split(f, a, data)
                    #Subset dataset
                    x1, x2 = data.loc[idx_pos], data.loc[idx_neg]
                    y1, y2 = targets.loc[idx_pos], targets.loc[idx_neg]
                    #Calculate score for each partition
                    score_pos = criterion.getSplitScore(y1)
                    score_neg = criterion.getSplitScore(y2)
                    #Get the current score and determine improvement
                    score_current = criterion.getGain(score_pos, score_neg, b_score)
                    # Update defaults if gain
                    if score_current is not None:
                        # Update
                        b_score = score_current
                        b_index = f
                        b_value = a
                        b_groups = [[x1,x2],[y1,y2]]
    return ({"index":b_index, "attr":b_value, "score":b_score, "groups":b_groups})

def random_under_sampler(targets, n_samples, replace=True, label="treatment"):
    """
    Performs random under sampling for majority class to match minority.
    Returns downsampled with adjusted class ratio.
    """
    bag = resample(targets, replace=replace, n_samples=n_samples) 
    maj_counts, min_counts = max(Counter(bag[label])),min(Counter(bag[label]))
    #Subset bag
    minority_subset = bag[bag[label] == min_counts]
    majority_subset = bag[bag[label] == maj_counts]
    #Down_sample majority class
    majority_ds = resample(majority_subset, replace=True, n_samples=len(minority_subset), random_state=123) 
    merged = pd.concat([minority_subset, majority_ds])
    #Return FIXED bag
    return(merged)
    
    
  



