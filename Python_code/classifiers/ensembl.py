#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 16:04:10 2019

@author: tschafers

This code has partially been adapted from 
https://machinelearningmastery.com/implement-random-forest-scratch-python/


Create multiple instances of class bftree that can be used for benefit prediction.

    - Predefined number of bftree (n_trees)
    - Random sampling is done for pre defined bag size (samples) and features (SNP's)
    - The majority class is down-samples for each bag.
    - The final prediction contain the majority vote for each trained tree.
"""
from collections import Counter  
import itertools
from classifiers.tree import Bftree
from classifiers.utils import random_under_sampler

import logging


logging.basicConfig(level=logging.DEBUG)

class RainForest():
    def __init__(self, max_depth=2, min_size=500, n_trees=1, bag_size=1000, n_features=1000, threshold=0.42, criterion=None):
        self.max_depth = max_depth
        self.min_size = min_size
        self.n_trees = n_trees
        self.bag_size = bag_size
        self.ttrees = []
        self.n_features = n_features
        self.threshold = threshold
        self.criterion = criterion
    #Training function for RF classifier
    def fit(self,X,y,replace=True):
        """
        Fits a bagged Xset onto multiple instances of bftree
        """
        try:
            # verify 
            assert len(X) == len(y)
            for i in range(self.n_trees):
                logging.info("Iteration:#%d" % (i+1))
                # Performs random under sampling for each bag
                train_y = random_under_sampler(y,self.bag_size)
                train_x = X[X.index.isin(train_y.index)]
                # Create bftree instance and fit bagged Xset.
                bf_tree = Bftree(max_depth=self.max_depth, 
                                 min_size=self.min_size,
                                 n_features=self.n_features,
                                 criterion=self.criterion)
                # Train tree
                bf_tree.fit(train_x,train_y)
                self.ttrees.append(bf_tree)   
                bf_tree.tree.show()
        except AssertionError as error:
            logging.exception(error)
        except Exception as exception:
            logging.exception(exception)
            
    def predictFinalClass(self,counts):
            p_pos = Counter(counts)[1]/len(self.ttrees) 
            p_neg = Counter(counts)[0]/len(self.ttrees) 
            if (p_pos > self.threshold or p_neg > self.threshold) and (p_pos != p_neg):
                if p_pos > p_neg:
                    return 1, round(p_pos,2), round(p_neg,2)
                else:
                    return 0, round(p_pos,2), round(p_neg,2)
            else:
                return "NA", p_pos, p_neg
        
    
    def ensembl_predict(self,X):
        """
        Performs benefit prediction for a data series (Pandas series)
        Returns posterior class probabilites, final class
        """
        temp = []
        for bft in self.ttrees:
            pred_class = bft.predict(feature=bft.tree.root,dataset=X)
            temp.append(pred_class)
        #Get most informatve SNPS      
        #Evaluate
        cl_final = None, 0, 0
        if(len(self.ttrees) > 0):
            cl_final = self.predictFinalClass(temp)
        return {"class":cl_final[0], "p_benefit":cl_final[1],"p_no_benefit": cl_final[2], "trained_trees": len(self.ttrees)}
    
    def get_tt_counts(self):
        snp_dicts = [list(tree.tree.nodes.keys()) for tree in self.ttrees]
        snps_flatten = list(itertools.chain.from_iterable(snp_dicts))
        return(Counter(snps_flatten))
