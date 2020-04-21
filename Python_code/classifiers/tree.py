#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 16:05:19 2019

@author: tschafers

This code has partially been adapted from 
https://machinelearningmastery.com/implement-random-forest-scratch-python/

Create a new class instance of Bftree.
Bftree fits a dataset containing SNP's(features) and patients (samples)
onto a treelib class instance. 
    - Each node contains the SNP and genotype that optimatlly divides he dataset according to the scoring criterium defined 
      in splitscorer. 
    - By default, the maximum abs. difference in Welch’s T-test statistic between the left partition subset(SNP=Genotype) 
      and right partition subset(SNP!=Genotype) is used as scoring criteria. 
    - New child nodes are recursively added if abs.difference is greater then root node.   
"""



from treelib import Tree
from classifiers import utils



class Bftree:
    def __init__(self, max_depth=2, min_size=30, n_features=500, criterion=None):    
        self.max_depth = max_depth
        self.min_size = min_size
        self.n_features = n_features
        self.tree = Tree()
        self.criterion = criterion
     
    def fit(self, X, y):
        """Fits a Trainingset 

        Calculates initial best scoring feature for root node
        Recurively adds child features to root node
        
        """
        root = utils.get_split(data=X,
                               targets=y,
                               tree=self.tree,
                               criterion=self.criterion,
                               n_features=self.n_features)

        n_tag = "%s=%d (%f)" % (root['index'], root['attr'], round(root['score'],2))
        n_dict = {"attr": root['attr'], "score": root['score'], "child": 0}
        self.tree.create_node(n_tag,root['index'], data=n_dict, parent=None)
        self.split(root,1)

    def split(self, res, current_depth):
        """Splits a node into l/r child
        Adds left/right node if split partitions exist
        Split partitions are added as child nodes if score child > score root
        Halts execution if termination criteria are meet:
                    1) get_split returns empty partitions
                    2) current_depth >= max_depth
                    3) partition (l/r) size <= min_size (rows)
        Returns current tree if termination criteria meet
        """
        # Check if partition(s) exist.
        if (res['groups'] != None):
            left, right = res['groups'][0] # L/R feature partitions
            left_y, right_y = res['groups'][1] # L/R scoring partitions
        else:
            return self.tree       
        # Return tree if max_depth reached.
        if current_depth >= self.max_depth:
            return self.tree
        
        # Return tree if left feature partition < min_size
        if left.shape[0] >= self.min_size:
            # Calculate split score for left feature partition.  
            current_left = utils.get_split(data=left, 
                                           targets=left_y,
                                           tree=self.tree, 
                                           criterion=self.criterion,
                                           n_features=self.n_features)
            # Check if feature is returned or None when not improving.
            if (current_left['index'] is not None):
                # Add new left child node
                n_tag = "%s=%d(L) (%f)" % (current_left['index'], current_left['attr'], round(current_left['score'],2))
                n_dict = {"attr":current_left['attr'], "score":current_left['score'],"child":1}
                self.tree.create_node(tag=n_tag, identifier=current_left['index'],data=n_dict,parent=res['index'])
                self.split(current_left, current_depth + 1)
                
        # Return tree if feature feature partition < min_size
        if right.shape[0] >= self.min_size:
            # Calculate split score for right feature partition.  
            current_right = utils.get_split(data=right,
                                            targets=right_y,
                                            tree=self.tree,
                                            criterion=self.criterion, 
                                            n_features=self.n_features)
            # Check if feature is returned or None when not improving.
            if (current_right['index'] is not None):
                # Add new right child node
                n_tag = "%s=%d(R) (%f)" % (current_right['index'], current_right['attr'], round(current_right['score'],2))
                n_dict = {"attr":current_right['attr'], "score":current_right['score'], "child": 2}
                self.tree.create_node(tag=n_tag, identifier=current_right['index'],data=n_dict,parent=res['index'])
                self.split(current_right, current_depth + 1)

                
    def getChild(self,feature,site):
        """Helper to return the L/R child of a specific node based on a
        Node's dict.
        L_child: 1
        R_child: 2         
        """
        children = self.tree.children(feature)
        if len(children) > 0:
            for ch in children:
                if ch.data['child'] is site:
                    return ch
        else:
            return None
    
    
    def predict(self, feature, dataset):
        """"Predicts the "Class" of a feature set (Pandas series).
        Left branches will be classified as class 0 (no-benedit)
        Right branches will be classified as class 1 (benedit)
        """
        current_node = self.tree.get_node(feature)
        if( dataset[feature].item() == current_node.data['attr']):
            # Root match. Check if left child (site=1) exists.
            if self.getChild(feature,1) is not None:
                # Recursive prediction for child node.
                return(self.predict(self.getChild(feature,1).identifier,dataset))
            else:
                # Stop iteration
                return(1)      
        else:
            # Root match. Check if right child (site=2) exists.
            if self.getChild(feature,2) is not None:
                # Recursive prediction for child node.
                return(self.predict(self.getChild(feature,2).identifier,dataset))
            else:
                # Stop iteration
                return(0)

     