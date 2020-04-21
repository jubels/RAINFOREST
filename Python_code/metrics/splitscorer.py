#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: tschafers

Defines a custom scoring criterion.
    - A scoring criterion must implement getSplitScore (partition) 
      and getGain (partition_left, partition_right, reference score) from
      abstract class Criterion.
    - getGain() determines if the current score improves upon the previous 
      reference score after splitting a dataset. 
"""

from abc import ABC, abstractmethod
from scipy import stats
import logging

class Criterion(ABC):
    def __init__(self):
        super().__init__()
    @abstractmethod
    def getSplitScore(self, partition):
        pass
    @abstractmethod
    def getGain(self, score_p1, score_p2, score_ref):
        pass
    

"""
Define a SurvDiff scoring criterion
   - Score: Calculates Welch's t statistic for treatment classes within partition.
   - Determines if abs. difference between left/right partition score is greater then
     reference score or the previous Node's score.

"""
class SurvDiffCriterion(Criterion):
    def __init__(self):
        super(Criterion, self).__init__()
    
    # Calculates split score
    def getSplitScore(self, partition, **kwargs):  
        tstat = 0
        try:
            a = partition.query("treatment == 0")
            b = partition.query("treatment == 1")     
            if a.shape[0] > 10 and b.shape[0] > 10:
                tstat,pval = stats.ttest_ind(a["pfs"], b["pfs"], equal_var=False)
        except Exception as exception:
            logging.exception(exception)
        return{'score':round(tstat,2)}
    
    # Determines gain
    def getGain(self, score_pos, score_neg, score_ref):
        #Check if abs diff in surv t-stats is greater than current score
        surv_diff =  abs(score_pos['score'] - score_neg['score'])
        if (surv_diff > score_ref):
            return surv_diff
        else:
            return None
       
        
    