# -*- coding: utf-8 -*-
import os
import unittest
import pandas as pd
from metrics.splitscorer import SurvDiffCriterion
from classifiers.ensembl import RainForest
from classifiers.tree import  Bftree
from sklearn.model_selection import train_test_split
import pickle
import csv

snps = pd.read_csv("data/simulatedSNPdata_Survival.csv",sep=";") 
labels =  pd.read_csv("data/simulatedPatientinfo_Survival.csv",sep=";")
X_train, X_test, y_train, y_test = train_test_split(snps.iloc[:600,:200], labels.iloc[:600,:], test_size=0.33, random_state=42)
X_train.set_index("ID"), y_train.set_index("ID"), X_test.set_index("ID")

#Define criterion
svd_criterion = SurvDiffCriterion()
##Check value
RF = RainForest(max_depth=2, min_size=10, n_trees=15, n_features=200, criterion=svd_criterion)
RF.fit(X_train, y_train)
#Save model to home directory
DIR_OUT = os.getenv('HOME')
pickle.dump(RF.ttrees,open(os.path.join(DIR_OUT,"RF_trees_TEST.bin"), "wb" ))
#Predicition
pred_test = X_test.apply(RF.ensembl_predict, axis=1)
df_out = pd.DataFrame([row[1] for row in pred_test.iteritems()])
df_out = df_out.set_index(pred_test.index) #Add index
     
#Write to path
df_out.to_csv(os.path.join(DIR_OUT,"RF_TEST.csv"),header=True)
print("Output written to: ",DIR_OUT)
#Write most common snps to file
snp_counts = RF.get_tt_counts()
with open(os.path.join(DIR_OUT,"RF_TEST_SNP_COUNTS.csv"),'w') as csvfile:
    writer=csv.writer(csvfile)
    for key, value in snp_counts.items():
        writer.writerow([key] + [value])

   
                
#Reload RF for prediction
pickle.dump(RF.ttrees,open( "RF_trees.bin", "wb" ))
RF2 = RainForest(max_depth=2, min_size=10, n_trees=15, n_features=200, criterion=svd_criterion)
RF2.ttrees = pickle.load( open( "RF_trees.bin", "rb" ) )
test_pred_2 = X_test.apply(RF2.ensembl_predict, axis=1)
test_pred.to_csv(path="../test_pred.csv",index=True)

