source("trainForest.R")
source("RainForestResults.R")

#### train a RAINFOREST model 

# gep.train = matrix with SNP data with the SNPs in columns and patient in rows
# surv.train = matrix with patients in rows. Needs to contain column 'pfs' with survival data and column 'treatment' coded as A and B
# nTrees = number of trees to be trained 
# fractionSample = which fraction of the feature should be sampled each split 
# seed = set random seed
# minSize = what is the minimum number of patients at which a node should still be split further 

model = RainForest(gep.train, surv.train,  nTrees, fractionSample, seed, minSize)

#### predict on other data 
#### returns probability per sample and votes of all trees 
# gep= validation data, colnames must match training data colnames 
# model = result from RainForest function
prediction = predictRainForest(gep, model)

#### evaluate which SNPs are uesed
#### returns SNP names and count of how often they were selected 
# model = result from RainForest function
snps = snpCount(model)