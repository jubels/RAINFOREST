# gep.train = matrix with SNP data with the SNPs in columns and patient in rows
# surv.train = matrix with patients in rows. Needs to contain column 'pfs' with survival data and column 'treatment' coded as A and B
# nTrees = number of trees to be trained 
# fractionSample = which fraction of the feature should be sampled each split 
# seed = set random seed
# minSize = what is the minimum number of patients at which a node should still be split further 

RainForest = function(gep.train, surv.train,  nTrees, fractionSample, seed, minSize){
trees =  vector("list", nTrees)
set.seed(seed)
OOBpatients = vector("list", nTrees)
pick = c(0,2)

  for(x in 1:length(trees)){
    
    
    decision = matrix("YES", 3,3)
    
    
    survA = surv.train[which(surv.train$treatment == "A"),]
    survB = surv.train[which(surv.train$treatment == "B"),]
    gepA = gep.train[which(surv.train$treatment == "A"),]
    gepB = gep.train[which(surv.train$treatment == "B"),]
    
    ##### Sample out of bag

    patientsA = sample(1:nrow(gepA), nrow(gepA), replace = T)
    patientsB = sample(1:nrow(gepB), nrow(gepB), replace = T)
    OOBpatients[[x]] = c(rownames(gepA)[setdiff(1:nrow(gepA), patientsA)], rownames(gepB)[setdiff(1:nrow(gepB), patientsB)])
       
    survA = survA[patientsA,]
    survB = survB[patientsB,]
    gepA = gepA[patientsA,]
    gepB = gepB[patientsB,]
    
    useSNPs = colnames(gepA)[sample(1:ncol(gepA), fractionSample, replace = F)]
    tDiff = matrix(NA, 2,length(useSNPs))
    colnames(tDiff) = useSNPs
    for(i in 1:length(useSNPs)){
      if(length(which(gepA[,useSNPs[i]] == 0)) > 10 & length(which(gepB[,useSNPs[i]] == 0)) > 10 & length(which(gepA[,useSNPs[i]] != 0)) > 10 & length(which(gepB[,useSNPs[i]] != 0)) > 10 &
         length(which(gepA[,useSNPs[i]] == 2)) > 10 & length(which(gepB[,useSNPs[i]] ==2)) > 10 & length(which(gepA[,useSNPs[i]] != 2)) > 10 & length(which(gepB[,useSNPs[i]] != 2)) > 10){
        tDiff[1,i] = abs(t.test(survA$pfs[which(gepA[,useSNPs[i]] != 0)], survB$pfs[which(gepB[,useSNPs[i]] != 0)])[[1]] - t.test(survA$pfs[which(gepA[,useSNPs[i]] == 0)], survB$pfs[which(gepB[,useSNPs[i]] == 0)])[[1]])
        tDiff[2,i] = abs(t.test(survA$pfs[which(gepA[,useSNPs[i]] != 2)], survB$pfs[which(gepB[,useSNPs[i]] != 2)])[[1]] - t.test(survA$pfs[which(gepA[,useSNPs[i]] == 2)], survB$pfs[which(gepB[,useSNPs[i]] == 2)])[[1]])
      }
      
    }
    
    
    
    best = apply(tDiff,2,max)
    names(best) = useSNPs
    
    decision[1,1] = names(best)[which(best == max(best,na.rm = T))[1]]
    decision[1,2] = pick[which(tDiff[,decision[1,1]] == max(tDiff[,decision[1,1]]))[1]]
    model2 = t.test(survA$pfs[which(gepA[,decision[1,1]] != decision[1,2])], survB$pfs[which(gepB[,decision[1,1]] != decision[1,2])])[[1]] 
    if(model2 < 0 ){decision[1,3] = "NO"}
    if(length(which(gep.train[,decision[1,1]] == decision[1,2])) > minSize){
      survA2 = survA[which(gepA[,decision[1,1]] == decision[1,2]),]
      survB2 = survB[which(gepB[,decision[1,1]] == decision[1,2]),]
      gepA2 = gepA[which(gepA[,decision[1,1]] == decision[1,2]),]
      gepB2 = gepB[which(gepB[,decision[1,1]] == decision[1,2]),]
      useSNPs = colnames(gepA)[sample(1:ncol(gepA), fractionSample, replace = F)]
      tDiff = matrix(NA, 2,length(useSNPs))
      colnames(tDiff) = useSNPs
      
      for(i in 1:length(useSNPs)){
        if(length(which(gepA2[,useSNPs[i]] == 0)) > 10 & length(which(gepB2[,useSNPs[i]] == 0)) > 10 & length(which(gepA2[,useSNPs[i]] != 0)) > 10 & length(which(gepB2[,useSNPs[i]] != 0)) > 10 &
           length(which(gepA2[,useSNPs[i]] == 2)) > 10 & length(which(gepB2[,useSNPs[i]] ==2)) > 10 & length(which(gepA2[,useSNPs[i]] != 2)) > 10 & length(which(gepB2[,useSNPs[i]] != 2)) > 10){
          tDiff[1,i] = abs(t.test(survA2$pfs[which(gepA2[,useSNPs[i]] != 0)], survB2$pfs[which(gepB2[,useSNPs[i]] != 0)])[[1]] - t.test(survA2$pfs[which(gepA2[,useSNPs[i]] == 0)], survB2$pfs[which(gepB2[,useSNPs[i]] == 0)])[[1]])
          tDiff[2,i] = abs(t.test(survA2$pfs[which(gepA2[,useSNPs[i]] != 2)], survB2$pfs[which(gepB2[,useSNPs[i]] != 2)])[[1]] - t.test(survA2$pfs[which(gepA2[,useSNPs[i]] == 2)], survB2$pfs[which(gepB2[,useSNPs[i]] == 2)])[[1]])
        }
        
      }
      
      
      best = apply(tDiff,2,max)
      names(best) = useSNPs
      
      decision[2,1] = names(best)[which(best == max(best,na.rm = T))[1]]
      if(!is.na(decision[2,1])){
        decision[2,2] = pick[which(tDiff[,decision[2,1]] == max(tDiff[,decision[2,1]]))[1]]
        model2 = t.test(survA2$pfs[which(gepA2[,decision[2,1]] != decision[2,2])], survB2$pfs[which(gepB2[,decision[2,1]] != decision[2,2])])[[1]] 
        if(model2 < 0 ){decision[2,3] = "NO"}
      }  
    }
    if(length(which(gep.train[,decision[1,1]] != decision[1,2])) > minSize){
      survA2 = survA[which(gepA[,decision[1,1]] != decision[1,2]),]
      survB2 = survB[which(gepB[,decision[1,1]] != decision[1,2]),]
      gepA2 = gepA[which(gepA[,decision[1,1]] != decision[1,2]),]
      gepB2 = gepB[which(gepB[,decision[1,1]] != decision[1,2]),]
      useSNPs = colnames(gepA)[sample(1:ncol(gepA), fractionSample, replace = F)]
      tDiff = matrix(NA, 2,length(useSNPs))
      colnames(tDiff) = useSNPs
      
      for(i in 1:length(useSNPs)){
        if(length(which(gepA2[,useSNPs[i]] == 0)) > 10 & length(which(gepB2[,useSNPs[i]] == 0)) > 10 & length(which(gepA2[,useSNPs[i]] != 0)) > 10 & length(which(gepB2[,useSNPs[i]] != 0)) > 10 &
           length(which(gepA2[,useSNPs[i]] == 2)) > 10 & length(which(gepB2[,useSNPs[i]] ==2)) > 10 & length(which(gepA2[,useSNPs[i]] != 2)) > 10 & length(which(gepB2[,useSNPs[i]] != 2)) > 10){
          tDiff[1,i] = abs(t.test(survA2$pfs[which(gepA2[,useSNPs[i]] != 0)], survB2$pfs[which(gepB2[,useSNPs[i]] != 0)])[[1]] - t.test(survA2$pfs[which(gepA2[,useSNPs[i]] == 0)], survB2$pfs[which(gepB2[,useSNPs[i]] == 0)])[[1]])
          tDiff[2,i] = abs(t.test(survA2$pfs[which(gepA2[,useSNPs[i]] != 2)], survB2$pfs[which(gepB2[,useSNPs[i]] != 2)])[[1]] - t.test(survA2$pfs[which(gepA2[,useSNPs[i]] == 2)], survB2$pfs[which(gepB2[,useSNPs[i]] == 2)])[[1]])
        }
        
      }
      
      colnames(tDiff) = useSNPs
      
      best = apply(tDiff,2,max)
      names(best) = useSNPs
      
      decision[3,1] = names(best)[which(best == max(best,na.rm = T))[1]]
      if(!is.na(decision[3,1])){
        decision[3,2] = pick[which(tDiff[,decision[3,1]] == max(tDiff[,decision[3,1]]))[1]]
        model2 = t.test(survA2$pfs[which(gepA2[,decision[3,1]] != decision[3,2])], survB2$pfs[which(gepB2[,decision[3,1]] != decision[3,2])])[[1]] 
        if(model2 < 0 ){decision[3,3] = "NO"}
      }  
    }
    trees[[x]] = decision
    # write out trees 
    if(x %% 500 == 0){
      save(trees, file = paste("trees_", x, ".Rdata", sep = "", collapse = NULL ))
      save(OOBpatients, file = paste("OOBpatients_", x, ".Rdata", sep = "", collapse = NULL ))
      
    }
  }
  results = list(trees, OOBpatients)
  return(results)
}
