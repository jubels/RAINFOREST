predictRainForest = function(gep, RF){
  trees = RF[[1]]
  perTree = matrix(0,length(trees), nrow(gep))
  for(i in 1:length(trees)){
  
    tree = trees[[i]]
    if(any(is.na(tree))){
      tree[which(is.na(tree), arr.ind  = T)] = "YES"
    }
    if(tree[1,3] == "YES"){
      if(tree[2,3] == "YES" & tree[2,1] != "YES" & tree[2,2] != "YES" ){perTree[i,which(gep[,tree[1,1]] == tree[1,2] & gep[,tree[2,1]] == tree[2,2])] = 1}
      if(tree[2,3] == "NO"& tree[2,1] != "YES" & tree[2,2] != "YES"){perTree[i,which(gep[,tree[1,1]] == tree[1,2] & gep[,tree[2,1]] != tree[2,2])] = 1}
      if(tree[3,3] == "YES"& tree[3,1] != "YES"){perTree[i,which(gep[,tree[1,1]] !=tree[1,2] & gep[,tree[3,1]] == tree[3,2])] = 1 }
      if(tree[3,3] == "NO"& tree[3,1] != "YES"){perTree[i,which(gep[,tree[1,1]] != tree[1,2] & gep[,tree[3,1]] != tree[3,2])] = 1 }
      if(tree[2,2] == "YES"){perTree[i,which(gep[,tree[1,1]] == tree[1,2])] = 1 }
    }
    if(tree[1,3] == "NO"){
      if(tree[2,3] == "YES"& tree[2,1] != "YES"){perTree[i,which(gep[,tree[1,1]] == tree[1,2] & gep[,tree[2,1]] == tree[2,2])] = 1}
      if(tree[2,3] == "NO"& tree[2,1] != "YES"){perTree[i,which(gep[,tree[1,1]] == tree[1,2] & gep[,tree[2,1]] != tree[2,2])] = 1}
      if(tree[3,3] == "YES"& tree[3,1] != "YES" & tree[3,2] != "YES"){perTree[i,which(gep[,tree[1,1]] != tree[1,2] & gep[,tree[3,1]] == tree[3,2])] = 1 }
      if(tree[3,3] == "NO"& tree[3,1] != "YES" & tree[3,2] != "YES"){perTree[i,which(gep[,tree[1,1]] != tree[1,2] & gep[,tree[3,1]] != tree[3,2])] = 1 }
      if(tree[3,2] == "YES"){perTree[i,which(gep[,tree[1,1]] != tree[1,2])] = 1 }
    }
  }

  
  score = apply(labels,2,sum)/length(trees)
  results = list(score, perTree)
  names(results) = c("prob", "votesPerTree")
  return(results)
} 

snpCount = function(RF){
  trees = RF[[1]]
  SNPs = matrix(NA, length(trees),3)
  for(i in 1:length(trees)){
    SNPs[i,] = trees[[i]][,1]
  }
  snps = as.vector(SNPs)
  snpsUsed =  snps[which(!is.na(snps)& snps != "YES")]
  snpFreq = table(snpsUsed)
  
  snpP = as.numeric(snpFreq)
  names(snpP) = names(snpFreq)
  return(snpP)
}
  