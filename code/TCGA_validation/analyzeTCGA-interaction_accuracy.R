##########################################################################################
##########################################################################################
#### load data and libraries
rm(list=ls())
require(binomialRF)
require(randomForest)
load('~/Dropbox/Samir/binomialRF_study/code/rhinoVirus/TCGA_all.RData')
require(TCGA2STAT)
setwd("~/Dropbox/Samir/binomialRF_study/code/TCGA_validation")

#### 
##########################################################################################
##########################################################################################

##########################################################################################
##########################################################################################

##########################################################################################
##########################################################################################
# Create a label column

prepareData <- function(d0){
  tumor_group <- rep("tumor", nrow(d0$primary.tumor))
  normal_group <- rep("normal", nrow(d0$normal))
  
  tumor <- data.frame(d0$primary.tumor)
  normal <- data.frame(d0$normal)
  
  tumor <- data.frame(t(tumor))
  normal <- data.frame(t(normal))
  
  tumor$group <- 'tumor'
  normal$group <- 'normal'
  
  # Combine the tumor and normal datasets 
  d <- rbind(tumor, normal)
  dim(d)
  
  # Remove columns with all zero values: 198 columns have all zero values
  d <- d[, colSums(d != 0) > 10]
  dim(d)
  
  X = d[,!names(d) %in% 'group']
  y = d$group
  
  return(list(X=X,y=y))
}

  
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################





##########################################################################################
##########################################################################################
#### run binomialRF algorithm

#### BRCA

get_accuracy <- function(d0, name, interactions, featurelist, interactions.mat, seed=1234){
      set.seed(seed)
  
      dlist = prepareData(d0)
      X = dlist$X
      idx = sample(nrow(X), size=nrow(X)/1.5)
      
      X.train = X[idx,]
      X.test  = X[-idx,]
      
      y = as.factor(dlist$y)
      ytrain = y[idx]
      ytest  = y[-idx]
      
      ntrees= 500
      
      
      interaction_set <- data.frame(matrix(rnorm(nrow(X.train)*nrow(interactions)), ncol=nrow(interactions)))
      interaction_set.test <- data.frame(matrix(rnorm(nrow(X.test)*nrow(interactions)), ncol=nrow(interactions)))
      
      for( i in 1:nrow(interactions)){
        
        inter.mat = X.train[, interactions[i,]]
        interaction_set[,i] <- inter.mat[,1] * inter.mat[,2]
        
       
        inter.mat.test = X.test[, interactions[i,]]
        interaction_set.test[,i] <- inter.mat.test[,1] * inter.mat.test[,2]
        
      }
      
      names(interaction_set) <- interactions.mat$Interaction
      names(interaction_set.test) <- interactions.mat$Interaction
      
      
      rf.objectOriginal = randomForest(X.train, ytrain)
      rf.objectMainEffects =  randomForest(X.train[, featurelist], ytrain)
      rf.objectInteraction = randomForest(interaction_set, ytrain)
      rf.objectMainEffectsAndInteractions = randomForest(cbind(X.train[, featurelist], interaction_set), ytrain)
      
      rf_pred = predict(rf.objectOriginal, X.test)
      binomRF.main = predict(rf.objectMainEffects, X.test[, featurelist])
      binomRF.inter= predict(rf.objectInteraction, interaction_set.test)
      binomRF.main_inter= predict(rf.objectMainEffectsAndInteractions, cbind( X.test[, featurelist],interaction_set.test))
      
      prec_rec <- function(pred, label){
        tp = sum((pred=='tumor') == (label=='tumor'))
        fp = sum((pred=='tumor') == (label=='normal'))
        fn = sum((pred=='normal') == (label=='tumor'))
        tn = sum((pred=='normal') == (label=='normal'))
        
        prec= tp/(tp + fp)
        rec = tp/(tp + fn)
        return(c(prec,rec))
      }
      
      
      scores = do.call(rbind, lapply(list(rf_pred,
                                          binomRF.main,
                                          binomRF.inter,
                                          binomRF.main_inter
                                          ), 
                                     function(x) prec_rec(x, ytest)))
      
      scores = round(scores,2)
      
      brca_accs = data.frame(
        Model = c('RF','RF_genes','RF_interactions'),
        NumGenes = c(ncol(X.train), length(featurelist), nrow(interactions.mat) ),
        Precision = scores[1:3,1],
        Recall = scores[1:3,2]
      )
      

      return(brca_accs)
}



## brca
featurelist <- read.csv('results/BRCACancerList.csv', stringsAsFactors = F)
featurelist.canc <- featurelist$variable

interactions.mat.canc <- read.csv('results/BRCACancerList_Interactions.csv', stringsAsFactors = F)
interactions.canc <- do.call(rbind,lapply(1:nrow(interactions.mat.canc), function(x) unlist(strsplit(interactions.mat.canc$Interaction[x],' \\| '))))

## kidney
featurelist <- read.csv('results/kidneyCancerList.csv', stringsAsFactors = F)
featurelist.kidn <- featurelist$variable

interactions.mat.kidn <- read.csv('results/KidneyCancerList_Interactions.csv', stringsAsFactors = F)
interactions.kidn <- do.call(rbind,lapply(1:nrow(interactions.mat.kidn), function(x) unlist(strsplit(interactions.mat.kidn$Interaction[x],' \\| '))))


d0 <- TumorNormalMatch(brca0$dat)
d1 <- TumorNormalMatch(kipan0$dat)



#### get results 
require(parallel)
require(data.table)
require(ggplot2)
mat1 = mclapply(1:10, function(x) get_accuracy(d0, 'BRCA',   interactions = interactions.canc, featurelist = featurelist.canc, interactions.mat = interactions.mat.canc, seed=x))
mat2 = mclapply(1:10, function(x) get_accuracy(d1, 'Kidney', interactions = interactions.kidn, featurelist = featurelist.kidn, interactions.mat = interactions.mat.kidn, seed=x))

test1 = melt(data.table(do.call(rbind,mat1)))
fin.mat1 = dcast(test1[, mean(value), by=list(Model, variable)], formula = Model ~ variable)

test2 = melt(data.table(do.call(rbind,mat2)))
fin.mat2 = dcast(test2[, mean(value), by=list(Model, variable)], formula = Model ~ variable)

plot(x = test1$value[test1$variable=='Recall'], y = test1$value[test1$variable=='Precision'], xlim=c(0,1), ylim=c(0,1), col='red', pch=2, size=3,cex =1.4)
points(x = test1$value[test1$variable=='Recall'], y = test1$value[test1$variable=='Precision'], xlim=c(0,1), ylim=c(0,1), col='orange', pch=3, size=3,cex =1.4)
points(x = test1$value[test1$variable=='Recall'], y = test1$value[test1$variable=='Precision'], xlim=c(0,1), ylim=c(0,1), col='green', pch=4, size=3,cex =1.4)



write.csv(fin.mat1, file='results/interaction_accuracies.csv', row.names = F )
write.csv(fin.mat2, file='results/interaction_accuracies.csv', row.names = F )
