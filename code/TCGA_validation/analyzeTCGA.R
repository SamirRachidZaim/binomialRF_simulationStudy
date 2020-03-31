##########################################################################################
##########################################################################################
#### load data and libraries
rm(list=ls())
require(binomialRF)
require(randomForest)
load('~/Dropbox/Samir/binomialRF_study/code/rhinoVirus/TCGA_all.RData')
require(TCGA2STAT)

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
set.seed(1234)

d0 <- TumorNormalMatch(brca0$dat)
name = 'BRCA'

dlist = prepareData(d0)
X = dlist$X
idx = sample(nrow(X), size=nrow(X)/1.1)

X.train = X[idx,]
X.test  = X[-idx,]

y = as.factor(dlist$y)
ytrain = y[idx]
ytest  = y[-idx]

rho = .2
ntrees= 500

cbinom = correlbinom::correlbinom(rho, successprob = calculateBinomialP(ncol(X),.5), 
                                  trials=ntrees, model='kuk', precision = ntrees)

bin.rf <- binomialRF::binomialRF(X.train,ytrain, ntrees = ntrees, fdr.threshold = 0.05, percent_features = .9,
                                 user_cbinom_dist = cbinom, sampsize = round(nrow(X)*rho))

rf.objectOriginal = randomForest(X.train, ytrain)

binomRF.vars <- bin.rf[bin.rf$adjSignificance < .2,]

if(length(binomRF.vars$variable) > 2){
  rf.object = randomForest(X.train[,binomRF.vars$variable], ytrain)
  binomRF.acc = mean(predict(rf.object, X.test) == ytest)
  cbinom2 = correlbinom::correlbinom(rho, successprob = calculateBinomialP_Interaction(ncol(X.train[,binomRF.vars$variable]),.5, K=2), 
                                     trials=ntrees, model='kuk', precision = ntrees)
  
  k.binom.RF <- k_binomialRF(X=X.train[,binomRF.vars$variable], y = ytrain, 0.05, 'BY', ntrees = 1000, 
                             percent_features = .8, K=2, cbinom_dist = cbinom2, sampsize = .25*nrow(X.train))
  
  k.binom.RF.vars <- k.binom.RF[k.binom.RF$adjSignificance < .2,]
  
  
} else{
  binomRF.acc = NA
  k.binom.RF.vars = data.frame()
  
}
org.acc = mean(predict(rf.objectOriginal, X.test) == ytest)

write.csv(binomRF.vars, file=paste(name,'CancerList.csv',sep=''), row.names=F)
write.csv(k.binom.RF.vars, file=paste(name,'CancerList_Interactions.csv',sep=''), row.names=F)

err.mat.brca = data.frame(Cancer = name,
                     RF_Accuracy = org.acc,
                     binomialRF_Accuracy=binomRF.acc,
                     BinomialRF_Genes = nrow(binomRF.vars),
                     BinomialRF_Interactions = nrow(k.binom.RF.vars),
                     NumberSamples = nrow(X)
                     )



#### 
##########################################################################################
##########################################################################################

##########################################################################################
##########################################################################################
#### run binomialRF algorithm

#### Kidney / KIPAN

d0 <- TumorNormalMatch(kipan0$dat)
name = 'Kidney'

dlist = prepareData(d0)
X = dlist$X
idx = sample(nrow(X), size=nrow(X)/1.1)

X.train = X[idx,]
X.test  = X[-idx,]

y = as.factor(dlist$y)
ytrain = y[idx]
ytest  = y[-idx]

rho = .2
ntrees= 500

cbinom = correlbinom::correlbinom(rho, successprob = calculateBinomialP(ncol(X),.5), 
                                  trials=ntrees, model='kuk', precision = ntrees)

bin.rf <- binomialRF::binomialRF(X.train,ytrain, ntrees = ntrees, fdr.threshold = 0.05, percent_features = .9,
                                 user_cbinom_dist = cbinom, sampsize = round(nrow(X)*rho))

rf.objectOriginal = randomForest(X.train, ytrain)

binomRF.vars <- bin.rf[bin.rf$adjSignificance < .2,]

if(length(binomRF.vars$variable) > 2){
  rf.object = randomForest(X.train[,binomRF.vars$variable], ytrain)
  binomRF.acc = mean(predict(rf.object, X.test) == ytest)
  cbinom2 = correlbinom::correlbinom(rho, successprob = calculateBinomialP_Interaction(ncol(X.train[,binomRF.vars$variable]),.5, K=2), 
                                     trials=ntrees, model='kuk', precision = ntrees)
  
  k.binom.RF <- k_binomialRF(X=X.train[,binomRF.vars$variable], y = ytrain, 0.05, 'BY', ntrees = 1000, 
                             percent_features = .8, K=2, cbinom_dist = cbinom2, sampsize = .25*nrow(X.train))
  
  k.binom.RF.vars <- k.binom.RF[k.binom.RF$adjSignificance < .2,]
  
  
} else{
  binomRF.acc = NA
  k.binom.RF.vars = data.frame()
  
}
org.acc = mean(predict(rf.objectOriginal, X.test) == ytest)

write.csv(binomRF.vars, file=paste(name,'CancerList.csv',sep=''), row.names=F)
write.csv(k.binom.RF.vars, file=paste(name,'CancerList_Interactions.csv',sep=''), row.names=F)

err.mat.kipan = data.frame(Cancer = name,
                          RF_Accuracy = org.acc,
                          binomialRF_Accuracy=binomRF.acc,
                          BinomialRF_Genes = nrow(binomRF.vars),
                          BinomialRF_Interactions = nrow(k.binom.RF.vars),
                          NumberSamples = nrow(X)
)



#### 
##########################################################################################
##########################################################################################

# ##########################################################################################
# ##########################################################################################
# #### run binomialRF algorithm
# 
# #### Stomach Adenocarcinoma
# 
# d0 <- TumorNormalMatch(stad0$dat)
# name = 'StomachAdeno'
# 
# set.seed(1234)
# dlist = prepareData(d0)
# X = dlist$X
# idx = sample(nrow(X), size=nrow(X)/1.1)
# 
# X.train = X[idx,]
# X.test  = X[-idx,]
# 
# y = as.factor(dlist$y)
# ytrain = y[idx]
# ytest  = y[-idx]
# 
# rho = .2
# ntrees= 500
# 
# cbinom = correlbinom::correlbinom(rho, successprob = calculateBinomialP(ncol(X),.5), 
#                                   trials=ntrees, model='kuk', precision = ntrees)
# 
# bin.rf <- binomialRF::binomialRF(X.train,ytrain, ntrees = ntrees, fdr.threshold = 0.05, percent_features = .9,
#                                  user_cbinom_dist = cbinom, sampsize = round(nrow(X)*rho))
# 
# rf.objectOriginal = randomForest(X.train, ytrain)
# 
# binomRF.vars <- bin.rf[bin.rf$adjSignificance < .2,]
# 
# if(length(binomRF.vars$variable) > 2){
#   rf.object = randomForest(X.train[,binomRF.vars$variable], ytrain)
#   binomRF.acc = mean(predict(rf.object, X.test) == ytest)
#   cbinom2 = correlbinom::correlbinom(rho, successprob = calculateBinomialP_Interaction(ncol(X.train[,binomRF.vars$variable]),.5, K=2), 
#                                      trials=ntrees, model='kuk', precision = ntrees)
#   
#   k.binom.RF <- k_binomialRF(X=X.train[,binomRF.vars$variable], y = ytrain, 0.05, 'BY', ntrees = 1000, 
#                              percent_features = .8, K=2, cbinom_dist = cbinom2, sampsize = .25*nrow(X.train))
#   
#   k.binom.RF.vars <- k.binom.RF[k.binom.RF$adjSignificance < .2,]
#   
#   
# } else{
#   binomRF.acc = NA
#   k.binom.RF.vars = data.frame()
#   
# }
# org.acc = mean(predict(rf.objectOriginal, X.test) == ytest)
# 
# write.csv(binomRF.vars, file=paste(name,'CancerList.csv',sep=''), row.names=F)
# write.csv(k.binom.RF.vars, file=paste(name,'CancerList_Interactions.csv',sep=''), row.names=F)
# 
# err.mat.stad = data.frame(Cancer = name,
#                           RF_Accuracy = org.acc,
#                           binomialRF_Accuracy=binomRF.acc,
#                           BinomialRF_Genes = nrow(binomRF.vars),
#                           BinomialRF_Interactions = nrow(k.binom.RF.vars),
#                           NumberSamples = nrow(X)
# )
# 
# 
# 
# #### 
# ##########################################################################################
# ##########################################################################################

# ##########################################################################################
# ##########################################################################################
# #### run binomialRF algorithm
# 
# #### kidney renal
# 
# d0 <- TumorNormalMatch(kirc0$dat)
# name = 'KidneyRenal'
# 
# set.seed(1234)
# dlist = prepareData(d0)
# X = dlist$X
# idx = sample(nrow(X), size=nrow(X)/1.1)
# 
# X.train = X[idx,]
# X.test  = X[-idx,]
# 
# y = as.factor(dlist$y)
# ytrain = y[idx]
# ytest  = y[-idx]
# 
# rho = .2
# ntrees= 500
# 
# cbinom = correlbinom::correlbinom(rho, successprob = calculateBinomialP(ncol(X),.5), 
#                                   trials=ntrees, model='kuk', precision = ntrees)
# 
# bin.rf <- binomialRF::binomialRF(X.train,ytrain, ntrees = ntrees, fdr.threshold = 0.05, percent_features = .9,
#                                  user_cbinom_dist = cbinom, sampsize = round(nrow(X)*rho))
# 
# rf.objectOriginal = randomForest(X.train, ytrain)
# 
# binomRF.vars <- bin.rf[bin.rf$adjSignificance < .2,]
# 
# if(length(binomRF.vars$variable) > 2){
#   rf.object = randomForest(X.train[,binomRF.vars$variable], ytrain)
#   binomRF.acc = mean(predict(rf.object, X.test) == ytest)
#   cbinom2 = correlbinom::correlbinom(rho, successprob = calculateBinomialP_Interaction(ncol(X.train[,binomRF.vars$variable]),.5, K=2), 
#                                      trials=ntrees, model='kuk', precision = ntrees)
#   
#   k.binom.RF <- k_binomialRF(X=X.train[,binomRF.vars$variable], y = ytrain, 0.05, 'BY', ntrees = 1000, 
#                              percent_features = .8, K=2, cbinom_dist = cbinom2, sampsize = .25*nrow(X.train))
#   
#   k.binom.RF.vars <- k.binom.RF[k.binom.RF$adjSignificance < .2,]
#   
#   
# } else{
#   binomRF.acc = NA
#   k.binom.RF.vars = data.frame()
#   
# }
# org.acc = mean(predict(rf.objectOriginal, X.test) == ytest)
# 
# write.csv(binomRF.vars, file=paste(name,'CancerList.csv',sep=''), row.names=F)
# write.csv(k.binom.RF.vars, file=paste(name,'CancerList_Interactions.csv',sep=''), row.names=F)
# 
# err.mat.kirc = data.frame(Cancer = name,
#                           RF_Accuracy = org.acc,
#                           binomialRF_Accuracy=binomRF.acc,
#                           BinomialRF_Genes = nrow(binomRF.vars),
#                           BinomialRF_Interactions = nrow(k.binom.RF.vars),
#                           NumberSamples = nrow(X)
# )
# 
# 
# 
# #### 
# ##########################################################################################
# ##########################################################################################

# ##########################################################################################
# ##########################################################################################
# #### run binomialRF algorithm
# 
# #### HeadNeck
# 
# d0 <- TumorNormalMatch(hnsc0$dat)
# name = 'HeadNeck'
# 
# set.seed(1234)
# dlist = prepareData(d0)
# X = dlist$X
# idx = sample(nrow(X), size=nrow(X)/1.1)
# 
# X.train = X[idx,]
# X.test  = X[-idx,]
# 
# y = as.factor(dlist$y)
# ytrain = y[idx]
# ytest  = y[-idx]
# 
# rho = .2
# ntrees= 500
# 
# cbinom = correlbinom::correlbinom(rho, successprob = calculateBinomialP(ncol(X),.5), 
#                                   trials=ntrees, model='kuk', precision = ntrees)
# 
# bin.rf <- binomialRF::binomialRF(X.train,ytrain, ntrees = ntrees, fdr.threshold = 0.05, percent_features = .9,
#                                  user_cbinom_dist = cbinom, sampsize = round(nrow(X)*rho))
# 
# rf.objectOriginal = randomForest(X.train, ytrain)
# 
# binomRF.vars <- bin.rf[bin.rf$adjSignificance < .2,]
# 
# if(length(binomRF.vars$variable) > 2){
#   rf.object = randomForest(X.train[,binomRF.vars$variable], ytrain)
#   binomRF.acc = mean(predict(rf.object, X.test) == ytest)
#   cbinom2 = correlbinom::correlbinom(rho, successprob = calculateBinomialP_Interaction(ncol(X.train[,binomRF.vars$variable]),.5, K=2), 
#                                      trials=ntrees, model='kuk', precision = ntrees)
#   
#   k.binom.RF <- k_binomialRF(X=X.train[,binomRF.vars$variable], y = ytrain, 0.05, 'BY', ntrees = 1000, 
#                              percent_features = .8, K=2, cbinom_dist = cbinom2, sampsize = .25*nrow(X.train))
#   
#   k.binom.RF.vars <- k.binom.RF[k.binom.RF$adjSignificance < .2,]
#   
#   
# } else{
#   binomRF.acc = NA
#   k.binom.RF.vars = data.frame()
#   
# }
# org.acc = mean(predict(rf.objectOriginal, X.test) == ytest)
# 
# write.csv(binomRF.vars, file=paste(name,'CancerList.csv',sep=''), row.names=F)
# write.csv(k.binom.RF.vars, file=paste(name,'CancerList_Interactions.csv',sep=''), row.names=F)
# 
# err.mat.hdnck = data.frame(Cancer = name,
#                           RF_Accuracy = org.acc,
#                           binomialRF_Accuracy=binomRF.acc,
#                           BinomialRF_Genes = nrow(binomRF.vars),
#                           BinomialRF_Interactions = nrow(k.binom.RF.vars),
#                           NumberSamples = nrow(X)
# )
# 
# 
# 
# #### 
# ##########################################################################################
# ##########################################################################################

finalSummary = rbind(err.mat.brca,
      err.mat.kipan,
      err.mat.kirc,
      err.mat.stad,
      err.mat.hdnck)

write.csv(finalSummary, file='tcgaSummary.csv', row.names = F)






