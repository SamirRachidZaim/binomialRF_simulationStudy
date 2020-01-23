################################
##
## test RF VIS v. LASSO
##
################################
set.seed(10)

require(parallel)
source('simulateData.R')
source('helper_functions.R')
source('variableModelFunctions.R')
#source('~/Desktop/shiny/test/helper.R')
require(data.table)
require(Boruta)
require(VSURF)
require(vita)
require(perm)
require(randomForest)
require(compiler)
require(varSelRF)
require(devtools)
ntrees = 500
percent_features = .4
percent_data=.5
pval_thresh = .05

## GenerateModel
trainModels.time <- function(ntrees=20, percent_features=.4, fdr.method='BY', pval_thresh =.05, dimX = 10, numSamples=100){

    require(data.table)
    require(Boruta)
    require(VSURF)
    require(vita)
    require(perm)
    require(randomForest)
    require(compiler)
    require(varSelRF)
    require(glmnet)
    require(devtools)
    load_all("~/Dropbox/Samir/binomialRF/")
    
    ##### Generate simulated data

    X = matrix(rnorm(numSamples*dimX),ncol=dimX)
    z = X%*% c(rep(1,5), rep(0,dimX-5))
    y = factor(rbinom(numSamples, size=1, prob= 1/(1+exp(-z))))
    X = data.frame(X)


    numVars = ncol(X)
    
    trueFeatures = colnames(X)
    
    candidateModels <- list(
      m1=sample(trueFeatures, size=floor(dimX * .3)),
      m2=sample(trueFeatures, size=floor(dimX * .3)),
      m3=sample(trueFeatures, size=floor(dimX * .3)),
      m4=sample(trueFeatures, size=floor(dimX * .5)),
      
      m5= sample(trueFeatures,  size=floor(dimX * .5)),
      m6= sample(trueFeatures,  size=floor(dimX * .5)),
      m7= sample(trueFeatures,  size=floor(dimX * .7)),
      m8= sample(trueFeatures,  size=floor(dimX * .7)),
      m9= sample(trueFeatures,  size=floor(dimX * .7)),
      m10=trueFeatures
    )

    modelAve.time <- system.time(tab <- evaluateCandidateModels(candidateModels, X,y,percent_features = 0.3, ntrees = 500))

    
    TimeTracking = data.frame(elapsed = modelAve.time['elapsed'],
                              Model = 'modelAvg')
      

    return(TimeTracking)
}


### Parallelize using
### Foreach package
fdr.method ='BY'

args <- commandArgs(TRUE)
ntrees = 500
nsim   = 500

for(dimX in c(10, 100, 1000)){
  
  numSamples = 100
  
  require(doSNOW)
  cl<-makeCluster(detectCores(), type = 'SOCK' ) # I have two cores
  registerDoSNOW(cl)
  
  modelList = foreach(j=1:nsim,
                        .packages=c('data.table','compiler','parallel')) %dopar%
      try(trainModels.time(ntrees = ntrees, fdr.method = fdr.method, dimX = dimX, numSamples = numSamples))
  
  stopCluster(cl)
  mod = do.call(rbind, modelList)
  
  write.csv(mod,file=paste('../results/ModelAvg_time_', dimX,'_features','_numSamples',numSamples,'.csv', sep=''), row.names = F)
}
