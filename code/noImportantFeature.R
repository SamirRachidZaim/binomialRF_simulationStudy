################################
##
## test RF VIS v. LASSO
##
################################
rm(list=ls())
set.seed(10)

require(parallel)
source('variableModelFunctions.R')
source('helper_functions.R')
#source('~/Desktop/shiny/test/helper.R')
require(data.table)
require(Boruta)
require(VSURF)
require(vita)
require(randomForest)
require(varSelRF)
require(devtools)
require(AUCRF)
require(EFS)



## GenerateModel
trainModels <- function(ntrees=1000, percent_features=.4, fdr.method='BY', pval_thresh =.05,
                         numSamples = 100, dim.X=100, numGenesSeeded=5, sampsize = 0.33, cbinom=cbinom){
  
  require(data.table)
  require(Boruta)
  require(VSURF)
  require(vita)
  require(perm)
  require(randomForest)
  require(compiler)
  require(varSelRF)
  require(devtools)
  require(AUCRF)
  require(EFS)
  require(binomialRF)
  require(correlbinom)
  


    load_all('~/Dropbox/Samir/binomialRF/')
    
    ########################################################################
    ##### Generate simulated data
    X = matrix(rnorm(numSamples*dim.X),ncol=dim.X)
    
    y = factor(rbinom(numSamples, size=1, prob= .3))
    X = data.frame(X)
    
    sampsize = sampsize * nrow(X)
    
    ## repeat for test data
    test.X= matrix(rnorm(numSamples*dim.X), ncol=dim.X)
    test.y = factor(rbinom(numSamples, size=1, prob= .3))
    test.X = data.frame(test.X)
    
    
    ########################################################################
    ########################################################################
    
    #### binomRF
    binom.rf <- binomialRF::binomialRF(X,y, fdr.method = fdr.method, ntrees=ntrees, percent_features = .3, user_cbinom_dist = cbinom)
    binom.rf$Significant <- binom.rf$adjSignificance < pval_thresh
    binom.rf.vars <- binom.rf$variable[which(binom.rf$adjSignificance < pval_thresh)]
    
    
    if(length(binom.rf.vars) < 2){
      binom.rf.object= rf.object
    } else {
      binom.rf.object <- randomForest(X[,binom.rf.vars], y)
    }
    
    binom.rf.size <- length(binom.rf.vars)

    ########################################################################
    ########################################################################

    return(
                Eval.Mat = data.frame(
                    ModelSize = length(binom.rf.vars),
                    Model = 'binomialRF',
                    DimX = dim.X
                )
               )

}


dim.X= 1000
sampsize <- rho <- 0.25
successprob=  calculateBinomialP(dim.X,.5)
ntrees <- trials <-  500
fdr.method='BY'

pval_thresh =0.05
fdr.method='BY'
numSamples=500

cbinom500 <- correlbinom::correlbinom(rho, successprob = successprob, trials = 500, model = 'kuk')

nsim = 50

require(doSNOW)
cl<-makeCluster(detectCores(), type = 'SOCK' ) # I have two cores
registerDoSNOW(cl)

modelList = foreach(j=1:nsim,
                    .packages=c('data.table','compiler','parallel')) %dopar%
  try(trainModels(ntrees=ntrees, percent_features=.4, fdr.method=fdr.method, pval_thresh =pval_thresh,
                  numSamples = numSamples, dim.X=dim.X, numGenesSeeded=numGenesSeeded, sampsize = rho, cbinom=cbinom500))

stopCluster(cl)
final.mat = do.call(rbind,modelList )
#
write.csv(final.mat,file=paste('../results/ntrees_analysis/simulation_', dim.X,'_dimX','.csv', sep=''))
print(summary(final.mat))

