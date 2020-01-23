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
                         numSamples = 100, dim.X=100, numGenesSeeded=5, sampsize = 0.33){
  
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
  



    ########################################################################
    ##### Generate simulated data
    X = matrix(rnorm(numSamples*dim.X),ncol=dim.X)
    #z = X%*% c(rep(numGenesSeeded,5), rep(0,dim.X-numGenesSeeded))
    
    z = 3* X[,1] + 3* X[,2] + 3* X[,3] + 3* X[,4] +3 * (X[,1] * X[,2]) + 3 * (X[,1] * X[,3]) + 3 * (X[,1] * X[,4]) + 
      3 * (X[,2] * X[,3]) +  3 * (X[,2] * X[,4]) +  3 * (X[,3] * X[,4])
    
    y = factor(rbinom(numSamples, size=1, prob= 1/(1+exp(-z))))
    X = data.frame(X)
    
    sampsize = .25 * nrow(X)
    
    ## trueVariable vector
    trueInterVars = c('X1_X2','X1_X3','X1_X4','X2_X3','X2_X4','X3_X4')
    #####
    ########################################################################
    ########################################################################
    
    ########################################################################
    ########################################################################
    
    binomRF = binomialRF(X,y, .05, 'BY', 1000, percent_features = .4, user_cbinom_dist = cbinom2)
    binomRFvars <- binomRF$variable[binomRF$adjSignificance < .05]
    
    #### binomRF
    binom.rf <- binomialRF::k_binomialRF(X[,binomRFvars],y, fdr.method = 'bonferroni', ntrees=ntrees, percent_features = .6, cbinom_dist  = cbinom)
    binom.rf$Significant <- binom.rf$adjSignificance < .05
    binom.rf.vars <- binom.rf$Interaction[which(binom.rf$Significant ==T)]
    
    
    binomRF.varlist <- strsplit(binom.rf.vars, ' | ')
    binomRF.varlist.idx <- sapply(1:length(binomRF.varlist), function(x) binomRF.varlist[[x]][1] == binomRF.varlist[[x]][3])
    idx = which(binomRF.varlist.idx !=T)
    
    binomRF_TP <- sapply(idx, function(x) paste(binomRF.varlist[[x]][c(1,3)],collapse='_') %in% trueInterVars |paste(binomRF.varlist[[x]][c(3,1)],collapse='_') %in% trueInterVars)
    
    binom.rf.vars = stringr::str_replace(binom.rf.vars, ' \\| ', '_')

    binom.rf.size <- length(idx)
    binom.rf.precision <- sum(binomRF_TP)/length(idx)
    binom.rf.recall <- sum(binomRF_TP)/ length(trueInterVars)
    
    return(
                Eval.Mat = data.frame(
                   Precision = binom.rf.precision,
                   Recall=binom.rf.recall,
                   ModelSize = binom.rf.size,
                    Model = 'binomialRF'
                )
               )

}

dim.X= 30
sampsize <- rho <- 0.25
ntrees <- trials <-  1000
fdr.method='BY'

pval_thresh =0.05
fdr.method='BY'
numSamples=1000

successProb = calculateBinomialP_Interaction(L=dim.X, percent_features = .5, K=2)

cbinom2 = correlbinom(.25, successprob = 1/dim.X, trials = ntrees,model = 'kuk')


cbinom <- correlbinom::correlbinom(rho, successprob = successProb, trials = trials, model = 'kuk')
# write.csv(cbinom, '../cbinom_distributions/interaction_cbinom.csv')
# plot(cbinom)

### Parallelize using
### Foreach package
nsim   = 50

require(doSNOW)
cl<-makeCluster(detectCores(), type = 'SOCK' ) # I have two cores
registerDoSNOW(cl)

modelList = foreach(j=1:nsim,
                    .packages=c('data.table','compiler','parallel')) %dopar%
  try(trainModels(ntrees=ntrees, percent_features=.4, fdr.method=fdr.method, pval_thresh =pval_thresh,
                  numSamples = numSamples, dim.X=dim.X, numGenesSeeded=numGenesSeeded, sampsize = rho))

stopCluster(cl)
final.mat = do.call(rbind,modelList )
#
write.csv(final.mat,file='../results/binomialRF1000_interaction.csv')


