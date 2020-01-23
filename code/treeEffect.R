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
    
    trueBeta = c(rep(numGenesSeeded,5), rep(0,dim.X-numGenesSeeded))
    trueVars <- paste('X', 1:numGenesSeeded, sep='')
    
    X = matrix(rnorm(numSamples*dim.X),ncol=dim.X)
    z = X%*% c(rep(numGenesSeeded,5), rep(0,dim.X-numGenesSeeded))
    y = factor(rbinom(numSamples, size=1, prob= 1/(1+exp(-z))))
    X = data.frame(X)
    
    sampsize = sampsize * nrow(X)
    
    ## repeat for test data
    test.X= matrix(rnorm(numSamples*dim.X), ncol=dim.X)
    z.test= 1+ test.X %*% trueBeta
    pr.test= 1/(1+exp(-z.test))
    test.y = factor(rbinom(numSamples,1,pr.test))
    test.X = data.frame(test.X)
    
    
    ########################################################################
    ########################################################################
    
    #### binomRF
    binom.rf <- binomialRF::binomialRF(X,y, fdr.method = fdr.method, ntrees=ntrees, percent_features = .4, user_cbinom_dist = cbinom)
    binom.rf$Significant <- binom.rf$adjSignificance < pval_thresh
    binom.rf.vars <- binom.rf$variable[which(binom.rf$adjSignificance < pval_thresh)]
    
    
    if(length(binom.rf.vars) < 2){
      binom.rf.object= rf.object
    } else {
      binom.rf.object <- randomForest(X[,binom.rf.vars], y)
    }
    
    binom.rf.size <- length(binom.rf.vars)
    binom.rf.precision <- sum(binom.rf.vars %in% trueVars)/length(binom.rf.vars)
    binom.rf.recall <- sum(binom.rf.vars %in% trueVars)/ length(trueVars)
    
    binom.rf.test.err  <- mean(predict(binom.rf.object, test.X)!= test.y)
    
    ########################################################################
    ########################################################################

    return(
                Eval.Mat = data.frame(
                    Precision = binom.rf.precision,
                    Recall=binom.rf.recall,
                    TestError = binom.rf.test.err,
                    ModelSize = length(binom.rf.vars),
                    Model = 'binomialRF',
                    Ntree = ntrees
                )
               )

}

numGenesSeeded=5
dim.X= 100
sampsize <- rho <- 0.25
successprob=  1/dim.X
ntrees <- trials <-  500
fdr.method='BY'

pval_thresh =0.05
fdr.method='BY'
numSamples=500

cbinom100 <- correlbinom::correlbinom(rho, successprob = 1/100, trials = 100, model = 'kuk')
cbinom250 <- correlbinom::correlbinom(rho, successprob = 1/100, trials = 250, model = 'kuk')
cbinom500 <- correlbinom::correlbinom(rho, successprob = 1/100, trials = 500, model = 'kuk')
cbinom1000 <- correlbinom::correlbinom(rho, successprob = 1/100, trials = 1000, model = 'kuk')
cbinom2000 <- correlbinom::correlbinom(rho, successprob = 1/100, trials = 2000, model = 'kuk')


pdf(file='../plots/exploringnTrees.pdf', width=8, height=8)
par(mfrow=c(2,2))
plot(cbinom100, type='b', main='100 Trees'); plot(cbinom250, type='b', main='250 Trees'); 
plot(cbinom500, type='b', main='500 Trees'); plot(cbinom1000, type='b', main='1000 Trees')
dev.off()

findAlphaQuantile <- function(cbinom100, alpha=0.95){
  cmf = cumsum(cbinom100)
  first(which(cmf > alpha))
}

quantileAnalysis = cbind(
findAlphaQuantile(cbinom100, .95),
findAlphaQuantile(cbinom250, .95),
findAlphaQuantile(cbinom500, .95),
findAlphaQuantile(cbinom1000, .95),
findAlphaQuantile(cbinom2000, .95))

print(quantileAnalysis)

nsim   = 50



############################################################################################################
############################################################################################################

#### 100 trees 

ntrees = 100

require(doSNOW)
cl<-makeCluster(detectCores(), type = 'SOCK' ) # I have two cores
registerDoSNOW(cl)

modelList = foreach(j=1:nsim,
                    .packages=c('data.table','compiler','parallel')) %dopar%
  try(trainModels(ntrees=ntrees, percent_features=.4, fdr.method=fdr.method, pval_thresh =pval_thresh,
                  numSamples = numSamples, dim.X=dim.X, numGenesSeeded=numGenesSeeded, sampsize = rho, cbinom=cbinom100))

stopCluster(cl)
final.mat = do.call(rbind,modelList )
#
write.csv(final.mat,file=paste('../results/ntrees_analysis/simulation_', ntrees,'_ntrees','.csv', sep=''))
print(summary(final.mat))


############################################################################################################
############################################################################################################

############################################################################################################
############################################################################################################

#### 250 trees 

ntrees = 250

require(doSNOW)
cl<-makeCluster(detectCores(), type = 'SOCK' ) # I have two cores
registerDoSNOW(cl)

modelList = foreach(j=1:nsim,
                    .packages=c('data.table','compiler','parallel')) %dopar%
  try(trainModels(ntrees=ntrees, percent_features=.4, fdr.method=fdr.method, pval_thresh =pval_thresh,
                  numSamples = numSamples, dim.X=dim.X, numGenesSeeded=numGenesSeeded, sampsize = rho, cbinom=cbinom250))

stopCluster(cl)
final.mat = do.call(rbind,modelList )
#
write.csv(final.mat,file=paste('../results/ntrees_analysis/simulation_', ntrees,'_ntrees','.csv', sep=''))
print(summary(final.mat))


############################################################################################################
############################################################################################################



############################################################################################################
############################################################################################################

#### 500 trees 

ntrees = 500

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
write.csv(final.mat,file=paste('../results/ntrees_analysis/simulation_', ntrees,'_ntrees','.csv', sep=''))
print(summary(final.mat))


############################################################################################################
############################################################################################################

############################################################################################################
############################################################################################################

#### 1000 trees 

ntrees = 1000

require(doSNOW)
cl<-makeCluster(detectCores(), type = 'SOCK' ) # I have two cores
registerDoSNOW(cl)

modelList = foreach(j=1:nsim,
                    .packages=c('data.table','compiler','parallel')) %dopar%
  try(trainModels(ntrees=ntrees, percent_features=.4, fdr.method=fdr.method, pval_thresh =pval_thresh,
                  numSamples = numSamples, dim.X=dim.X, numGenesSeeded=numGenesSeeded, sampsize = rho, cbinom=cbinom1000))

stopCluster(cl)
final.mat = do.call(rbind,modelList )
#
write.csv(final.mat,file=paste('../results/ntrees_analysis/simulation_', ntrees,'_ntrees','.csv', sep=''))
print(summary(final.mat))


############################################################################################################
############################################################################################################

############################################################################################################
############################################################################################################

#### 2000 trees 

ntrees = 2000

require(doSNOW)
cl<-makeCluster(detectCores(), type = 'SOCK' ) # I have two cores
registerDoSNOW(cl)

modelList = foreach(j=1:nsim,
                    .packages=c('data.table','compiler','parallel')) %dopar%
  try(trainModels(ntrees=ntrees, percent_features=.4, fdr.method=fdr.method, pval_thresh =pval_thresh,
                  numSamples = numSamples, dim.X=dim.X, numGenesSeeded=numGenesSeeded, sampsize = rho, cbinom=cbinom2000))

stopCluster(cl)
final.mat = do.call(rbind,modelList )
#
write.csv(final.mat,file=paste('../results/ntrees_analysis/simulation_', ntrees,'_ntrees','.csv', sep=''))
print(summary(final.mat))


############################################################################################################
############################################################################################################

# aggregate all results 
setwd('../results/ntrees_analysis/')
t100 = fread('simulation_100_ntrees.csv')
t250 = fread('simulation_250_ntrees.csv')
t500 = fread('simulation_500_ntrees.csv')
t1000= fread('simulation_1000_ntrees.csv')
t2000= fread('simulation_2000_ntrees.csv')

res = rbind(t100, t250, t500, t1000, t2000)
res = data.table(res)
res = res[,-1]
res = melt(res, id.vars = c('Model','Ntree'))

res2 = res[, median(value), by=list(Model,Ntree, variable)]
res2 = dcast(res2, Model+Ntree ~variable)
write.csv(res2, file='ntree_summary.csv', row.names = F)




