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
dim.X = 100

## GenerateModel
trainModels.time <- function(ntrees=20, percent_features=.4, fdr.method='BY', pval_thresh =.05, dim.X = 10, numSamples=100){

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
    load_all("~/Dropbox/Samir/SMART_VIS/SMART_VIS/SMARTVIS")
    
    ##### Generate simulated data

    X = matrix(rnorm(numSamples*dim.X),ncol=dim.X)
    z = X%*% c(rep(1,5), rep(0,dim.X-5))
    y = factor(rbinom(numSamples, size=1, prob= 1/(1+exp(-z))))
    X = data.frame(X)


    numVars = ncol(X)

    rf.breiman.time <- system.time(rf.object <- randomForest(X, y, keep.forest = T, keep.inbag = T, ntree = ntrees))

    #### binomRF
    binom.rf.time = system.time(binom.rf <- SMARTVIS::binomialRF(X,y, fdr.method = fdr.method, ntrees=ntrees, percent_features = .4))

    ## Grow a logistic regression model
    ## and compare LASSO

    lasso.time <- system.time(lasso <- cv.glmnet(as.matrix(X), y=y, alpha=1, family="binomial"))

    ### Random Forest SOIL
    rf.soil.time <- system.time(RFSoil <- try(RF_Soil(data.frame(X), as.numeric(y)-1, data.frame(test.X),test.y, ntrees)))

    ### VSURF
    vsurf.time <- system.time(vsurf.object <- VSURF(X, y, parallel = T))

    ### Boruta
    boruta.time <- system.time(boruta.object <-  Boruta(X, y, verbose=1))
    ### Atlmann PIMP
    ### use Kolgomorov-Smirnoff test instead of
    ### specifying dist'n

    pimp.time <- system.time(pimp.test <- data.frame(PimpTest(PIMP(X,y, rf.object, S=100, parallel = T, ncores=4), para=T)$p.ks.test, stringsAsFactors = F))

    ### varSelRF
    varselrf.time <- system.time(varselrf.object <- varSelRF(xdata = X, Class = y))
    
    
    TimeTracking = data.frame(rbind(
                                    binom.rf.time,
                                    boruta.time,
                                    vsurf.time,
                                    varselrf.time,
                                    pimp.time,
                                    lasso.time))
    TimeTracking$Model = c('BinomialRF','Boruta','VSURF','VarSelRF','PIMP','Oracle')

    return(TimeTracking)
}


### Parallelize using
### Foreach package
fdr.method ='BY'

args <- commandArgs(TRUE)
ntrees = 500
nsim   = 20

dim.X =100
numSamples = 100

require(doSNOW)
cl<-makeCluster(detectCores(), type = 'SOCK' ) # I have two cores
registerDoSNOW(cl)

modelList = foreach(j=1:nsim,
                      .packages=c('data.table','compiler','parallel')) %dopar%
    try(trainModels.time(ntrees = ntrees, fdr.method = fdr.method, dim.X = dim.X, numSamples = numSamples))

stopCluster(cl)
mod = do.call(rbind, modelList)

write.csv(mod,file=paste('../results/time_', dim.X,'_features','_numSamples',numSamples,'.csv', sep=''))
