################################
##
## test RF VIS v. LASSO
##
################################
set.seed(2000)
require(parallel)
source('simulateData.R')
source('variableModelFunctions.R')
source('helper_functions.R')
source('../../binomialRF/R/k.binomialRF.R')
#source('~/Desktop/shiny/test/helper.R')
require(data.table)
require(randomForest)
require(devtools)
#    load_all("~/Dropbox/Samir/SMART_VIS/SMART_VIS/SMARTVIS")
load_all('~/Dropbox/Samir/binomialRF/')


numSamples=100
dimX=10
BetaSignal=3
##### Generate simulated data
### Generate multivariate normal data in R10
X = matrix(rnorm(numSamples*dimX), ncol=dimX)

### let half of the coefficients be 0, the other be 10
trueBeta= c(rep(BetaSignal,5), rep(0,dimX-5))

### do logistic transform and generate the labels
z =  X %*% trueBeta    
pr = 1/(1+exp(-z))        
y = factor(rbinom(numSamples,1,pr))

candidateModels <- list(
      m1=c(paste('X',c(1,4,5,6),sep='')),
      m2=c(paste('X',c(2,3,5,6,9),sep='')),
      m3=c(paste("X", c(2,5:10),sep='')),
      m4= c(paste("X", 3:9,sep='')),
      m5= c(paste("X", c(1,4,5,7,9),sep='')),
      m6= c(paste("X", c(2,5,6,8,9),sep='')),
      m7= c(paste("X", c(5,7,9,10),sep='')),
      m8= c(paste("X", c(3,4,5,8,10),sep='')),
      m9= c(paste("X", c(1,4,5,7,9,10),sep='')),
      m10= c(paste("X",1:10,sep=''))
      )

tab = evaluateCandidateModels(candidateModels, X,y,percent_features = 0.3, ntrees = 2000)
tab



