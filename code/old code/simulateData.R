################################################
################################################
#####
##### simulate binomial data
#####
##### Simulates datasets with 5 relevant and 5
##### irrelevant coefficients/variables in one
##### dataset with no correlation and one with
##### extremeley high correlation
#####
##### Author: Samir Rachid Zaim
##### Date : 10/11/2018
################################################
################################################

require(magic)
require(Matrix)
require(MASS)

set.seed(666)

simUncorrData <- function(pathwaySize=250, numPathways=10, numSamples=10){
  
  X = matrix(rnorm(numPathways*pathwaySize*numSamples), ncol=numSamples)           # some continuous variables
  trueBeta= rep(c(rep(1,2), rep(0,pathwaySize-2)), numPathways)
  
  z = t(X) %*% trueBeta        # linear combination with a bias
  pr = 1/(1+exp(-z))         # pass through an inv-logit function
  y = rbinom(numSamples,1,pr)      # bernoulli response variable
  
  return(list(X=X, y=y,Betas =trueBeta))
}


generateBlockCovariance <- function(pathwaySize, correlation=.9, numPathways=10){
  Sigma_list=  lapply(1:numPathways, function(x) matrix(correlation, nrow=pathwaySize, ncol=pathwaySize) +diag(pathwaySize)*(1-correlation))
  block_diag <- Reduce(adiag, Sigma_list)
  return(block_diag)
}

simCorrData <- function(pathwaySize=250, correlation=.9, numPathways=10, numSamples=10){
  
  require(magic)
  require(MASS)
  
  Sigma = generateBlockCovariance(pathwaySize = pathwaySize, correlation = correlation, numPathways=numPathways)
  Sigma=Matrix(Sigma, sparse = T)
  X = t(mvrnorm(n=numSamples, mu=rep(0,pathwaySize *numPathways), Sigma=Sigma))
  
  trueBeta= rep(c(rep(1,2), rep(0,pathwaySize-2)), numPathways)
  
  z =  t(X) %*% trueBeta      # linear combination with a bias
  pr = 1/(1+exp(-z))         # pass through an inv-logit function
  y = rbinom(numSamples,1,pr)      # bernoulli response variable
  
  
  return(list(X=X, y=y,Betas =trueBeta))
}

