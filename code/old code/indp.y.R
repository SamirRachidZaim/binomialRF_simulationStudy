
## testing binomialRF when features not significant 
require(devtools)
load_all("~/Dropbox/Samir/SMART_VIS/SMART_VIS/SMARTVIS")
source('~/Dropbox/Samir/SMART_VIS/code/simulateData.R')

run_indp.sim <- function(){
  ### Generate multivariate normal data in R10
  X = matrix(rnorm(1000), ncol=10)
  y = rbinom(100,1,.4)
  binom.rf <- SMARTVIS::binomialRF(X,factor(y), fdr.threshold = .05,
                                   ntrees = 1000,percent_features = .6,
                                   fdr.method = 'BY')
  return(binom.rf)
}

nsim=1000
cv.bigMat = data.table(do.call(rbind,lapply(1:nsim, function(i) run_indp.sim())))
setkey(cv.bigMat, Variable)
cv.bigMat=cv.bigMat[, list(cv.Freq =mean(Frequency), cv.Pvalue =median(Pvalue), cv.AdjPvalue=median(AdjPvalue), Avg.TimesSignificant=mean(Significant)), by=Variable]
print(cv.bigMat)
write.csv(cv.bigMat, file='../results/indp_sim.csv')



run_dpnd.sim <- function(){
  
  X = matrix(rnorm(1000), ncol=10)           # some continuous variables
  trueBeta= c(rep(10,5), rep(0,5))
  
  z = X %*% trueBeta        # linear combination with a bias
  pr = 1/(1+exp(-z))         # pass through an inv-logit function
  y = rbinom(100,1,pr)      # bernoulli response variable
  
  binom.rf <- SMARTVIS::binomialRF(X,factor(y), fdr.threshold = .05,
                                   ntrees = 1000,percent_features = .6,
                                   fdr.method = 'BY')
  return(binom.rf)
  
}

nsim=1000
cv.bigMat2 = data.table(do.call(rbind,lapply(1:nsim, function(i) run_dpnd.sim())))
setkey(cv.bigMat2, Variable)
cv.bigMat2=cv.bigMat2[, list(cv.Freq =mean(Frequency), cv.Pvalue =median(Pvalue), cv.AdjPvalue=median(AdjPvalue), Avg.TimesSignificant=mean(Significant)), by=Variable]
print(cv.bigMat2)
write.csv(cv.bigMat2, file='../results/uncorr_sim.csv')

run_corr.sim <- function(){
  
  simData = simCorrData(pathwaySize = 10, .7, numSamples = 100, numPathways = 1)
  
  X = t(simData$X)
  y = simData$y
  
  binom.rf <- SMARTVIS::binomialRF(X,factor(y), fdr.threshold = .05,
                                   ntrees = 1000,percent_features = .6,
                                   fdr.method = 'BY')
  
}

nsim=1000
cv.bigMat3 = data.table(do.call(rbind,lapply(1:nsim, function(i) run_corr.sim())))
setkey(cv.bigMat3, Variable)
cv.bigMat3=cv.bigMat3[, list(cv.Freq =mean(Frequency), cv.Pvalue =median(Pvalue), cv.AdjPvalue=median(AdjPvalue), Avg.TimesSignificant=mean(Significant)), by=Variable]
print(cv.bigMat3)
write.csv(cv.bigMat3, file='../results/corr_sim.csv')


library(graphics)
single.modelViz <- function(cv.bigMat){
  mat <- t(matrix(sort(cv.bigMat$Avg.TimesSignificant, decreasing=T)))
  plot(y=seq(.1,1,by=.1), x=rep(1,10), pch='', xaxt='n', xlab='Predictors', yaxt='n', ylab='')
  image(mat, zlim=c(0,1), col=heat.colors(10, alpha=.5), add=T)
  axis(labels=paste('X',1:10), at=seq(.1,1,by=.1), side=2)
}

## rf1 
candidateModels <- list(m1=c('X1','X4','X6','X9'),
                        m2=c(paste('X',c(1,2,3,4,6,9),sep='')),
                        m3=c(paste("X", 3:10,sep='')),
                        m4= c(paste("X", 5:9,sep='')),
                        m5= c(paste("X", c(4,5,8,10),sep='')),
                        m6= c(paste("X",1:10,sep='')))

X = matrix(rnorm(1000), ncol=10)

### let half of the coefficients be 0, the other be 10
trueBeta= c(rep(10,5), rep(0,5))

### do logistic transform and generate the labels
z = 1 + X %*% trueBeta    
pr = 1/(1+exp(-z))        
y = rbinom(100,1,pr)


evaluateCandidateModels <- function(candidateModels, X, y){
  require(dplyr)
  
  X = data.frame(X)
  y = factor(y)
  
  f <- function(i, candidateModels, X,y){
    d= binomialRF(X[, candidateModels[[i]]], factor(y), ntrees = 5000)
    d$model = names(candidateModels)[i]
    return(d)
  }

  candidateList = lapply(1:length(candidateModels), function(i) f(i, candidateModels, X,y) )
  numModels <- length(candidateList)
  
  err.mat= data.frame(Variable= character(10*numModels),
             Significant=logical(10*numModels),
             weight=numeric(10*numModels),
             model=t(do.call(cbind, lapply(1:numModels, function(x) t(rep(paste('m',x,sep=''),10))))),
             stringsAsFactors = F)
  

  i=0
  for(model in candidateList){
    a = ((i*10)+1)
    err.mat[a:(a+nrow(model)-1),1:3] <- model[, c('Variable','Significant','weight')]
    i=i+1
  }
  
  err.mat$Significant = as.numeric(err.mat$Significant)
  err.mat$Norm.Weight = err.mat$weight / sum(unique(err.mat$weight))
  err.mat$Significant = as.logical(err.mat$Significant)
  err.mat= err.mat[err.mat$Variable!='',]
  
  return(err.mat)
}

build_weighted_viz <- function(err.mat){
  
  x = err.mat$model
  y = as.numeric(image.mat$Variable)
  z = image.mat$Significant
  
  M =data.frame(x,y,z)
  M = M[order(M$x),c(1,3)]
  m1 = M[M$x==1, c(1,3)]
  
  M = data.table(err.mat)
  M[, list(Variable, Significant), by=model]
  
  
  M <- t(as.matrix(M))
  #heatmap(image.mat2)
  heatmap(M, Rowv = NA, Colv = NA)
  
  heatmap(M,Rowv = NA, Colv = NA, scale='row', RowSideColors = y)
  
  heatmap(t(as.matrix(m1)))
  
  
  
  
}

require(ggplot2)
err.mat = evaluateCandidateModels(candidateModels = candidateModels, X,y)
err.mat_weighted <- err.mat[order(err.mat$Norm.Weight,decreasing=TRUE), ]
err.mat_weighted$model <- factor(err.mat_weighted$model, levels = unique(err.mat_weighted$model))
err.mat_weighted$Variable <- factor(err.mat_weighted$Variable, levels= names(X))
# ggplot(data=err.mat_weighted, aes(x=model, y=Variable,  fill=Significant, width=weight)) +
#   geom_raster(inherit.aes = F , aes(x=model,y=Variable,fill=Significant, width=Norm.Weight))+
#   theme_minimal()+ labs(title='Feature Selection by Model',
#                    x='Likeliest Candidate Model',
#                    y='Feature',
#                    fill="Significant")+
#   theme(
#     plot.title = element_text(color="red", size=20, face="bold.italic"),
#     axis.title.x = element_text(color="blue", size=14, face="bold"),
#     axis.title.y = element_text(color="#993333", size=14, face="bold")
#   )  


ggplot(data=err.mat_weighted, aes(x=model, y=Variable,  fill=Significant, width=Norm.Weight)) +
  geom_tile(inherit.aes = F , aes(x=model,y=Variable,fill=Significant, width=Norm.Weight), position = position_identity())+
  theme_minimal()+ labs(title='Feature Selection by Model',
                        x='Likeliest Candidate Model (OOB Error)',
                        y='Feature',
                        fill="Significant")+
  theme(
    plot.title = element_text(color="red", size=20, face="bold.italic"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="#993333", size=14, face="bold")
  )  

err.mat= data.table(err.mat)
final.err.mat = err.mat[, list(AvgSignificant= mean(Significant), Avg.Weight=mean(weight)), by =Variable]
final.err.mat = final.err.mat[order(final.err.mat$AvgSignificant, decreasing = T),]
