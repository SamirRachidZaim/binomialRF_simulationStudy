##########################################################################################
##########################################################################################
#### load data and libraries
rm(list=ls())
require(binomialRF)
require(randomForest)
load('~/Dropbox/Samir/binomialRF_study/code/rhinoVirus/TCGA_all.RData')
require(TCGA2STAT)
gobp <- read.csv('~/Dropbox/Lab-Tools-Files/GeneOntology/2019-03-20/Human/gene2go_Human_BP_filter_15-500.txt', sep='\t',stringsAsFactors = F)

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


