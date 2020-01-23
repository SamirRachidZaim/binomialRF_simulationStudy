################################
##
## test RF VIS v. LASSO
##
################################
rm(list=ls())
set.seed(10)

require(parallel)
source('helper_functions.R')
source('variableModelFunctions.R')

#install.packages('parallelRandomForest')
# devtools::install_github("bert9bert/ParallelForest")
# library(ParallelForest)
# devtools::install_bitbucket("mkuhn/parallelRandomForest", ref="parallelRandomForest")
# install.packages('r2VIM2/', repos=NULL, type='source')
# devtools::install_github("sivarajankumar/rf-ace")



ntrees = 500
percent_features = .4
percent_data=.5
pval_thresh = .05
dim.X = 10
numSamples = 100
fdr.method= 'BY'

rho = 0.33
successprob=  binomialRF::calculateBinomialP_Interaction(L = dim.X, percent_features = .4, K=2)
trials = ntrees

cbinom <- correlbinom::correlbinom(rho, successprob = successprob, trials = trials, model = 'kuk')



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
  require(devtools)
  require(AUCRF)
  require(EFS)
  require(binomialRF)
  require(correlbinom)
  source('~/Dropbox/Samir/binomialRF/R/k_binomialRF.R')


    X = matrix(rnorm(numSamples*dim.X),ncol=dim.X)
    z = X[,1] + X[,2] +  3* (X[,1] * X[,2]) + 1 
    y = factor(rbinom(numSamples, size=1, prob= 1/(1+exp(-z))))
    X = data.frame(X)
    
    plainX = X 
    
    n=ncol(X)
    combb=combn(n,2)
    names_interactions = paste(paste('X', combb[1,],sep=''),'_', paste('X', combb[2,],sep=''), sep='')
    
    res=apply(X, 1, function(x) { apply(combb, 2, function(y) prod(x[y])) })
    interacts= data.frame(t(res)); names(interacts) <- names_interactions
    
    X = cbind(X,interacts)

    numVars = ncol(X)

    ##### Generate simulated data
    rf.object <- randomForest(X, y, keep.forest = T, keep.inbag = T, ntree = ntrees)
    
    #### binomRF
    binom.rf.time = system.time(binom.rf <- .k_binomialRF(plainX,y, fdr.method = fdr.method, ntrees=ntrees, percent_features = .4,K=2, cbinom=cbinom))

    ### VSURF
    vsurf.time <- system.time(vsurf.object <- VSURF(X, y, parallel = T))

    ### Boruta
    boruta.time <- system.time(boruta.object <-  Boruta(X, y, verbose=1))
    
    ### Altmann PIMP
    ### use Kolgomorov-Smirnoff test instead of
    ### specifying dist'n

    pimp.time <- system.time(pimp.test <- data.frame(PimpTest(PIMP(X,y, rf.object, S=100, parallel = T, ncores=4), para=T)$p.ks.test, stringsAsFactors = F))
    
    vita_function <- function(X,y, ntrees){
      ##################################################################
      # cross-validated permutation variable importance
      mtry  = round( dim.X * .4)
      cv_vi = CVPVI(X,y,k = 5,mtry =mtry,ntree = ntrees,ncores = 4)

      # Novel Test approach
      cv_p = NTA(cv_vi$cv_varim)
      summary(cv_p,pless = 0.1)

    }
    
    vita.time <- system.time(vita.test <-vita_function(X,y,ntrees) )

    ### varSelRF
    varselrf.time <- system.time(varselrf.object <- varSelRF(xdata = X, Class = y))
    
    ### AUCRF 
    joint.df <- cbind(X,y)
    
    aucrf.time <- system.time(aucrf.object <- AUCRF(y~., data = joint.df))
    
    ### perm time
    
    perm_function<- function(X,y,ntrees){
      rf.object <- randomForest(X, y, keep.forest = T, keep.inbag = T, ntree = ntrees)
      imp_df <- importance(rf.object)
      
      permuted_y <- sample(y)
      
      rf.object2 <- randomForest(X, permuted_y, keep.forest = T, keep.inbag = T, ntree = ntrees)
      imp_df2 <- importance(rf.object2)
      
      selected_vars <- imp_df[which(imp_df- imp_df2 >0),1]
      
    }
    
    perm.time <- system.time(perm.object <- perm_function(X,y, ntrees))
    
    ### EFS 
    joint.df$y <- as.numeric(joint.df$y)-1
    classY= dim.X+choose(dim.X,2)+1
    efs.time <- system.time(EFS::ensemble_fs(data=joint.df, classnumber = classY, runs = 50))
    
    ### RFE
    
    rfe_function <- function(X,y,ntrees){
      rf.object <- randomForest(X, y, keep.forest = T, keep.inbag = T, ntree = ntrees)
      
      idx <- order(importance(rf.object), decreasing = T)
      ranked_imp <- data.frame(Variable = row.names(importance(rf.object))[idx],
                               Importance= importance(rf.object)[idx])
      
      last_obs = nrow(ranked_imp) *.8
      ranked_imp <- ranked_imp[1:last_obs, ]
      
      rf.object2 <- randomForest(X[,ranked_imp$Variable], y, keep.forest = T, keep.inbag = T, ntree = ntrees)
      idx <- order(importance(rf.object2), decreasing = T)
      ranked_imp2 <- data.frame(Variable = row.names(importance(rf.object2))[idx],
                               Importance= importance(rf.object2)[idx])
      
      last_obs = nrow(ranked_imp2) *.8
      ranked_imp2 <- ranked_imp2[1:last_obs, ]
      
      rf.object3 <- randomForest(X[,ranked_imp2$Variable], y, keep.forest = T, keep.inbag = T, ntree = ntrees)
      idx <- order(importance(rf.object3), decreasing = T)
      ranked_imp3 <- data.frame(Variable = row.names(importance(rf.object3))[idx],
                                Importance= importance(rf.object3)[idx])
      
      last_obs = nrow(ranked_imp3) *.8
      ranked_imp3 <- ranked_imp3[1:last_obs, ]
    }
    
    rfe.time <- system.time(rfe.object <- rfe_function(X, y, ntrees=ntrees))
    
    
    
    TimeTracking = data.frame(rbind(
                                    binom.rf.time,
                                    boruta.time,
                                    vsurf.time,
                                    varselrf.time,
                                    pimp.time,
                                    aucrf.time,
                                    efs.time,
                                    rfe.time,
                                    vita.time,
                                    perm.time
                                    
                                    ))
    TimeTracking$Model = c('BinomialRF','Boruta','VSURF','VarSelRF','PIMP','AUCRF','EFS','RFE', 'Vita','Perm')

    return(TimeTracking)
}


### Parallelize using
### Foreach package
dim.X =100

numSamples = 100
ntrees = 500
nsim   = 50
fdr.method ='BY'

require(doSNOW)
cl<-makeCluster(detectCores(), type = 'SOCK' ) # I have two cores
registerDoSNOW(cl)

modelList = foreach(j=1:nsim,
                    .packages=c('data.table','compiler','parallel')) %dopar%
  try(trainModels.time(ntrees = ntrees, fdr.method = fdr.method, dim.X = dim.X, numSamples = numSamples))

stopCluster(cl)
mod = do.call(rbind, modelList)

write.csv(mod,file=paste('../results/interactions_TimeProfile_', dim.X,'_features.csv', sep=''))
print(mod)

