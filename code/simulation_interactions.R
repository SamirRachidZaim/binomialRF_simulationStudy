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
    
    z = X[,1] + X[,2] + X[,3] + X[,4] +3 * (X[,1] * X[,2]) +  X[,3] + X[,4] + 3 * (X[,1] * X[,3]) + 3 * (X[,1] * X[,4]) + 
      3 * (X[,2] * X[,3]) +  3 * (X[,2] * X[,4]) +  3 * (X[,3] * X[,4])
    
    y = factor(rbinom(numSamples, size=1, prob= 1/(1+exp(-z))))
    X = data.frame(X)
    
    sampsize = .25 * nrow(X)
    
    ## trueVariable vector
    trueInterVars = c('X1_X2','X1_X3','X1_X4','X2_X3','X2_X4','X3_X4')
    
    
    
    plainX = X 
    
    n=ncol(X)
    combb=combn(n,2)
    names_interactions = paste(paste('X', combb[1,],sep=''),'_', paste('X', combb[2,],sep=''), sep='')
    
    res=apply(X, 1, function(x) { apply(combb, 2, function(y) prod(x[y])) })
    interacts= data.frame(t(res)); names(interacts) <- names_interactions
    
    X = cbind(X,interacts)
    
    
    ## repeat for test data
    test.X= matrix(rnorm(numSamples*dim.X), ncol=dim.X)
    z.test = test.X[,1] + test.X[,2] + test.X[,3] + test.X[,4] +3 * (test.X[,1] * test.X[,2]) +  test.X[,3] + test.X[,4] + 3 * (test.X[,1] * test.X[,3]) + 3 * (test.X[,1] * test.X[,4]) + 
      3 * (test.X[,2] * test.X[,3]) +  3 * (test.X[,2] * test.X[,4]) +  3 * (test.X[,3] * test.X[,4])
    pr.test= 1/(1+exp(-z.test))
    test.y = factor(rbinom(numSamples,1,pr.test))
    test.X = data.frame(test.X)
    
    plainX.text = test.X
    
    n=ncol(test.X)
    combb=combn(n,2)
    names_interactions = paste(paste('X', combb[1,],sep=''),'_', paste('X', combb[2,],sep=''), sep='')
    
    res=apply(test.X, 1, function(x) { apply(combb, 2, function(y) prod(x[y])) })
    interacts= data.frame(t(res)); names(interacts) <- names_interactions
    
    test.X = cbind(test.X,interacts)
    
    
    #####
    ########################################################################
    ########################################################################
    


    numVars = ncol(X)
    rf.object <- randomForest(X, y, keep.forest = T, keep.inbag = T, ntree = ntrees, sampsize = sampsize)
  
    
    ########################################################################
    ########################################################################
    
    #### binomRF
    binom.rf <- binomialRF::k_binomialRF(plainX,y, fdr.method = 'bonferroni', ntrees=ntrees, percent_features = .6, cbinom_dist  = cbinom)
    binom.rf$Significant <- binom.rf$adjSignificance < .001
    binom.rf.vars <- binom.rf$Interaction[which(binom.rf$Significant ==T)]
    
    
    # if(length(binom.rf.vars) < 2){
    #   binom.rf.object= rf.object
    # } else {
    #   binom.rf.object <- randomForest(X[,binom.rf.vars], y)
    # }
    
    
    binom.rf.vars = stringr::str_replace(binom.rf.vars, ' \\| ', '_')

    binom.rf.size <- length(binom.rf.vars)
    binom.rf.precision <- sum(binom.rf.vars %in% trueInterVars)/length(binom.rf.vars)
    binom.rf.recall <- sum(binom.rf.vars %in% trueInterVars)/ length(trueInterVars)
    
    # binom.rf.test.err  <- mean(predict(binom.rf.object, test.X)!= test.y)
    
    ########################################################################
    ########################################################################
    


    ########################################################################
    ########################################################################
    ### VSURF
    vsurf.object = VSURF(X, y, parallel = T, ntree = ntrees/2, sampsize=sampsize)
    vsurf.vars = vsurf.object$varselect.pred
    vsurf.vars <- paste('X', vsurf.vars, sep='')

    if(length(vsurf.vars) < 2){
	vsurf.rf.object= rf.object
    } else {
	vsurf.rf.object <- randomForest(X[,vsurf.vars], y)
    }
    
    vsurf.size <- length(vsurf.vars)
    vsurf.precision <- sum(vsurf.vars %in% trueInterVars) / length(vsurf.vars)
    vsurf.recall <- sum(vsurf.vars %in% trueInterVars) / length(trueInterVars)
    # vsurf.test.err <-  mean(predict(vsurf.rf.object, test.X)!=test.y)
    
    ########################################################################
    ########################################################################
    
    
    
    ########################################################################
    ########################################################################
    ### Boruta
    boruta.object = Boruta(X, y, verbose=1)
    boruta.vars <- names(X)[which(boruta.object$finalDecision=='Confirmed')]
    #boruta.rf.object <- randomForest(X[, boruta.vars], y)
    boruta.vars <- boruta.vars[grep('_', boruta.vars)]
    
    
    boruta.size <- length(boruta.vars)
    boruta.precision <- sum(boruta.vars %in% trueInterVars) / length(boruta.vars)
    boruta.recall <- sum(boruta.vars%in% trueInterVars) / length(trueInterVars)
    # boruta.test.err <- mean(predict(boruta.rf.object, test.X)!=test.y)
    
    ########################################################################
    ########################################################################
    
    
    
    ########################################################################
    ########################################################################
    
    ### Atlmann PIMP
    ### use Kolgomorov-Smirnoff test instead of
    ### specifying dist'n

    pimp = PIMP(X,y, rf.object, S=100, parallel = T, ncores=4)
    pimp.test = data.frame(PimpTest(pimp, para=T)$p.ks.test, stringsAsFactors = F)
    pimp.test$adj.ks.test <- p.adjust(pimp.test$ks.test, method=fdr.method)

    pimp.test.vars <- row.names(pimp.test)[which(pimp.test$adj.ks.test < pval_thresh)]

    # if(length(pimp.test.vars) > 1){
    # pimp.object <- randomForest(X[, pimp.test.vars], y)
    # } else {
    # pimp.object <- rf.object
    # }

    pimp.size <- length(pimp.test.vars)
    pimp.precision <- sum(pimp.test.vars %in% trueInterVars) / length(pimp.test.vars)
    pimp.recall <- sum(pimp.test.vars %in% trueInterVars) / length(trueInterVars)
    # pimp.test.err <- mean(predict(pimp.object, test.X)!=test.y)
    
    ########################################################################
    ########################################################################
    
    
    
    ########################################################################
    ########################################################################
    
    ### aucRF
    ### technique that chooses variable based on maximizing
    ### AUC

    newData <- as.data.frame(cbind(X,y))
    colnames(newData) <- c(paste('X', 1:dim.X, sep=''), 'y')
    newData$y <- as.factor(newData$y)

    aucrf.object <- AUCRF(y~., data=newData)

    aucrf.vars <- aucrf.object$Xopt

    # aucrf.object <- randomForest(X[, aucrf.vars], y)
    
    
    aucrf.size = length(aucrf.vars)
    aucrf.precision <- sum(aucrf.vars %in% trueInterVars) / aucrf.size
    aucrf.recall <- sum(aucrf.vars %in% trueInterVars) / length(trueInterVars)
    # aucrf.test.err <- mean(predict(aucrf.object, test.X)!=test.y)
    
    ########################################################################
    ########################################################################
    


    ########################################################################
    ########################################################################
    ### varSelRF
    varselrf.object = varSelRF(xdata = X, Class = y)
    varselrf.var = varselrf.object$selected.vars

    # varselrf.object <- randomForest(X[, varselrf.var], y)
    
    
    varselrf.size <- length(varselrf.var)
    varselrf.precision <- sum(varselrf.var %in% trueInterVars) / varselrf.size
    varselrf.recall <- sum(varselrf.var %in% trueInterVars) / length(trueInterVars)
    # varselrf.test.err <- mean(predict(varselrf.object, test.X)!=test.y)

    ########################################################################
    ########################################################################
    
    
    ########################################################################
    ########################################################################
    ### Vita 
    
    vita_function <- function(X,y, ntrees){
      # cross-validated permutation variable importance
      mtry  = round( dim.X * .4)
      cv_vi = CVPVI(X,y,k = 5,mtry =mtry,ntree = ntrees,ncores = 4)
      
      # Novel Test approach
      cv_p = NTA(cv_vi$cv_varim)
      summary(cv_p,pless = 0.1)
      
    }
    
    vita.object <-vita_function(X,y,ntrees) 
    idx.vita <- which(vita.object$cmat[,2] < pval_thresh)
    vita.vars <- row.names(vita.object$cmat)[idx.vita]
    # vita.object <- randomForest(X[, vita.vars], y)
    
    
    vita.size <- length(vita.vars)
    vita.precision <- sum(vita.vars %in% trueInterVars) / length(vita.vars)
    vita.recall <- sum(vita.vars %in% trueInterVars) / length(trueInterVars)
    # vita.test.err <- mean(predict(vita.object, test.X)!=test.y)
    
    ########################################################################
    ########################################################################
    
    ########################################################################
    ########################################################################
    ### perm 
    
    perm_function<- function(X,y,ntrees){
      rf.object <- randomForest(X, y, keep.forest = T, keep.inbag = T, ntree = ntrees)
      imp_df <- importance(rf.object)
      
      permuted_y <- sample(y)
      
      rf.object2 <- randomForest(X, permuted_y, keep.forest = T, keep.inbag = T, ntree = ntrees)
      imp_df2 <- importance(rf.object2)
      
      selected_vars <- imp_df[which(imp_df- imp_df2 >0),1]
      return(selected_vars)
    }
    
    perm.object <- perm_function(X,y, ntrees)
    perm.vars <- names(perm.object)
    
    # perm.object <- randomForest(X[, perm.vars], y)
    
    
    perm.size <- length(perm.vars)
    perm.precision <- sum(perm.vars %in% trueInterVars) / perm.size
    perm.recall <- sum(perm.vars %in% trueInterVars) / length(trueInterVars)
    # perm.test.err <- mean(predict(perm.object, test.X)!=test.y)
    
    ########################################################################
    ########################################################################

    ########################################################################
    ########################################################################
    ### EFS 
    
    newData$y <- as.numeric(newData$y)-1
    select.idx <- c(T,T,T,T,F,T,F,F)
    
    efs.object <- EFS::ensemble_fs(data=newData, classnumber = dim.X+1, runs = 20,
                                  selection = select.idx)
    

    efs.vars <- names(which(colSums(efs.object) > 0.5))
    # efs.object <- randomForest(X[, efs.vars], y)
    
    
    efs.size <- length(efs.vars)
    efs.precision <- sum(efs.vars %in% trueInterVars) / length(efs.vars)
    efs.recall <- sum(efs.vars %in% trueInterVars) / length(trueInterVars)
    # efs.test.err <- mean(predict(efs.object, test.X)!=test.y)
    
        
    ########################################################################
    ########################################################################
    
    
    
    ########################################################################
    ########################################################################
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
      return(ranked_imp3)
    }
    
    rfe.object <- rfe_function(X, y, ntrees=ntrees)
    rfe.vars <- as.character(rfe.object$Variable[which(rfe.object$Importance > mean(rfe.object$Importance))])
    # rfe.object <-  randomForest(X[, rfe.vars], y)
    
    
    
    rfe.size <- length(rfe.vars)
    rfe.precision <- sum(rfe.vars %in% trueInterVars) / rfe.size
    rfe.recall <- sum(rfe.vars %in% trueInterVars) / length(trueInterVars)
    # rfe.test.err <- mean(predict(rfe.object, test.X)!=test.y)
    
    ########################################################################
    ########################################################################
    
    
    
    

    return(
                Eval.Mat = data.frame(
                   TestError = c(binom.rf.test.err,
                                 vsurf.test.err,
                                 boruta.test.err,
                                 varselrf.test.err,
                                 pimp.test.err,
                                 aucrf.test.err,
                                 perm.test.err,
                                 efs.test.err,
                                 rfe.test.err,
                                 vita.test.err
                                 ),

                   ModelSize = c(length(binom.rf.vars),
                                 length(vsurf.vars),
                                 length(boruta.vars),
                                 length(varselrf.var),
                                 length(pimp.test.vars),
                                 length(aucrf.vars),
                                 length(perm.vars),
                                 length(efs.vars),
                                 length(rfe.vars),
                                 length(vita.vars)),

                    Model = c('binomialRF','VSURF','Boruta','VarSelRF','PIMP', 'AUCRF', 'Perm','EFS','RFE','Vita')
                )
               )

}

dim.X= 30
sampsize <- rho <- 0.25
ntrees <- trials <-  500
fdr.method='BY'

pval_thresh =0.05
fdr.method='BY'
numSamples=1000

successProb = calculateBinomialP_Interaction(L=dim.X, percent_features = .5, K=2)

cbinom <- correlbinom::correlbinom(rho, successprob = successProb, trials = trials, model = 'kuk')
write.csv(cbinom, '../cbinom_distributions/interaction_cbinom.csv')
plot(cbinom)

### Parallelize using
### Foreach package
args <- commandArgs(TRUE)
nsim   = 5

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
write.csv(final.mat,file=paste('InteractionSimulation_', dimX,'_P','.csv', sep=''))


