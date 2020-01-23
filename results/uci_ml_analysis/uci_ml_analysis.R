##################################################
##################################################
##########
##########
rm(list=ls())

##################################################


##################################################
##################################################

## load libraries
library(binomialRF)
require(data.table)
require(parallel)
source('fast_correlbinom.R')
set.seed(324)
setwd("~/Dropbox/Samir/binomialRF_study/code/uci_ml_analysis")

fdr.method = 'BY'
pval_thresh = 0.001


readInData <- function(directory){
  
  Xvalid <- fread(paste('data/',directory,'/',directory,'_valid.data',sep=''), data.table = F)
  Xtrain <- fread(paste('data/',directory,'/',directory,'_train.data',sep=''), data.table = F)
  
  ytrain <- fread(paste('data/',directory,'/',directory,'_train.labels',sep=''), data.table = F)
  yvalid <- fread(paste('data/',directory,'/',directory,'_valid.labels',sep=''), data.table = F)
  
  #params <- fread(paste('data/',directory,'/',directory,'.param',sep=''))
  dataType = paste(unique(sapply(1:ncol(Xvalid), function(i) class(Xvalid[,i]))),collapse = ' | ')
  
  params <- c(dim(Xvalid),dataType, directory, 'Classification')
  
  return(list(Xvalid=Xvalid,
              Xtrain=Xtrain,
              yvalid=yvalid,
              ytrain=ytrain,
              params=params))
}

directories <- dir('data')
directories <- directories[!directories %in% c('otherFormat')]

data_list <- list()
for(direc in directories){
  data_list[[direc]] <- readInData(direc)
}

DataDict <- as.data.frame(do.call(rbind, lapply(directories, function(x) data_list[[x]]$params)))
colnames(DataDict) <- c('Instances','Attributes','AttributeType','Dataset', 'Task')
# write.csv(DataDict, '../results/dataDict.csv', row.names=F)


analyzeData <- function(directory){
  
  dt <- data_list[[directory]]
  
  rho = 0.25
  X <- dt$Xtrain
  y <- as.factor(dt$ytrain$V1)
  X_test = dt$Xvalid
  y_test = as.factor(dt$yvalid$V1)
  
  dim.X = ncol(X)
  
  ntrees=1000
  
  cbinom_dist <- fast_correlbinom(rho=rho, successprob = 1/ncol(X),trials = ntrees,model = 'kuk')
  
  run_replicates <- function(X,y, X_test,y_test, cbinom_dist = cbinom_dist){
    
    idx <- sample.int(nrow(X), size=nrow(X) * .25)
    
    rf.object <- randomForest(X[idx,], y[idx], ntrees=ntrees)
      
    ## CorBinomialRF
    binomRF = binomialRF(X[idx,],y[idx] ,
                                  fdr.threshold=.05,
                                  fdr.method='BH',
                                  ntrees=ntrees,
                                  percent_features=.2,
                                  user_cbinom_dist=cbinom_dist,
                                  sampsize= nrow(X[idx,])*rho)
    
    cor.binomRF_vars <- as.character(binomRF$variable[ binomRF$adjSignificance < .05])
    cor.binomRF_vars <- stringr::str_replace(string = cor.binomRF_vars,pattern = 'X','V')
    
    cor.rf <- randomForest::randomForest(X[idx,cor.binomRF_vars], y[idx])
    cor.rf_err <- mean(!predict(cor.rf, X_test) ==y_test)
    cor.rf_modelSize <- length(cor.binomRF_vars)
    
    
    ########################################################################
    ########################################################################
    ### VSURF
    
    sampsize = rho * nrow(X[idx,])
    vsurf.object = VSURF(X[idx,],y[idx], parallel = T, ntree = ntrees/2, sampsize=sampsize)
    vsurf.vars = vsurf.object$varselect.pred
    vsurf.vars <- paste('V', vsurf.vars, sep='')
    
    if(length(vsurf.vars) < 2){
      vsurf.rf.object= rf.object
    } else {
      vsurf.rf.object <- randomForest(X[,vsurf.vars], y)
    }
    
    vsurf.size <- length(vsurf.vars)
    vsurf.test.err <-  mean(predict(vsurf.rf.object, X_test)!=y_test)
    
    ########################################################################
    ########################################################################
    
    
    
    ########################################################################
    ########################################################################
    ### Boruta
    boruta.object = Boruta(X[idx,],y[idx], verbose=1)
    boruta.vars <- names(X)[which(boruta.object$finalDecision=='Confirmed')]
    boruta.rf.object <- randomForest(X[, boruta.vars], y)
    boruta.size <- length(boruta.vars)
    boruta.test.err <- mean(predict(boruta.rf.object, X_test)!=y_test)
    
    ########################################################################
    ########################################################################
    
    
    
    ########################################################################
    ########################################################################
    
    ### Atlmann PIMP
    ### use Kolgomorov-Smirnoff test instead of
    ### specifying dist'n
    
    pimp = PIMP(X[idx,],y[idx], rf.object, S=50, parallel = T, ncores=4)
    pimp.test = data.frame(PimpTest(pimp, para=T)$p.ks.test, stringsAsFactors = F)
    pimp.test$adj.ks.test <- p.adjust(pimp.test$ks.test, method=fdr.method)
    
    pimp.test.vars <- row.names(pimp.test)[which(pimp.test$adj.ks.test < pval_thresh)]

    if(length(pimp.test.vars) > 1){
    pimp.object <- randomForest(X[, pimp.test.vars], y)
    } else {
    pimp.object <- rf.object
    }
    
    pimp.size <- length(pimp.test.vars)
    pimp.test.err <- mean(predict(pimp.object, X_test)!=y_test)
    
    ########################################################################
    ########################################################################
    
    
    
    ########################################################################
    ########################################################################
    
    ### aucRF
    ### technique that chooses variable based on maximizing
    ### AUC
    
    newData <- as.data.frame(cbind(X[idx,],y[idx]))
    colnames(newData) <- c(paste('V', 1:ncol(X), sep=''), 'y')
    newData$y <- as.factor(newData$y)
    
    aucrf.object <- AUCRF(y~., data=newData)
    
    aucrf.vars <- aucrf.object$Xopt
    
    aucrf.object <- randomForest(X[, aucrf.vars], y)
    
    
    aucrf.size = length(aucrf.vars)
    aucrf.test.err <- mean(predict(aucrf.object, X_test)!=y_test)
    
    ########################################################################
    ########################################################################
    
    
    
    ########################################################################
    ########################################################################
    ### varSelRF
    varselrf.object = varSelRF(xdata = X[idx,], Class = y[idx])
    varselrf.var = varselrf.object$selected.vars
    
    varselrf.object <- randomForest(X[, varselrf.var], y)
    varselrf.size <- length(varselrf.var)
    varselrf.test.err <- mean(predict(varselrf.object, X_test)!=y_test)
    
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
    
    vita.object <-vita_function(X[idx,],y[idx],ntrees) 
    idx.vita <- which(vita.object$cmat[,2] < pval_thresh)
    vita.vars <- row.names(vita.object$cmat)[idx.vita]
    vita.object <- randomForest(X[, vita.vars], y)
    vita.size <- length(vita.vars)
    vita.test.err <- mean(predict(vita.object, X_test)!=y_test)
    
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
    
    perm.object <- perm_function(X[idx,],y[idx], ntrees)
    perm.vars <- names(perm.object)
    
    perm.object <- randomForest(X[, perm.vars], y)
    perm.size <- length(perm.vars)
    perm.test.err <- mean(predict(perm.object, X_test)!=y_test)
    
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
    efs.object <- randomForest(X[, efs.vars], y)
    efs.size <- length(efs.vars)
    efs.test.err <- mean(predict(efs.object, X_test)!=y_test)
    
    
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
    
    rfe.object <- rfe_function(X[idx,],y[idx], ntrees=ntrees)
    rfe.vars <- as.character(rfe.object$Variable[which(rfe.object$Importance > mean(rfe.object$Importance))])
    rfe.object <-  randomForest(X[, rfe.vars], y)
    rfe.size <- length(rfe.vars)
    rfe.test.err <- mean(predict(rfe.object, X_test)!=y_test)
    
    ########################################################################
    ########################################################################
    Eval.Mat = data.frame(
      TestError = c(cor.rf_err,
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
      
      ModelSize = c(length(cor.binomRF_vars),
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
    return(Eval.Mat)
  }
  
  results <- mclapply(1:50, function(zz) run_replicates(X,y, X_test,y_test, cbinom_dist = cbinom_dist))
  return(results) 
}


analyzeData_cmp <- compiler::cmpfun(analyzeData)
resultsMadelon <- do.call(rbind, analyzeData(directory = 'madelon'))
# resultsArcene <- do.call(rbind, analyzeData_cmp(directory = 'arcene'))


resultsMadelon = data.frame(resultsMadelon)
backup = resultsMadelon

str(resultsMadelon)
clean_resultsMadelon = resultsMadelon[-grep('Error in socket', resultsMadelon$TestError),]

clean_resultsMadelon = data.table(clean_resultsMadelon)
clean_resultsMadelon$TestError <- as.numeric(clean_resultsMadelon$TestError)
clean_resultsMadelon$ModelSize <- as.numeric(clean_resultsMadelon$ModelSize)
clean_resultsMadelon_dt = clean_resultsMadelon[, list(median(TestError), sd(TestError), median(ModelSize), sd(ModelSize)), by=Model ]
clean_resultsMadelon_dt[, 2:5] <- round(clean_resultsMadelon_dt[, 2:5],2)
clean_resultsMadelon_dt$Error <- paste(clean_resultsMadelon_dt$V1, ' (', clean_resultsMadelon_dt$V2, ')')
clean_resultsMadelon_dt$Size <- paste(clean_resultsMadelon_dt$V3, ' (', clean_resultsMadelon_dt$V4, ')')
summarytable <- clean_resultsMadelon_dt[, c('Model','Error','Size')]
write.csv(summarytable, file='summarytable.csv', row.names = F)

dt= melt(clean_resultsMadelon, id.vars = 'Model')
ggplot(dt[variable=='TestError']) + geom_boxplot(aes(x = reorder(Model, value, median), y = value)) + ylim(c(0,1)) +title('')











resultsGisette <- do.call(rbind, analyzeData(directory = 'gisette'))
# resultsArcene











