################################
##
## test RF VIS v. LASSO
##
################################
set.seed(10)

require(parallel)
source('simulateData.R')
source('variableModelFunctions.R')
source('helper_functions.R')
source('../../binomialRF/R/k.binomialRF.R')
#source('~/Desktop/shiny/test/helper.R')
require(data.table)
require(Boruta)
require(VSURF)
require(vita)
require(randomForest)
require(varSelRF)
require(devtools)

## GenerateModel
trainModels <- function(ntrees=1000, percent_features=.4, percent_data=.3, fdr.method='BY', pval_thresh =.05,
                         numSamples = 100, dimX=100, BetaSignal=1, interaction.study=1){

    require(data.table)
    require(Boruta)
    require(VSURF)
    require(vita)
    require(randomForest)
    require(varSelRF)
    require(glmnet)
    require(devtools)
#    load_all("~/Dropbox/Samir/SMART_VIS/SMART_VIS/SMARTVIS")
    load_all('~/Dropbox/Samir/binomialRF/')
    ##### Generate simulated data

    ### Generate multivariate normal data in R10
    X = matrix(rnorm(numSamples*dimX), ncol=dimX)

    if(interaction.study==1){
          ### do logistic transform and generate the labels
          ## int 1
          z = 1 + X[,1]*BetaSignal + X[,3] *BetaSignal + (X[,1] * X[,3]) * BetaSignal
          K=2
    } else if (interaction.study==2){
          ## int 2
          z = 1 + X[,1]*BetaSignal + X[,2] *BetaSignal + X[,3]*BetaSignal + X[,4]*BetaSignal+ (X[,1] * X[,3]) * BetaSignal+ (X[,2] * X[,4]) * BetaSignal
          K=2
    } else if (interaction.study==3){
          ## int 3
          z = 1 + X[,1]*BetaSignal + X[,2] *BetaSignal + X[,3]*BetaSignal + (X[,1] * X[,3] * X[,3]) * BetaSignal
          K=3
    }
    ## int 3

    pr = 1/(1+exp(-z))
    y = factor(rbinom(numSamples,1,pr))
    X = data.frame(X)

# 
#     ## repeat for test data
#     test.X= matrix(rnorm(numSamples*dimX), ncol=dimX)
#     z.test= 1 + test.X[,1]*BetaSignal + test.X[,3]*BetaSignal + (test.X[,1] * test.X[,3]) * BetaSignal
#     pr.test= 1/(1+exp(-z.test))
#     test.y = factor(rbinom(numSamples,1,pr.test))
#     test.X = data.frame(test.X)



    numVars = ncol(X)

    rf.object <- randomForest(X, y, keep.forest = T, keep.inbag = T, ntree = ntrees)

    #### binomRF

    k.binom.rf <- k.binomialRF(X,y, K = K , fdr.method = fdr.method, ntrees=ntrees, percent_features = .5)

    return(k.binom.rf)

}

for(BetaSignal in c(1,3)){

  for(interaction.study in c(3)){

    dimX=10
    ### Parallelize using
    ### Foreach package
    fdr.method ='BY'

    ntrees = 500
    fdr    = 'BY'
    nsim   = 500

    require(doSNOW)
    cl<-makeCluster(detectCores(), type = 'SOCK' ) # I have two cores
    registerDoSNOW(cl)

    modelList = foreach(j=1:nsim,
                        .packages=c('data.table','compiler','parallel')) %dopar%
      try(trainModels(ntrees = ntrees, fdr.method = fdr, dimX = dimX, interaction.study = interaction.study, BetaSignal = BetaSignal))

    stopCluster(cl)
    final.mat = do.call(rbind,modelList )

    summary.mat <- final.mat[, lapply(.SD, mean), by=Interaction, .SDcols=2:4]
    summary.mat <- summary.mat[order(summary.mat$Adj.Pval),]

    summary.mat$Significant <- summary.mat$Adj.Pval < .05
    #
    write.csv(summary.mat,file=paste('../results/beta', BetaSignal,'/interaction', interaction.study,'_', dimX,'_P','.csv', sep=''))

  }

}






#
# fsr_simulations = mclapply(1:NVars, function(modelSize)
#   cbind(data.frame(do.call(rbind, mclapply(1:length(modelList), function(x)
#     calculateFSR_PerSample(betaNames, ntrees=ntrees, modelList = modelList[[x]], modelSize = modelSize) ))), modelSize))
#
# discovery_simulations = mclapply(1:NVars, function(modelSize)
#   cbind(data.frame(do.call(rbind, mclapply(1:length(modelList), function(x)
#     calculateDiscoveryRate_PerSample(betaNames, ntrees=ntrees, modelList = modelList[[x]], modelSize = modelSize) ))), modelSize))
#
#
# trainErr_simulations = data.frame(do.call(rbind, mclapply(1:length(modelList), function(x) modelList[[x]]$TrainError)))
# testErr_simulations = data.frame(do.call(rbind, mclapply(1:length(modelList), function(x) modelList[[x]]$TestError)))
#
# ## Permuted Matrix Noise
# permutedX = generatePermutedMatrix(X)
# permutedXtest = generatePermutedMatrix(Xtest)
#
# permutedModelList = mclapply(1:NSimulations, function(i) trainModels(permutedX,YTrainlist[[i]] , permutedXtest, ytest=YTestlist[[i]], ntrees=ntrees))
#
# fsr_permuted_simulations = mclapply(1:(2*NVars), function(modelSize)
#   cbind(data.frame(do.call(rbind, mclapply(1:length(permutedModelList), function(x)
#     calculateFSR_PerSample(betaNames, ntrees=ntrees, modelList = permutedModelList[[x]], modelSize=modelSize) ))), modelSize))
#
# discovery_permuted_simulations = mclapply(1:(2*NVars), function(modelSize)
#   cbind(data.frame(do.call(rbind, mclapply(1:length(permutedModelList), function(x)
#     calculateDiscoveryRate_PerSample(betaNames, ntrees=ntrees, modelList = permutedModelList[[x]], modelSize = modelSize) ))), modelSize))
#
#
# trainErr_permuted_simulations = data.frame(do.call(rbind, mclapply(1:length(permutedModelList), function(x) permutedModelList[[x]]$TrainError)))
# testErr_permuted_simulations = data.frame(do.call(rbind, mclapply(1:length(permutedModelList), function(x) permutedModelList[[x]]$TestError)))
#
#
#
#
#
#
#
#
#
#
#
#
#
# # plotFSR <- function(fsr_simulations){
# #
# #   J = 20
# #
# #   plot(x=1:10, y=1:10, xlim=c(1,J),ylim=c(0,1) ,type="n",ylab="False Selection Rate",xlab="Model Size", main='Controlling FSR Through Model Size')
# #
# #   for(l in 1:J)  {
# #     points(y=jitter(fsr_simulations[[l]]$SMART_VIS), x=rep(l, NSimulations), pch=0, col=1, cex=.4)
# #     points(y=mean(fsr_simulations[[l]]$SMART_VIS), x=l, pch=0, col=1, cex=1)
# #
# #     points(y=jitter(fsr_simulations[[l]]$GiniCoefficient),x= rep(l, NSimulations), pch=1, col=2, cex=.4)
# #     points(y=mean(fsr_simulations[[l]]$GiniCoefficient),x= l, pch=1, col=2, cex=1)
# #
# #     points(y=jitter(fsr_simulations[[l]]$Lasso), x=rep(l, NSimulations), pch=2, col=3, cex=.4)
# #     points(y=mean(fsr_simulations[[l]]$Lasso),x= l, pch=2, col=3, cex=1)
# #
# #
# #     points(y=jitter(fsr_simulations[[l]]$RFSoil), x=rep(l, NSimulations), pch=3, col=4, cex=.4)
# #     points(y=mean(fsr_simulations[[l]]$RFSoil),x= l, pch=3, col=4, cex=1)
# #
# #
# #     #segments( l, rank_summary$`Lower 95% CI`[j] ,l, rank_summary$`Upper 95% CI`[j] ,lwd=2)
# #   }
# #   legend('topleft', c('SMART_VIS','GINI','LASSO','RFSOIL'), col=c(1:4), pch=0:3)
# #   # legend('topright', c('FSR=5%','FSR=10%','FSR=20%'), col=c('green','orange','red'), lty=2, lwd=2)
# #
# #   abline(h=.05, lwd=2, col='green', lty=2)
# #   # abline(h=.1, lwd=2, col='orange', lty=2)
# #   abline(h=.2, lwd=2, col='red', lty=2)
# #
# #   plot(x=1:10, y=1:10, xlim=c(1,J),ylim=c(0,1) ,type="n",ylab="False Selection Rate",xlab="Model Size", main='Mean FSR ~ Model Size')
# #
# #   for(l in 1:J)  {
# #     points(y=mean(fsr_simulations[[l]]$SMART_VIS), x=l, pch=0, col=1, cex=1)
# #     points(y=mean(fsr_simulations[[l]]$GiniCoefficient),x= l, pch=1, col=2, cex=1)
# #     points(y=mean(fsr_simulations[[l]]$Lasso),x= l, pch=2, col=3, cex=1)
# #     points(y=mean(fsr_simulations[[l]]$RFSoil),x= l, pch=3, col=4, cex=1)
# #
# #     #segments( l, rank_summary$`Lower 95% CI`[j] ,l, rank_summary$`Upper 95% CI`[j] ,lwd=2)
# #   }
# #
# #
# #   #plot(test$SMART_VIS, type='l')
# #
# #   legend('topleft', c('SMART_VIS','GINI','LASSO','RFSOIL'), col=c(1:4), pch=0:3)
# #   # legend('topright', c('FSR=5%','FSR=10%','FSR=20%'), col=c('green','orange','red'), lty=2, lwd=2)
# #
# #   abline(h=.05, lwd=2, col='green', lty=2)
# #   # abline(h=.1, lwd=2, col='orange', lty=2)
# #   abline(h=.2, lwd=2, col='red', lty=2)
# # }
#
# plotConfidenceBandsFSR <- function(fsr_list, permuted=F, yTitle="False Selection Rate"){
#   test= do.call(rbind, fsr_list)
#   test2 =  data.table(melt(test, id.vars = 'modelSize',measure.vars = c('SMART_VIS','Lasso','GiniCoefficient','RFSoil')))
#   test3 = test2[, mean(value), by=list(variable, modelSize)]
#
#   if(permuted){
#     ggplot(test3, aes( modelSize, V1 ) ) +
#       geom_point(aes( colour = variable )) +
#       geom_smooth(aes( group = variable , colour=variable)) +
#       geom_line(aes( modelSize, y=rep(0.2,160) ),col='red', lty=2) +
#       geom_line(aes( modelSize, y=rep(0.05,160) ),col='green', lty=2) + ylim(0,1) +
#       ylab(yTitle) + xlab('Model Size') + theme_minimal()+
#       theme(axis.text=element_text(size=12),
#             axis.title=element_text(size=14,face="bold"))
#
#   } else{
#     ggplot(test3, aes( modelSize, V1 ) ) +
#       geom_point(aes( colour = variable )) +
#       geom_smooth(aes( group = variable , colour=variable)) +
#       geom_line(aes( modelSize, y=rep(0.2,80) ),col='red', lty=2) +
#       geom_line(aes( modelSize, y=rep(0.05,80) ),col='green', lty=2) +
#       ylab(yTitle) + xlab('Model Size') + theme_minimal()+ ylim(0,1) +
#       theme(axis.text=element_text(size=12),
#             axis.title=element_text(size=14,face="bold"))
#     }
#
#
#
# }
#
# #pdf('../plots/FSR_grid.pdf')
# require(ggpubr)
# p = plotConfidenceBandsFSR(fsr_simulations)
# g = plotConfidenceBandsFSR(fsr_permuted_simulations, permuted=T)
#
# ggarrange(p, g ,
#           labels = c("Normal Variables", "Normal + Phony Variables"),
#           ncol = 1, nrow = 2)
# #dev.off()
#
#
# require(ggpubr)
# p = plotConfidenceBandsFSR(discovery_simulations, yTitle = 'Discovery Rate')
# g = plotConfidenceBandsFSR(discovery_permuted_simulations, permuted=T, yTitle = 'Discovery Rate')
#
# ggarrange(p, g ,
#           labels = c("Normal Variables", "Normal + Phony Variables"),
#           ncol = 1, nrow = 2)
#
# #dev.off()
#
#
#
#
#
# #pdf('../plots/boxplot_grid.pdf')
#
# err1 =  data.table(melt(trainErr_simulations, id.vars = NULL ,measure.vars = c('SMART_VIS','LASSO','Gini','RFSoil')))
# err2 =  data.table(melt(testErr_simulations, id.vars = NULL ,measure.vars = c('SMART_VIS','LASSO','Gini','RFSoil')))
# err3 =  data.table(melt(trainErr_permuted_simulations, id.vars = NULL ,measure.vars = c('SMART_VIS','LASSO','Gini','RFSoil')))
# err4 =  data.table(melt(testErr_permuted_simulations, id.vars = NULL ,measure.vars = c('SMART_VIS','LASSO','Gini','RFSoil')))
#
# a1 = ggplot(err1, aes(x=variable, y=value)) + geom_boxplot() + xlab('Method') + ylim(c(0,.5)) +theme_minimal()+
#   ylab("Misclassification Rate") +  theme(axis.text=element_text(size=6), axis.title=element_text(size=10,face="bold"))
# a2 = ggplot(err2, aes(x=variable, y=value)) + geom_boxplot()+ xlab('Method') + ylim(c(0,.5)) +theme_minimal()+
#   ylab("Misclassification Rate") +  theme(axis.text=element_text(size=6),axis.title=element_text(size=10,face="bold"))
# a3 = ggplot(err3, aes(x=variable, y=value)) + geom_boxplot()+ xlab('Method') + ylim(c(0,.5)) +theme_minimal()+
#   ylab("Misclassification Rate") +  theme(axis.text=element_text(size=6),axis.title=element_text(size=10,face="bold"))
# a4 = ggplot(err4, aes(x=variable, y=value)) + geom_boxplot()+ xlab('Method') + ylim(c(0,.5)) +theme_minimal()+
#   ylab("Misclassification Rate") +  theme(axis.text=element_text(size=6), axis.title=element_text(size=10,face="bold"))
#
# ggarrange(a1, a2, a3, a4 ,
#           labels = c("Normal Vars: Training Error", "Normal Vars: Test Error", "Phony Vars: Training Error", "Phony Vars: Test Error"),
#           ncol = 2, nrow = 2, font.label = list(size=10))
# #dev.off()
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# # pdf('../plots/test.FSR.pdf')
# #
# # dev.new()
# # par( mfrow = c( 1, 2 ), oma = c( 8, 1, 4, 1 ) )
# # tname = paste('Y = ', 'I(',paste('X', c(1:NTrueBetas), sep='', collapse=' + '), ') > 3')
# #
# # boxplot(fsr_simulations, xlab='Method', ylab='FSR', main='FSR by Technique', ylim=c(0,.5))
# # boxplot(fsr_permuted_simulations, xlab='Method', ylab='FSR', main='FSR by Technique in Permuted Matrix', ylim=c(0,.5))
# #
# # title(main='False Selection Rate by Variable Selection Method',
# #       sub = paste(strwrap(tname), collapse='\n'),
# #       outer = TRUE , cex.main=2, cex.sub=1)
# # par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
# # plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
# #
# # legend( 'bottomright', c(paste('NTrain=', NTrain)
# #                          ,paste("Ntest=", NTest)
# #                          ,paste("|Total B|=", NVars)
# #                          ,paste("|True B|=",NTrueBetas)
# #                          ,paste("NSimulation=",NSimulations)),
# #         cex=1.3)
# #
# # dev.off()
# #
# # pdf('../plots/test.error.pdf')
# # dev.new()
# # par( mfrow = c( 2, 2 ), oma = c( 8, 1, 4, 1 ) )
# #
# # boxplot(trainErr_simulations, xlab='Method', ylab='Error', main='Training Error by Technique', ylim=c(0,.5))
# # boxplot(trainErr_permuted_simulations, xlab='Method', ylab='Error', main='Training Error by Technique in Permuted Matrix', ylim=c(0,.5))
# # boxplot(testErr_simulations, xlab='Method', ylab='Error', main='Test Error by Technique', ylim=c(0,.3))
# # boxplot(testErr_permuted_simulations, xlab='Method', ylab='Error', main='Test Error by Technique in Permuted Matrix', ylim=c(0,.3))
# #
# # title(main='Training & Test Error by Variable Selection Method',
# #       sub = paste(strwrap(tname)),
# #       outer = TRUE , cex.main=2, cex.sub=1)
# #
# # par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
# # plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
# # # legend("bottom", c("IM", "IBD", "1R", "2R"), xpd = TRUE, horiz = TRUE, inset = c(0,
# # #                                                                                  0), bty = "n", pch = c(4, 2, 15, 19), col = 1:4, cex = 2)
# # legend( 'bottomright', c(paste('NTrain=', NTrain)
# #                          ,paste("Ntest=", NTest)
# #                          ,paste("|Total B|=", NVars)
# #                          ,paste("|True B|=",NTrueBetas)
# #                          ,paste("NSimulation=",NSimulations)),
# #         cex=1.3)
# #
# #
# # dev.off()














