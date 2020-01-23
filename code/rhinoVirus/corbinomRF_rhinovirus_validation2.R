# binomialRF virogram 
require(data.table)
require(randomForest)
require(foreign)
# training.rhv <-  fread('~/Dropbox/Virogram/Java/ClassificationVirogram/HRV/summary.binary.signed.txt', data.table = F); row.names(rhv) <- rhv$V1; rhv= rhv[,-1]
# 

require(binomialRF)



as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}


# asthma.validation. <-  fread('~/Dropbox/Virogram/Java/ClassificationVirogram/Ex Vivo/summary.binary.signed.txt', data.table = F); row.names(rhv) <- rhv$V1; rhv= rhv[,-1]
# flu <-  fread('~/Dropbox/Virogram/Java/ClassificationVirogram/FLU/summary.binary.signed.txt', data.table = F); row.names(flu) <- flu$V1; flu= flu[,-1]
rhv.pheno <- foreign::read.arff('~/Dropbox/Virogram/Java/ClassificationVirogram/HRV/summary.binary.signed.txt.arff')
exvivo.pheno <- foreign::read.arff('~/Dropbox/Virogram/Java/ClassificationVirogram/Ex Vivo/summary.binary.signed.txt.arff')
flu.pheno <-  foreign::read.arff('~/Dropbox/Virogram/Java/ClassificationVirogram/FLU/summary.binary.signed.txt.arff')
true.pathways <- fread('~/Dropbox/Samir/binomialRF_study/old submission folder/Table2_pathways.csv')


go.bp.description <- fread("~/Dropbox/Virogram/Pathways/go.bp.description.txt")
go.bp.description$Pathway <- go.bp.description$path_id


X = rhv.pheno[,! names(rhv.pheno) %in% c('Class', 'Name')]
y = rhv.pheno[, names(rhv.pheno) %in% 'Class']
Xtest <- exvivo.pheno[, !names(exvivo.pheno) %in% c("Class",'Name')]
ytest <-  exvivo.pheno[, names(exvivo.pheno) %in% c("Class")]


full.intersect = intersect(colnames(X), colnames(Xtest))
X = X[, full.intersect]
Xtest = Xtest[, full.intersect]

newX <- data.frame(do.call(cbind,lapply(1:ncol(X), function(zz) as.numeric(as.character(X[,zz]))))); names(newX) <- full.intersect
newXtest <- data.frame(do.call(cbind,lapply(1:ncol(Xtest), function(zz) as.numeric(as.character(Xtest[,zz]))))); names(newXtest) <- full.intersect

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### 
#### Main Effects STUDY
####  
#### look at Pathway-based classifiers
#### 
#### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
set.seed(1)
rho=0.2
ntrees = 1000

cbinom = correlbinom::correlbinom(rho, successprob = 1/ncol(newX), trials=1000, model='kuk')
bin.rf <- binomialRF::binomialRF(newX,y, ntrees = 1000, fdr.threshold = 0.05, user_cbinom_dist = cbinom)

trueFeatures =bin.rf$variable[bin.rf$adjSignificance < 0.1]
sum(trueFeatures %in% true.pathways$Pathway)

predicted.features =bin.rf[bin.rf$adjSignificance < 0.05,c(1:4)]
predicted.features$Pathway = predicted.features$variable
predicted.features <-  plyr::join(predicted.features, go.bp.description[, c('Pathway','description')])

write.csv(predicted.features, file='~/Dropbox/Samir/binomialRF_study/results/HRV/sample_binomialRF.csv',row.names=F)

predicted.features = data.frame(Pathway=trueFeatures)
predicted.features <-  plyr::join(predicted.features, go.bp.description[, c('Pathway','description')])

final.RF2 <- randomForest(newXtest[, names(newXtest) %in% trueFeatures], ytest )
final.RF2
OOB.error <- last(final.RF2$err.rate[,1])
OOB.error

naiveRF <- randomForest(newXtest, ytest)
OOB.error.naive <- last(naiveRF$err.rate[,1])
OOB.error.naive

compareModels <- function(){
  final.RF2 <- randomForest(newXtest[, names(newXtest) %in% trueFeatures], ytest )
  final.RF2
  OOB.error <- last(final.RF2$err.rate[,1])
  OOB.error
  
  
  naiveRF <- randomForest(newXtest, ytest)
  OOB.error.naive <- last(naiveRF$err.rate[,1])
  OOB.error.naive
  return(c(OOB.error,OOB.error.naive))
}

test = data.frame(do.call(rbind, lapply(1:100, function(x) compareModels())))
colnames(test) <- c('binomialRF', 'naive RF')
summary(test)

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### 
#### INTERACTION STUDY
####  
#### look at Pathway-Pathway Interaction 
#### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
require(reshape)

prob = calculateBinomialP_Interaction(ncol(newX[,trueFeatures]), .5, 2)
cbinom = correlbinom::correlbinom(rho, successprob =prob, trials=1000, model='kuk')

k.bin.rf <- k_binomialRF(X=newX[,trueFeatures],y=y, ntrees =2000, percent_features = .4, fdr.threshold = 0.05, K=2,cbinom_dist=cbinom)
k.bin.rf <- k.bin.rf$Interaction[k.bin.rf$adjSignificance < .1]
k.bin.rf <- data.frame(k.bin.rf, stringsAsFactors = F)

k.bin.rf <- transform(k.bin.rf, Interaction = colsplit(k.bin.rf, pattern = " | ", names = c('pathway1', 'pathway2')))

interaction.matrix <- data.frame(Interaction = k.bin.rf$k.bin.rf,
                                 Pathway1 = as.character(k.bin.rf$Interaction.pathway1),
                                 Pathway2 = as.character(k.bin.rf$Interaction.pathway2),
                                 stringsAsFactors = F)

interaction.matrix$Pathway2 <- sapply(interaction.matrix$Pathway2, function(x) gsub(x, pattern = '\\| ', replacement =''))

go.bp.description <- data.frame(go.bp.description)

interaction.matrix <- plyr::join(interaction.matrix, go.bp.description ,by=c('Pathway1'=='path_id'))

final.interaction.matrix <- cbind(k.bin.rf, interaction.matrix[, c('description.x','description.y')])
write.csv(final.interaction.matrix, '~/Dropbox/Samir/binomialRF_study/results/HRV case study/pathway_pathway_interactions.csv', row.names = F)














































