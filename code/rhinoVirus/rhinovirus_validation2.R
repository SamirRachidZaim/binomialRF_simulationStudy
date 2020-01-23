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
ntrees = 2000

cbinom = correlbinom::correlbinom(rho, successprob = 1/ncol(newX), trials=2000, model='kuk')
bin.rf <- binomialRF::binomialRF(newX,y, ntrees = 2000, fdr.threshold = 0.05, user_cbinom_dist = cbinom)

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

test = data.frame(do.call(rbind, lapply(1:50, function(x) compareModels())))
colnames(test) <- c('binomialRF', 'naive RF')


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### 
#### Model Averaging Effects STUDY
####  
#### look at Pathway-based classifiers
#### 
#### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

trueFeatures = colnames(newX)

candidateModels <- list(
  m1=sample(trueFeatures, size=200),
  m2=sample(trueFeatures, size=200),
  m3=sample(trueFeatures, size=300),
  m4=sample(trueFeatures, size=500),

  m5= sample(trueFeatures,  size=500),
  m6= sample(trueFeatures,  size=1000),
  m7= sample(trueFeatures,  size=1000),
  m8= sample(trueFeatures,  size=1500),
  m9= sample(trueFeatures,  size=1500),
  m10=trueFeatures
)

tab = evaluateCandidateModels(candidateModels, newX,y,percent_features = 0.8, ntrees = 3000)
tab$Prop.Selected <- as.numeric(as.character(tab$Prop.Selected))

write.csv(tab, file='~/Dropbox/Samir/binomialRF_study/results/HRV case study/model_averaging_table.csv',row.names = F)

finalFeatures <- as.character(tab$Variable[which(tab$Prop.Selected > 0.5)])

compareModels2 <- function(){
    final.RF2 <- randomForest(newXtest[, names(newXtest) %in% finalFeatures], ytest )
    final.RF2
    OOB.error <- last(final.RF2$err.rate[,1])
    OOB.error

    naiveRF <- randomForest(newXtest, ytest)
    OOB.error.naive <- last(naiveRF$err.rate[,1])
    OOB.error.naive
    return(c(OOB.error,OOB.error.naive))
}

test2 = data.frame(do.call(rbind, lapply(1:50, function(x) compareModels2())))
colnames(test2) <- c('modelAveraging', 'naive RF')
test2
test$modelAveraging <- test2$modelAveraging

test <- test[, c('modelAveraging','binomialRF','naive RF')]

boxplot(test, ylim=c(0,1) )
abline(h=median(test$modelAveraging)[1], lwd=3, lty=2, col='green')

final.mat <- melt(test)

pdf('~/Dropbox/Samir/binomialRF_study/results/HRV case study/binomialRFmodel_averaging.pdf')
ggplot(final.mat) + geom_boxplot(aes(x=variable, y=value)) + ylim(c(0,1)) +
  theme_linedraw()+
  theme(plot.title = element_text(color="black", size=30, face="bold",hjust = 0.5),
        axis.title.x = element_text(color="black", size=40, face="bold"),
        strip.text = element_text(size=25),
        axis.title.y = element_text(color="black", size=20, face="bold"),
        axis.text.x = element_text(color="black", size=20, face="bold", angle = 90))+
  xlab('')  +ylab('Average Validation Error') + ggtitle('HRV Validation Error')

dev.off()

final.mat = data.table(final.mat)
setkey(final.mat, variable)
final.mat.summary =final.mat[, list(MeanValidationError = round(mean(value),3),
                 MedianValidationError = round(median(value),3),
                 SD =round(sd(value),2)), by=variable]

write.csv(final.mat.summary, file='~/Dropbox/Samir/binomialRF_study/results/HRV case study/binomialRFmodel_averaging.csv')
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

k.bin.rf <- k.binomialRF(newX[,trueFeatures],y, ntrees =3000, percent_features = .4, fdr.threshold = 0.05, K=2)
k.bin.rf <- k.bin.rf[Significant==T]
k.bin.rf <- data.frame(k.bin.rf, stringsAsFactors = F)

k.bin.rf <- transform(k.bin.rf, Interaction = colsplit(Interaction, split = " ⊗ ", names = c('pathway1', 'pathway2')))

interaction.matrix <- k.bin.rf$Interaction
interaction.matrix$pathway1 <- as.character(interaction.matrix$pathway1)
interaction.matrix$pathway2 <- as.character(interaction.matrix$pathway2)
colnames(interaction.matrix)[1] <- 'path_id'

go.bp.description <- data.frame(go.bp.description)

interaction.matrix <- plyr::join(interaction.matrix, go.bp.description[,c('path_id','description')])
interaction.matrix <- dplyr::left_join(interaction.matrix, go.bp.description[,c('path_id','description')], c('pathway2'='path_id'))

final.interaction.matrix <- cbind(k.bin.rf, interaction.matrix[, c('description.x','description.y')])
write.csv(final.interaction.matrix, '~/Dropbox/Samir/binomialRF_study/results/HRV case study/pathway_pathway_interactions.csv', row.names = F)
