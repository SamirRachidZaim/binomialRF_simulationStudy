# binomialRF virogram 
require(data.table)
require(randomForest)
require(foreign)
# training.rhv <-  fread('~/Dropbox/Virogram/Java/ClassificationVirogram/HRV/summary.binary.signed.txt', data.table = F); row.names(rhv) <- rhv$V1; rhv= rhv[,-1]
# 



as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}


# asthma.validation. <-  fread('~/Dropbox/Virogram/Java/ClassificationVirogram/Ex Vivo/summary.binary.signed.txt', data.table = F); row.names(rhv) <- rhv$V1; rhv= rhv[,-1]
# flu <-  fread('~/Dropbox/Virogram/Java/ClassificationVirogram/FLU/summary.binary.signed.txt', data.table = F); row.names(flu) <- flu$V1; flu= flu[,-1]
rhv.pheno <- foreign::read.arff('~/Dropbox/Virogram/Java/ClassificationVirogram/HRV/summary.binary.signed.txt.arff')
exvivo.pheno <- foreign::read.arff('~/Dropbox/Virogram/Java/ClassificationVirogram/Ex Vivo/summary.binary.signed.txt.arff')
go.bp.description <- fread("~/Dropbox/Virogram/Pathways/go.bp.description.txt")
go.bp.description$Pathway <- go.bp.description$path_id


X = rhv.pheno[,! names(rhv.pheno) %in% c('Class', 'Name')]
y = rhv.pheno[, names(rhv.pheno) %in% 'Class']
Xtest <- exvivo.pheno[, !names(exvivo.pheno) %in% c("Class",'Name')]
ytest <-  exvivo.pheno[, names(exvivo.pheno) %in% c("Class")]

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

require(binomialRF)
devtools::load_all('~/Dropbox/Samir/binomialRF')

get.table6 <- function(X,y, s=.8){
    bin.rf <- binomialRF(X,y, ntrees = 1000, percent_features = s, fdr.threshold = 0.05)
    true.pathways <- fread('~/Dropbox/Samir/binomialRF_study/Table2_pathways.csv')
    
    rf.object = randomForest(X[,bin.rf$Variable[bin.rf$Significant==T]],y)
    
    capt <- sum(true.pathways$Pathway %in% bin.rf$Variable[bin.rf$Significant==T] ) / length(true.pathways$Pathway)
    
    predicted.features <- data.frame(Pathway = bin.rf$Variable[bin.rf$Significant==T],
                                     ValidatedInHRV = ifelse(bin.rf$Variable[bin.rf$Significant==T] %in% true.pathways$Pathway, "Yes","No"),
                                      stringsAsFactors = F)
    
    
    predicted.features <-  plyr::join(predicted.features, go.bp.description[, c('Pathway','description')])
    predicted.features <-  plyr::join(predicted.features, true.pathways[, c('Pathway','Class')])
    size <- nrow(predicted.features)
    
    write.csv(predicted.features, file='~/Dropbox/Samir/binomialRF_study/results/HRV case study/predicted.pathways.csv',row.names = F)
    
    return(data.frame( s = s,
                       Predicted_Pathways = size,
                       ConfirmedPathways=capt))

}

final.table <- 
  rbind(
    get.table6(X,y, .3),
    get.table6(X,y, .6),
    get.table6(X,y, .8)
    )
final.table

write.csv(final.table, file='~/Dropbox/Samir/binomialRF_study/results/HRV case study/pathway_analysis_percent_features.csv',row.names = F)


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### 
#### Model Averaging STUDY
####  
#### 
#### 
#### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
bin.rf <- binomialRF(X,y, ntrees = 1000, percent_features = 0.3 , fdr.threshold = 0.05)

trueFeatures =bin.rf$Variable[bin.rf$Significant==T]

candidateModels <- list(
  m1=sample(trueFeatures, size=20),
  m2=sample(trueFeatures, size=20),
  m3=sample(trueFeatures, size=30),
  m4=sample(trueFeatures, size=30),
  
  m5=sample(trueFeatures, size=40),
  m6= sample(trueFeatures, size=40),
  m7=sample(trueFeatures, size=40),
  m8= sample(trueFeatures, size=40),
  m9=sample(trueFeatures, size=40),
  m10= trueFeatures
)

tab = evaluateCandidateModels(candidateModels, X,y,percent_features = 0.3, ntrees = 5000)
tab$Prop.Selected <- as.numeric(as.character(tab$Prop.Selected))

finalFeatures <- as.character(tab$Variable[which(tab$Prop.Selected > 0)])

final.RF <- randomForest(X[, trueFeatures], y , mtry = 9)
final.RF2 <- randomForest(Xtest[, names(Xtest) %in% trueFeatures], ytest , mtry = 15)
final.RF2







# for (p in finalFeatures) { 
#   if (class(X[[p]]) == "factor") { 
#     levels(Xtest[[p]]) <- levels(X[[p]]) 
#   } 
# }
# 
# preds <- ifelse(as.character(predict(final.RF, newdata = Xtest) )=='symptomatic', 1,0)
# truth <- ifelse(as.character(ytest)=='RE', 1,0)
# 
# cbind(preds,truth)
# mean(preds==truth)


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### 
#### INTERACTION STUDY
####  
#### look at Pathway-Pathway Interaction 
#### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

k.bin.rf <- k.binomialRF(X[,trueFeatures],y, ntrees =1000, percent_features = .4, fdr.threshold = 0.05, K=2)
k.bin.rf[Significant==T]


