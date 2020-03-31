##########################################################################################
##########################################################################################
#### load data and libraries
rm(list=ls())
require(binomialRF)
require(randomForest)
load('~/Dropbox/Samir/binomialRF_study/code/rhinoVirus/TCGA_all.RData')
require(TCGA2STAT)
setwd("~/Dropbox/Samir/binomialRF_study/code/TCGA_validation")


res.brca <- read.csv('BRCACancerList_Interactions.csv')
res.kidn <- read.csv('KidneyCancerList_Interactions.csv')

#### 
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

plot_interactions <- function(interaction_string, x.data,y.data){
  
  orange =rgb(0.9,0.52,0.14,0.8)
  green= rgb(0.1,0.8,0.1,0.5)
  
  inter.mat <- x.data[, interaction_string]
  inter.mat$y <- y.data
  
  inter.mat$ratio <- log(inter.mat[,1] / inter.mat[,2])
  
  rangeVals = c(min(inter.mat$ratio),max(inter.mat$ratio))
  
  hist(inter.mat$ratio[inter.mat$y=='normal'], breaks=20, col=orange, xlim=rangeVals, freq=F, main='', 
       xlab='Log(Ratio)',ylab='Density distribution')
  hist(inter.mat$ratio[inter.mat$y=='tumor'], breaks=20, add=T, col=green,  freq=F)
  
  legend('topright', c('Tumor', 'Normal'), col=c(green,orange), lty=1, lwd=10)
  title(paste('Log ratio',paste(interaction_string,collapse = ' | '), 'Interaction'))
  
}

##########################################################################################
##########################################################################################

brca_interactions <- as.character(res.brca$Interaction)
brca_interactions <- do.call(rbind, lapply(1:length(brca_interactions), function(x) unlist(strsplit(brca_interactions[x], ' \\| '))))

data.brca <- prepareData(TumorNormalMatch(brca0$dat))
x.data.brca <- data.brca$X
y.data.brca <- data.brca$y

# pdf('brca_interactions.pdf')
# par(mfrow=c(4,4))
# lapply(1:16, function(x) plot_interactions(brca_interactions[x,], x.data.brca, y.data.brca))
# dev.off()
# 
# ###### 
# ######
# ######
# 
kidn_interactions <- as.character(res.kidn$Interaction)
kidn_interactions <- do.call(rbind, lapply(1:length(kidn_interactions), function(x) unlist(strsplit(kidn_interactions[x], ' \\| '))))

data.kidn <- prepareData(TumorNormalMatch(kipan0$dat))
x.data.kidn <- data.kidn$X
y.data.kidn <- data.kidn$y
# 
# pdf('kipan_interactions.pdf')
# par(mfrow=c(4,3))
# lapply(1:11, function(x) plot_interactions(kidn_interactions[x,], x.data.kidn, y.data.kidn))
# dev.off()

##########################################################################################
##########################################################################################

##########################################################################################
##########################################################################################

#### Create plots for interaction figure 

#spyr col ---- brca
orange =rgb(0.9,0.52,0.14,0.8)
green= rgb(0.1,0.8,0.1,0.5)

pdf('brca-hist1.pdf')
hist(log(x.data.brca[y.data.brca=='tumor', 'SPRY2']), col =green, xlim=c(-1, 4), 
     freq=F , xlab='', ylab='', main='SPRY2 log(Expression)',breaks=20, ylim=c(0,3))
hist(log(x.data.brca[y.data.brca=='normal', 'SPRY2']), col = orange,add =T, freq=F, breaks = 20)
legend('topright', c('Tumor', 'Normal'), col=c(green,orange), lty=1, lwd=10)
dev.off()

pdf('brca-hist2.pdf')
hist(log(x.data.brca[y.data.brca=='tumor', 'COL10A1']), col =green, xlim=c(-6,6), freq=F , xlab='', ylab='', 
     main='COL10A1 log(Expression)', breaks=20)
hist(log(x.data.brca[y.data.brca=='normal', 'COL10A1']), col = orange,add =T, freq=F, breaks=20)
legend('topright', c('Tumor', 'Normal'), col=c(green,orange), lty=1, lwd=10)
dev.off()

pdf('brca-hist3.pdf')
plot_interactions(brca_interactions[10,], x.data.brca, y.data.brca)
dev.off()



# tfa2pa sgpp1 --- kidney

orange =rgb(0.9,0.52,0.14,0.8)
green= rgb(0.1,0.8,0.1,0.5)

pdf('kidn-hist1.pdf')
hist(log(x.data.kidn[y.data.kidn=='tumor', 'TFAP2A']), col =green, xlim=c(-4,4), 
     breaks=20, freq=F , xlab='',ylab='', main='TFAP2A log(Expression)', ylim=c(0,2))
hist(log(x.data.kidn[y.data.kidn=='normal', 'TFAP2A']), col = orange,add =T, freq=F, breaks=20)
legend('topright', c('Tumor', 'Normal'), col=c(green,orange), lty=1, lwd=10)
dev.off()

pdf('kidn-hist2.pdf')
hist(log(x.data.kidn[y.data.kidn=='tumor', 'SGPP1']), col =green, breaks=20,xlim=c(0,5),
     freq=F , xlab='',ylab='', main='SGPP1 log(Expression)', ylim=c(0,2.5))
hist(log(x.data.kidn[y.data.kidn=='normal', 'SGPP1']), col = orange,add =T, breaks=20, freq=F)
legend('topright', c('Tumor', 'Normal'), col=c(green,orange), lty=1, lwd=10)
dev.off()


pdf('kidn-hist3.pdf')
plot_interactions(kidn_interactions[1,], x.data.kidn, y.data.kidn)
dev.off()





