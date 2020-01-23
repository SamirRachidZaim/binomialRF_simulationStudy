#plot time results 
require(data.table)
require(ggplot2)
require(parallel)

sim10 = fread('~/Dropbox/Samir/binomialRF_study/results/timeProfile_10_features.csv')
sim100 = fread('~/Dropbox/Samir/binomialRF_study/results/timeProfile_100_features.csv')
sim1000 = fread('~/Dropbox/Samir/binomialRF_study/results/timeProfile_1000_features.csv')

sim10Interaction = fread('~/Dropbox/Samir/binomialRF_study/results/interactions_TimeProfile_10_features.csv')

removeOutliers <- function(sim10){
  
  new.mat = data.frame()
  
  for(model in unique(sim10$Model)){
    model.mat = sim10[sim10$Model==model][order(sim10[Model == model]$elapsed)]
    new.mat = rbind(new.mat,model.mat )
  }
  
  new.mat <- new.mat[, c('Model','elapsed')]
  new.mat <- new.mat[ ! grep('Error', new.mat$Model),]
  new.mat$elapsed <- as.numeric(new.mat$elapsed)
  
  return(new.mat)
}

sim10 = removeOutliers(sim10)
sim100 = removeOutliers(sim100)
sim1000 = removeOutliers(sim1000)
sim10Interaction = removeOutliers(sim10Interaction)

generate_boxplots <- function(sim100, numFeatures=100){

  ggplot(sim100 ) + geom_boxplot(aes(x=reorder(Model, elapsed, FUN=median), y=log10(elapsed+1)))  +  
    theme_linedraw()+ labs(title=paste(numFeatures, 'Features'), x='',y='Log10(Time (seconds))')+
    theme(
      plot.title = element_text(color="black", size=16, face="bold.italic"),
      axis.title.x = element_text(color="black", size=24, face="bold"),
      axis.title.y = element_text(color="black", size=20, face="bold"),
      axis.text.x = element_text(color="black", size=24, face="bold", angle = 90))+
    ylim(c(0, max(sim100$elapsed)['90%'])) + ylim(c(0,3))
  
}

m1 = generate_boxplots(sim10, numFeatures = 10)
m2 = generate_boxplots(sim100, numFeatures = 100)
m3 = generate_boxplots(sim1000, numFeatures = 1000)

m1
m2
m3

m1Inter = generate_boxplots(sim10Interaction, 10)

require(ggpubr)

pdf('~/Dropbox/BiniomialRF/Bioinformatics submission/Figures/tmp_Figure6_timeProfile.pdf', width = 15, height = 9)
figure= ggarrange(m1,m2,m3,
                  ncol=3, nrow=1)
annotate_figure(figure,
                top = text_grob("Computation Time Across Methods", 
                                color = "black", face = "bold", size = 24),
                bottom = text_grob("", color = "blue",
                                   hjust = 1, x = 1, face = "italic", size = 10),
                left = text_grob("", color = "green", rot = 90),
                right = "",
                fig.lab = "", fig.lab.face = "bold"
)

dev.off()



#### create time table 

X= sim10[, list(Avg.Time = mean(elapsed), Std.Dev =sd(elapsed)), by=Model]
X2=sim100[, list(Avg.Time = mean(elapsed), Std.Dev =sd(elapsed)), by=Model]
X3=sim1000[, list(Avg.Time = mean(elapsed), Std.Dev =sd(elapsed)), by=Model]

fin.mat = cbind(X,X2[,2:3], X3[, 2:3])
colnames(fin.mat) <- c('Model', 'Avg.Time (X10)', 'Std.Dev (X10)', 'Avg.Time (X100)','Std.Dev (X100)', 'Avg.Time (X1000)','Std.Dev (X1000)') 
fin.mat[, 2:7] <- round(fin.mat[,2:7],2)


fin.mat2 = data.frame(cbind(
  paste(fin.mat$`Avg.Time (X10)`,' (', fin.mat$`Std.Dev (X10)`, ')', sep = ''),
  paste(fin.mat$`Avg.Time (X100)`,' (', fin.mat$`Std.Dev (X100)`, ')', sep = ''),
  paste(fin.mat$`Avg.Time (X1000)`,' (', fin.mat$`Std.Dev (X1000)`, ')', sep = '')),
  row.names = fin.mat$Model); colnames(fin.mat2) <- c('10 features','100 features','1000 features')
fin.mat2 = fin.mat2[order(fin.mat$`Avg.Time (X10)`),]

write.csv(fin.mat2, '~/Dropbox/Samir/binomialRF_study/results/final.time.table.csv')
