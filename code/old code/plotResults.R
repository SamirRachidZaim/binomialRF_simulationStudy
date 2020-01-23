#plot simulation results error
require(data.table)
require(ggplot2)

load('~/results/results/simulation_FDR_BY_1000_trees_indpnFALSE.RData'); m1 = modelList
load('~/results/results/simulation_FDR_BY_1000_trees_indpnTRUE.RData'); m2 = modelList
load('~/results/results/simulation_FDR_BY_5000_trees_indpnFALSE.RData'); m3= modelList
load('~/results/results/simulation_FDR_BY_5000_trees_indpnTRUE.RData'); m4=modelList
load('~/results/results/simulation_FDR_BY_10000_trees_indpnFALSE.RData'); m5= modelList
load('~/results/results/simulation_FDR_BY_10000_trees_indpnTRUE.RData'); m6=modelList

format.data <- function(m1, type='TrainError'){
  m1.Train = lapply(1:1000, function(x) m1[[x]][type]) 
  m1.Train = lapply(1:1000, function(x) if(class(m1.Train[x][[1]]) =='list'){ m1.Train[x][[1]][[1]]})
  m1.Train = data.table(do.call(rbind, m1.Train))
  m1.Train =  melt(m1.Train)
  m1.Train$model = type
  m1.Train$variable = as.character(m1.Train$variable)
  return(m1.Train)
}

generate_boxplots <- function(m1, y.label = 'Independent Y'){
  m1.Test = format.data(m1, "TestError")
  m1.Train = format.data(m1, "TrainError")
  m1.final = rbind(m1.Test, m1.Train)
  m1.final$variable[m1.final$variable=='VarselRF'] <- 'VarSelRF'
  m1.final$variable[m1.final$variable=='SMART_VIS'] <- 'binomialRF'

  m1.final<- m1.final[!m1.final$variable%in% c('LASSO','Gini','RFSoil'),] 
  
  m1.final$model <- factor(m1.final$model, levels = c('TrainError','TestError'))
  
  p1 <- ggplot(m1.final ) + geom_boxplot(aes(x=reorder(variable, value, FUN=median), y=value)) + facet_grid(.~model) +  
    theme_dark()+ labs( x='Model',y=paste("Error:", y.label))+ggtitle('0-1 Error across Different Methods')+
  
    theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
          axis.title.x = element_text(color="black", size=14, face="bold"),
          axis.title.y = element_text(color="black", size=14, face="bold"),
          axis.text.x = element_text(color="black", size=14, face="bold", angle = 90))+ ylim(c(0,1))
  
  
  return(list(p1, m1.final))
}
# 
# p1 = generate_boxplots(m1)
# p2 = generate_boxplots(m2)
p3 = generate_boxplots(m3, y.label = '0-1 Error')
p4 = generate_boxplots(m4,y.label = '0-1 Error')
p5 = generate_boxplots(m5,y.label = '0-1 Error')
p6 = generate_boxplots(m6,y.label = '0-1 Error')

require(ggpubr)
figure= ggarrange(p3[[1]],p5[[1]] + rremove('y.title'),
          ncol=2, nrow=1)

pdf('~/Dropbox/Samir/SMART_VIS/plots/sim1.errors.pdf',
    width = 11, height = 8)
p3[[1]]
dev.off()
# figure2= ggarrange(p2,p4,p6, labels=c('500 Trees','1000 Trees','10000 Trees'),
#           ncol=3, nrow=1)
# 
# annotate_figure(figure2,
#                 top = text_grob("0-1 Error Across Deep Forests: Independent Y", color = "red", face = "bold", size = 14),
#                 bottom = text_grob("Data source: \n Simulation", color = "blue",
#                                    hjust = 1, x = 1, face = "italic", size = 10),
#                 left = text_grob("Figure arranged using ggpubr", color = "green", rot = 90),
#                 right = "",
#                 fig.lab = "Figure 1", fig.lab.face = "bold"
# )



# m1 = p1[[2]]
# m2 = p2[[2]]
m3 = p3[[2]]
m4 = p4[[2]]
m5 = p5[[2]]
m6 = p6[[2]]

m3 = m3[, list(Avg.Error=paste(round(mean(value),3),'(', round(sd(value),2),')')),
          by= list(variable, model)]
m3 = dcast(m3, formula = variable ~ model)

m4 = m4[, list(Avg.Error=paste(round(mean(value),3),'(', round(sd(value),2),')')),
   by= list(variable, model)]
m4 = dcast(m4, formula = variable ~ model)

m5 = m5[, list(Avg.Error=paste(round(mean(value),3),'(', round(sd(value),2),')')),
   by= list(variable, model)]
m5 = dcast(m5, formula = variable ~ model)

m6 = m6[, list(Avg.Error=paste(round(mean(value),3),'(', round(sd(value),2),')')),
   by= list(variable, model)]
m6 = dcast(m6, formula = variable ~ model)

dependent.y = cbind(m3, m5[, 2:3])
independent.y= cbind(m4, m6[, 2:3])

colnames(dependent.y) <- c('variable','TrainError (5k)','TestError (5k)',
                           'TrainError (10k)', 'TestError (10k)')

dependent.y = dependent.y[order(dependent.y$`TestError (10k)`),]
write.csv(dependent.y, file = '~/Dropbox/Samir/SMART_VIS/results/10_features_error_dpn.y.csv', row.names=F)

colnames(independent.y) <- c('variable','TrainError (5k)','TestError (5k)',
                           'TrainError (10k)', 'TestError (10k)')

independent.y= independent.y[order(independent.y$`TestError (10k)`),]
write.csv(independent.y, file = '~/Dropbox/Samir/SMART_VIS/results/10_features_error_indp.y.csv', row.names=F)



