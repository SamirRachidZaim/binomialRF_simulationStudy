## createPltos
require(data.table)
require(ggplot2)

setwd('~/Dropbox/Samir/binomialRF_study/results')

dim.X = 132000
ntrees=1000

create_plots <- function(dim.X){
  m100_5 <-  fread(paste('accuracy_1000Trees_simulation_',dim.X,'_P_5_genes.csv', sep=''), header = T,stringsAsFactors = F); m100_5$NumGenes <- 5
  m100_25 <- fread(paste('accuracy_1000Trees_simulation_',dim.X,'_P_25_genes.csv',sep=''), header = T, stringsAsFactors = F); m100_25$NumGenes <- 25
  m100_50 <- fread(paste('accuracy_1000Trees_simulation_',dim.X,'_P_50_genes.csv',sep=''), header = T, stringsAsFactors = F); m100_50$NumGenes <- 50
  
  if(dim.X > 100){
    m100_100 <- fread(paste('accuracy_1000Trees_simulation_',dim.X,'_P_100_genes.csv',sep=''), header = T, stringsAsFactors = F); m100_100$NumGenes <- 100
    m100 <- rbind(m100_5,m100_25, m100_50, m100_100)
  } else{
    m100 <- rbind(m100_5,m100_25, m100_50)
  }
  
    
  
  
  m100 <- cbind(as.data.frame(sapply(m100[, c(1:5,7)], as.numeric)), m100$Model)
  m100 = m100[, 2:7]
  m100 = data.table::melt(m100, id.vars=c('m100$Model','NumGenes'))
  colnames(m100)[1] <- 'Model'
  m100 = m100[!is.na(m100$value),]
  m100$Model = as.character(m100$Model)
  m100$variable = as.character(m100$variable)
  
  pdf(paste('Precision_',dim.X,'.pdf', sep=''), width = 11, height = 8)
  a = ggplot(m100[m100$variable=='Precision',]) + geom_boxplot(aes(x=reorder(Model, value, FUN=median), y = value)) + facet_wrap(~NumGenes, scales='free') + theme_linedraw()+
    theme(plot.title = element_text(color="black", size=14, face="bold.italic"),
          axis.title.x = element_text(color="black", size=24, face="bold"),
          axis.title.y = element_text(color="black", size=14, face="bold"),
          axis.text.y = element_text(color="black", size=24, face="bold"),
          
          axis.text.x = element_text(color="black", size=24, face="bold", angle = 90),
          plot.margin = margin(2,2,2,2, 'cm'))+ xlab('') + ylab('') + 
    ggtitle('Variable Precision')
  
  a = ggplot(m100[m100$variable=='Precision',]) + geom_boxplot(aes(x=reorder(Model, value, FUN=median), y = value, color=NumGenes)) +
    #facet_wrap(~NumGenes, scales='free') + theme_linedraw()+
    theme(plot.title = element_text(color="black", size=14, face="bold.italic"),
          axis.title.x = element_text(color="black", size=24, face="bold"),
          axis.title.y = element_text(color="black", size=14, face="bold"),
          axis.text.y = element_text(color="black", size=24, face="bold"),
          
          axis.text.x = element_text(color="black", size=24, face="bold", angle = 90),
          plot.margin = margin(2,2,2,2, 'cm'))+ xlab('') + ylab('') + 
    ggtitle('Variable Precision')
  print(a)
  dev.off()
  
  pdf(paste('Recall_',dim.X,'.pdf', sep=''), width = 11, height = 8)
  b = ggplot(m100[m100$variable=='Recall',]) + geom_boxplot(aes(x=reorder(Model, value, FUN=median), y = value)) + facet_wrap(~NumGenes, scales='free') + theme_linedraw()+
    theme(plot.title = element_text(color="black", size=14, face="bold.italic"),
          axis.title.x = element_text(color="black", size=14, face="bold"),
          axis.title.y = element_text(color="black", size=14, face="bold"),
          axis.text.y = element_text(color="black", size=24, face="bold"),
          axis.text.x = element_text(color="black", size=14, face="bold", angle = 90),
          plot.margin = margin(2,2,2,2, 'cm'))+ xlab('') + ylab('')+ ggtitle('Variable Recall')
  
  b = ggplot(m100[m100$variable=='Recall',]) + geom_boxplot(aes(x=reorder(Model, value, FUN=median), y = value)) + 
    #facet_wrap(~NumGenes, scales='free') + theme_linedraw()+
    theme(plot.title = element_text(color="black", size=14, face="bold.italic"),
          axis.title.x = element_text(color="black", size=24, face="bold"),
          axis.title.y = element_text(color="black", size=14, face="bold"),
          axis.text.y = element_text(color="black", size=24, face="bold"),
          axis.text.x = element_text(color="black", size=24, face="bold", angle = 90),
          plot.margin = margin(2,2,2,2, 'cm'))+ xlab('') + ylab('')+ ggtitle('Variable Recall')
  
  print(b)
  dev.off()
  
  pdf(paste('TestError',dim.X,'.pdf', sep=''), width = 11, height = 8)
  c = ggplot(m100[m100$variable=='TestError',]) + geom_boxplot(aes(x=reorder(Model, value, FUN=median), y = value)) + facet_wrap(~NumGenes, scales='free') + theme_linedraw()+
    theme(plot.title = element_text(color="black", size=14, face="bold.italic"),
          axis.title.x = element_text(color="black", size=24, face="bold"),
          axis.title.y = element_text(color="black", size=24, face="bold"),
          axis.text.x = element_text(color="black", size=24, face="bold", angle = 90),
          plot.margin = margin(2,2,2,2, 'cm'))+ xlab('') + ylab('')+
    ggtitle('Test Error')
  
  c = ggplot(m100[m100$variable=='TestError',]) + geom_boxplot(aes(x=reorder(Model, value, FUN=median), y = value)) +
    #facet_wrap(~NumGenes, scales='free') + theme_linedraw()+
    theme(plot.title = element_text(color="black", size=14, face="bold.italic"),
          axis.title.x = element_text(color="black", size=24, face="bold"),
          axis.title.y = element_text(color="black", size=14, face="bold"),
          axis.text.y = element_text(color="black", size=24, face="bold"),
          
          axis.text.x = element_text(color="black", size=24, face="bold", angle = 90),
          plot.margin = margin(2,2,2,2, 'cm'))+ xlab('') + ylab('')+ 
    ggtitle('Test Error')
  
  print(c)
  dev.off()
  
  pdf(paste('ModelSize',dim.X,'.pdf', sep=''), width = 11, height = 8)
  d = ggplot(m100[m100$variable=='ModelSize',]) + geom_boxplot(aes(x=reorder(Model, value, FUN=median), y = value)) + facet_wrap(~NumGenes, scales='free') + theme_linedraw()+
    theme(plot.title = element_text(color="black", size=14, face="bold.italic"),
          axis.title.x = element_text(color="black", size=14, face="bold"),
          axis.title.y = element_text(color="black", size=14, face="bold"),
          axis.text.x = element_text(color="black", size=14, face="bold", angle = 90),
          plot.margin = margin(2,2,2,2, 'cm'))+ xlab('') + ylab('')+ 
    ggtitle('Model Size')
  
  d = ggplot(m100[m100$variable=='ModelSize',]) + geom_boxplot(aes(x=reorder(Model, value, FUN=median), y = value)) + 
    #facet_wrap(~NumGenes, scales='free') + theme_linedraw()+
    theme(plot.title = element_text(color="black", size=14, face="bold.italic"),
          axis.title.x = element_text(color="black", size=24, face="bold"),
          axis.title.y = element_text(color="black", size=14, face="bold"),
          axis.text.y = element_text(color="black", size=24, face="bold"),
          axis.text.x = element_text(color="black", size=24, face="bold", angle = 90),
          
          plot.margin = margin(2,2,2,2, 'cm'))+ xlab('') + ylab('')+ 
    ggtitle('Model Size')
  
  print(d)
  dev.off()
  
  
}

create_plots(100)
create_plots(500)
create_plots(1000)
create_plots(2000)
create_plots(5000)


create_summary_tables <- function(dim.X){
  
  m100_5 <-  fread(paste('accuracy_1000Trees_simulation_',dim.X,'_P_5_genes.csv', sep=''), header = T,stringsAsFactors = F); m100_5$NumGenes <- 5
  m100_25 <- fread(paste('accuracy_1000Trees_simulation_',dim.X,'_P_25_genes.csv',sep=''), header = T, stringsAsFactors = F); m100_25$NumGenes <- 25
  m100_50 <- fread(paste('accuracy_1000Trees_simulation_',dim.X,'_P_50_genes.csv',sep=''), header = T, stringsAsFactors = F); m100_50$NumGenes <- 50
  
  if(dim.X > 100){
    m100_100 <- fread(paste('accuracy_1000Trees_simulation_',dim.X,'_P_100_genes.csv',sep=''), header = T, stringsAsFactors = F); m100_100$NumGenes <- 100
    m100 <- rbind(m100_5,m100_25, m100_50, m100_100)
  } else{
    m100 <- rbind(m100_5,m100_25, m100_50)
  }
  
  
  m100 <- cbind(as.data.frame(sapply(m100[, c(1:5,7)], as.numeric)), m100$Model)
  m100 = m100[, 2:7]
  m100 = data.table::melt(m100, id.vars=c('m100$Model','NumGenes'))
  colnames(m100)[1] <- 'Model'
  m100 = m100[!is.na(m100$value),]
  m100$Model = as.character(m100$Model)
  m100$variable = as.character(m100$variable)
  
  m100 = data.table(m100)
  summaryRes = m100[,  list(round(mean(value),2),round(sd (value),2)), by=list(Model, variable)] 
  summaryRes$Value = paste(summaryRes$V1, ' (',summaryRes$V2, ')', sep='')
  summaryRes = summaryRes[,c(1,2,5)]
  dcast(summaryRes, formula = Model~variable)
}

sum100  = create_summary_tables(100)
sum500  = create_summary_tables(500)
sum1000 =create_summary_tables(1000)
sum1000 =create_summary_tables(2000)


write.csv(sum100, file = '../agg100Results.csv', row.names = F)
write.csv(sum500, file = '../agg500Results.csv', row.names = F)
write.csv(sum1000, file = '../agg1000Results.csv', row.names = F)




create_summary_of_all_results <- function(){
  
  fnames= dir()
  idx = grep('accuracy', dir())
  
  files = lapply(fnames[idx], function(x) fread(x, header = T,stringsAsFactors = F))
  files = do.call(rbind, files)
  
  files = data.table(files)
  files = files[,c(2:4,6)]
  files = data.table::melt(files, id.vars=c('Model'))
  files$value = as.numeric(files$value)
  files = files[!is.na(files$value)]
  
  
  summaryRes = files[,  list(round(mean(value),2),round(sd (value),2)), by=list(Model, variable)] 
  summaryRes$Value = paste(summaryRes$V1, ' (',summaryRes$V2, ')', sep='')
  summaryRes = summaryRes[,c(1,2,5)]
  return(dcast(summaryRes, formula = Model~variable))
}
final_table <- create_summary_of_all_results()

write.csv(final_table, file = '../fullAggregatedResults.csv', row.names = F)
