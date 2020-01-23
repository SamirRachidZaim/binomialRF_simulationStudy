require(data.table)
require(ggplot2)

res = fread('~/Dropbox/Samir/binomialRF_study/code/InteractionSimulation.csv')
res_kbinomRF <- fread('~/Dropbox/Samir/binomialRF_study/results/binomialRF1000_interaction.csv')
res = res[res$Model != 'binomialRF']
res = rbind(res, res_kbinomRF)

res = res[, -1]
res2 = melt(res)
res2$value = round(res2$value,2)
res2$value[is.na(res2$value)] <- 0

ggplot(res2[res2$variable=='Precision']) + geom_boxplot(aes(x=reorder(Model,value,mean), y=value))
ggplot(res2[res2$variable=='Recall']) + geom_boxplot(aes(x=reorder(Model,value,mean), y=value))
ggplot(res2[res2$variable=='ModelSize']) + geom_boxplot(aes(x=reorder(Model,value,mean), y=value))

res3 = res2[, list(round(mean(value),2), round(sd(value),2)), by=list(Model, variable)]
m1 = dcast(res3, formula = Model~variable, value.var = 'V1')
m2 = dcast(res3, formula = Model~variable, value.var = 'V2')

m1$Precision = paste( m1$Precision, ' (', m2$Precision, ')',sep='')
m1$Recall = paste( m1$Recall, ' (', m2$Recall, ')',sep='')
m1$ModelSize = paste( m1$ModelSize, ' (', m2$ModelSize, ')',sep='')

m1
write.csv(m1, file='~/Dropbox/Samir/binomialRF_study/results/InteractionFinalTable.csv', row.names = F)
30