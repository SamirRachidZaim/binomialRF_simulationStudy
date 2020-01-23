### plot false selection & discovery rates
require(data.table)
require(ggplot2)
require(ggpubr)
require(ggpubr)

calculateFSR <- function(modelVars,trueVars, totalVars=10){
  
  Uninformative <- sum(!modelVars %in% trueVars)
  #selected <- length(modelVars) +1
  fsr <- Uninformative/10
  return(fsr)
  
}

calculateDiscoveryRate <- function(modelVars,trueVars){
  
  selected <- sum(modelVars %in% trueVars)
  NTrueVars <- length(trueVars) 
  discoveryRate <- selected/(NTrueVars+1)
  return(discoveryRate)
  
}

load('~/results/results/simulation_FDR_BY_1000_trees_indpnFALSE.RData'); m1 = modelList
load('~/results/results/simulation_FDR_BY_1000_trees_indpnTRUE.RData'); m2 = modelList
load('~/results/results/simulation_FDR_BY_5000_trees_indpnFALSE.RData'); m3= modelList
load('~/results/results/simulation_FDR_BY_5000_trees_indpnTRUE.RData'); m4=modelList
load('~/results/results/simulation_FDR_BY_10000_trees_indpnFALSE.RData'); m5= modelList
load('~/results/results/simulation_FDR_BY_10000_trees_indpnTRUE.RData'); m6=modelList

m1.vars = paste('X', 1:5, sep='')

format_data_na<- function(z){
  if(class(z)=='try-error'){
  } else{
    return(z)
  }
}

m3= plyr::compact(lapply(m3, function(z) format_data_na(z)))
m4= plyr::compact(lapply(m4, function(z) format_data_na(z)))
m5= plyr::compact(lapply(m5, function(z) format_data_na(z)))
m6= plyr::compact(lapply(m6, function(z) format_data_na(z)))


get_fsr_plots <- function(m1, m1.vars){
  fsr.mat = data.table(
    binomialRF = sapply(1:length(m1), function(x) calculateFSR(m1[[x]][[1]][1]['SMART_VIS'][[1]], m1.vars ) ),
    #Lasso = sapply(1:length(m1),     function(x) calculateFSR(names(m1[[x]][[1]][3]['LASSO'][[1]]), m1.vars ) ),
    #RF_Soil = sapply(1:length(m1),    function(x) calculateFSR(m1[[x]][[1]][4]['RFSoil'][[1]], m1.vars ) ),
    VSURF = sapply(1:length(m1),     function(x) calculateFSR(paste('X',m1[[x]][[1]][5]['VSURF'][[1]],sep=''), m1.vars ) ),
    Boruta = sapply(1:length(m1),    function(x) calculateFSR(m1[[x]][[1]][6]['Boruta'][[1]], m1.vars ) ),
    VarSelRF = sapply(1:length(m1),  function(x) calculateFSR(m1[[x]][[1]][7]['VarSelRF'][[1]], m1.vars ) ),
    PIMP = sapply(1:length(m1),      function(x) calculateFSR(m1[[x]][[1]][8]['PIMP'][[1]], m1.vars ) )
  )
  
  fsr.mat = melt(fsr.mat)
  
  p <- ggplot(fsr.mat) + geom_boxplot(aes(x=reorder(variable, value, FUN=median), y=value)) + xlab("Method")+
    ylab("False Discovery Rate") + ggtitle("False Discovery Rate by Method")+ 
    theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    axis.text.x = element_text(color="black", size=14, face="bold", angle = 90))+ ylim(c(0,1))

  return(list(p, fsr.mat))
}

p1 = get_fsr_plots(m3, m1.vars )
p2 = get_fsr_plots(m4, m1.vars = '' )
p3 = get_fsr_plots(m5, m1.vars )
p4 = get_fsr_plots(m6, m1.vars = '' )




# pdf('~/Dropbox/Samir/SMART_VIS/plots/sim1.fsr.pdf',
#     width = 11, height = 8)
# 
# 
# figure= ggarrange(p1[[1]], p2[[1]],p3[[1]], p4[[1]],
#                   labels = 'AUTO',
#                   ncol=2, nrow=2)
# annotate_figure(figure,
#                 top = text_grob("False Selection Rate", 
#                                 color = "black", face = "bold", size = 28),
#                 bottom = text_grob("", color = "blue",
#                                    hjust = 1, x = 1, face = "italic", size = 10),
#                 left = text_grob("", color = "green", rot = 90),
#                 right = "",
#                 fig.lab = "", fig.lab.face = "bold"
# )
# dev.off()


pdf('~/Dropbox/Samir/SMART_VIS/plots/sim1.fsr.pdf',
    width = 11, height = 8)

p3[[1]]


dev.off()

get_dscvr_plots <- function(m1, m1.vars){
  dscvr.mat = data.table(
    binomRF = sapply(1:length(m1), function(x) calculateDiscoveryRate(m1[[x]][[1]][1]['SMART_VIS'][[1]], m1.vars ) ),
   # Lasso = sapply(1:length(m1),     function(x) calculateDiscoveryRate(names(m1[[x]][[1]][3]['LASSO'][[1]]), m1.vars ) ),
    RF_Soil = sapply(1:length(m1),    function(x) calculateDiscoveryRate(m1[[x]][[1]][4]['RFSoil'][[1]], m1.vars ) ),
    VSURF = sapply(1:length(m1),     function(x) calculateDiscoveryRate(paste('X',m1[[x]][[1]][5]['VSURF'][[1]],sep=''), m1.vars ) ),
    Boruta = sapply(1:length(m1),    function(x) calculateDiscoveryRate(m1[[x]][[1]][6]['Boruta'][[1]], m1.vars ) ),
    VarSelRF = sapply(1:length(m1),  function(x) calculateDiscoveryRate(m1[[x]][[1]][7]['VarSelRF'][[1]], m1.vars ) ),
    PIMP = sapply(1:length(m1),      function(x) calculateDiscoveryRate(m1[[x]][[1]][8]['PIMP'][[1]], m1.vars ) )
  )
  
  dscvr.mat = melt(dscvr.mat)
  
  p <- ggplot(dscvr.mat) + geom_boxplot(aes(x=variable, y=value)) + ylim(c(0,1))
  return(list(p, dscvr.mat))
}

n1= get_dscvr_plots(m3, m1.vars )
n2= get_dscvr_plots(m4, m1.vars = '' )
n3= get_dscvr_plots(m5, m1.vars )
n4= get_dscvr_plots(m6, m1.vars = '' )

pdf('~/Dropbox/Samir/SMART_VIS/plots/sim1.dscvr.pdf',
    width = 11, height = 8)


figure= ggarrange(n1[[1]], n2[[1]],n3[[1]],n4[[1]],
                  labels = 'AUTO',
                  ncol=2, nrow=2)

annotate_figure(figure,
                top = text_grob("True Coverage/Discovery Rate", 
                                color = "black", face = "bold", size = 28),
                bottom = text_grob("", color = "blue",
                                   hjust = 1, x = 1, face = "italic", size = 10),
                left = text_grob("", color = "green", rot = 90),
                right = "",
                fig.lab = "", fig.lab.face = "bold"
)
dev.off()





