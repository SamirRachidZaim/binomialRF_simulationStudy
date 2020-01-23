#### Variable Importance Helper Functions 

########################################################
########################################################
### Run Random forest
###
### input
### @ rf = random forest object 
###
### output
### @ sortedImpMatrix = sorted Gini Variable Importance Matrix
getImportantRF <- function(rf){
  
  norm_importance <- 1- (rf$importance / sum(rf$importance))
  names(norm_importance) <- row.names(rf$importance)
  sortedImpMatrix = sort(norm_importance, decreasing=F)
  
  sortedImpMatrix
  
}

########################################################
########################################################

### get model size based on FSR
### @ fsr = vector of false selection rates
### @ alpha = p-to-enter in model

getModelSize <- function(fsr, alpha){
  model_size = which(fsr > alpha)[1] - 1
  return(model_size)
}


########################################################
########################################################

### Create phony variables
###
### input 
### @ design.matrix = input 
### 
### output 
### @ design.matrix = input 


generatePermutedMatrix <- function(design.matrix){
  
  N = nrow(design.matrix)
  p = ncol(design.matrix)
  permutedMatrix <- design.matrix[sample(N, size=N, replace=F),]
  names(permutedMatrix) <- paste('phony', c(1:p))
  design.matrix <- cbind(design.matrix, permutedMatrix)
  
  design.matrix
  
}