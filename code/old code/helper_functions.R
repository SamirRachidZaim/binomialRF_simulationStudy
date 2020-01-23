##################################################################
##################################################################
######
###### Samir Rachid Zaim Final Project
######
######
###### Helper Functions
######
##################################################################
##################################################################



##################################################################
##################################################################
##### find top K important variables in Random Forest

findTopGini <- function(rf.object){
  full_rf_imp <- data.frame(importance(rf.object))
  full_rf_imp$Var <- row.names(full_rf_imp)
  full_rf_imp<- full_rf_imp[order(full_rf_imp$MeanDecreaseGini, decreasing = T),]
  return(Importance =full_rf_imp)
}

##### rank_important_variables
#####
#####   desc: builds and plot rankings based on split choices
#####
#####
#####   input:#####     - RF                = random forest object | list of decision trees (split_matrix)
#####
#####   output:
#####     - barplot           = barplot comparison of all features rankings in first 3 node splits
#####                           vs. top M variables
#####     - top_ranked        = character vector of the top M ranked variables
#####                         where M is the minimum of the # of desired (input)
#####                         and p*, where p* < p are the variables that were chosen at least once

#####   future work
#####     - smarter/less naive ranking
#####     - maybe weight rankings by node appearance (i.e. top node selection = 3x 3rd node selection)
#####

###### Find function to identify top variables
rank_important_variables <- function(RF, top_ranked_m = 20, weighted_ranking =F, name = NULL){

  if(weighted_ranking){
    ranked_vars = get_weighted_variable_rank(RF)
    name = paste('Weighted',name)
  } else{
    ranked_vars = get_variable_rank(RF)
  }

  #### Weighted or Unweighted Ranking
  #### Produced 'RankedVars' List
  #### Used below to plot and assess

  pdf(paste('~/Desktop/classes/fall2017/MATH574/FinalProject/',name,'.pdf', sep=''))

  par(mfrow=c(2,1))
  barplot(ranked_vars, ylim = c(0,max(ranked_vars)+20),
          xlab = "Variable", ylab='Frequency',
          main = paste(name, 'Full set of Ranked Variables'))

  #### Top ranked M varialbes
  #### m = minimum (desired var length, total variable length)
  ####

  m = min(top_ranked_m, length(ranked_vars))
  top_ranked      = ranked_vars[1:m]
  barplot(top_ranked, ylim = c(0,max(ranked_vars)+20),
          xlab = "Variable", ylab='Frequency',
          main = paste(name, 'Top Ranked Variables'))
  dev.off()

  return(names(top_ranked))

}

##################################################################
##################################################################

##### Get weighted ranking 
get_weighted_variable_rank <- function(RF){
    RF_split_matrix = do.call(rbind, RF)
    RF_split_matrix$Split_Var = as.character(RF_split_matrix$Split_Var)
    root_node_idx         = grep('1|1.', row.names(RF_split_matrix))
    second_level_node_idx = grep('[2-3]|[2-3].', row.names(RF_split_matrix))
    third_level_node_idx = grep('[4-7]|[4-7].', row.names(RF_split_matrix))
    
    weighted_var_vec = c(rep(RF_split_matrix[root_node_idx,'Split_Var' ], 3),
                         rep(RF_split_matrix[second_level_node_idx,'Split_Var' ], 2),
                         rep(RF_split_matrix[third_level_node_idx,'Split_Var' ], 1))
    
    weighted_var_vec = weighted_var_vec[!weighted_var_vec %in% 'Leaf Node']
    
    ranked_vars <- sort(table(weighted_var_vec), decreasing = T)
}

get_variable_rank <- function(RF){
  
  RF_split_matrix = do.call(rbind, RF)$Split_Var
  RF_split_matrix = RF_split_matrix[!RF_split_matrix %in% 'Leaf Node']
  RF_split_matrix = droplevels(RF_split_matrix)
  ranked_vars     = sort(table(RF_split_matrix), decreasing = T)
  
}

##################################################################
##################################################################
##### Consider Expanding to weighted cost function

##### Classification Losses

ClassificationLossFunction <- function(preds, y, method='0-1'){
    if(method =='0-1'){
        loss =sum(preds!=y)/length(y)
    }
    return(loss)
}

##### Regression Losses
RegressionLossFunction <- function(preds, y, method='LeastSquares'){

    if(method =='LeastSquares'){
        loss = sum((y-preds)^2)
        names(loss) <- method

    } else if(method =='LAD'){
        loss = sum(abs(y-preds))
        names(loss) <- method

    }

    return(loss)

}

######
##################################################################
##################################################################

##################################################################
##################################################################
##### evaluate_var
#####
#####   desc: calculates quality of splitting variable
#####
#####   input:
#####     - chosen_var = 'ith splitting var' | character
#####     - grid       = 'ith splitting var possible values'| numeric
#####     - X          = 'subsampled (in n, p) design matrix'
#####     - y          = subsampled target
#####
#####   output:
#####     - Best Error
#####     - Best Split Value

evaluate_var.regression <- function(chosen_var, grid, X, y, loss='LeastSquares'){

  if( sum(y)==0 | sum(y) == length(y)   ){

    ## If y's are all 0s or 1 nothing to do.
    return(data.frame(Error = 0, Split_Value = 0, Split_Var = 'Leaf Node'))

  }else{
    ## Evaluate Loss Function for all possible levels
    evals = unlist(mclapply(grid, function(x) RegressionLossFunction(preds = X[, chosen_var]<x, y=y, method=loss)))

    ## Find best splitting point (argmin)
    best.split_point <- grid [which.min(evals)]
    best.split_error <- evals[which.min(evals)]

    ## Return argmin
    return(data.frame(Error = best.split_error, Split_Value = round(best.split_point,3), Split_Var = as.character(chosen_var)))
  }

}

evaluate_var.classification <- function(chosen_var, grid, X, y, loss='0-1'){

  if( sum(y)==0 | sum(y) == length(y)   ){

    ## If y's are all 0s or 1 nothing to do.
    return(data.frame(Error = 0, Split_Value = 0, Split_Var = 'Leaf Node'))

  }else{
    ## Evaluate Loss Function for all possible levels
    evals = unlist(mclapply(grid, function(x) ClassificationLossFunction(preds =X[, chosen_var]<x, y=y, method=loss)))

    ## Find best splitting point (argmin)
    best.split_point <- grid [which.min(evals)]
    best.split_error <- evals[which.min(evals)]

    ## Return argmin
    return(data.frame(Error = best.split_error, Split_Value = round(best.split_point,3), Split_Var = as.character(chosen_var)))
  }

}



##### choose_best_split
#####
#####   desc: calculates loss func. for all splitting var
#####         and chooses best splitting variable
#####
#####   input:
#####     - X          = 'subsampled (in n, p) design matrix'
#####     - y          = subsampled target
#####
#####   output:
#####     - Best Error
#####     - Best Split Value

choose_best_split.regression <- function(X, y, percent_features = .5, loss='LeastSquares'){

  ### Make Sure Valid Index is Given

  if(length(y)==0){

    return(data.frame(Error = numeric(1), Split_Value = numeric(1), Split_Var = 'Leaf Node'))

  } else{

    ### For each split pick 60% of the samples
    ### and 40% of the features
    samp_features <- sample(1:dim(X)[2], size = floor(dim(X)[2]*percent_features))

    X_new = X[, samp_features]
    y_new = y

    err.mat = data.frame(Error=numeric(0),
                         Split_Value=numeric(0),
                         Split_Var = character(0),
                         stringsAsFactors = F)

    ### This can be parallelized
    err.mat = mclapply(1:length(names(X_new)), function(i) rbind(err.mat,evaluate_var.regression(names(X_new)[i], unique(X_new[, chosen_var]), X_new, y_new, loss=loss)))
    err.mat$NodeSize <- length(y_new
                               )
    best_combo = which.min(err.mat$Error)
    return((err.mat[best_combo,]))
  }

}


### classification
choose_best_split.classification <- function(X, y, percent_features = .5, loss='0-1', minNode = 5){

    ### Only split Nodes with at least 5 observations
    ### that are not pure nodes

  if(length(y) < minNode | sum(y)==0 | sum(y) == length(y)){

    return(data.frame(Error = numeric(1), Split_Value = numeric(1), Split_Var = 'Leaf Node', NodeSize=length(y)))

  } else{

    ### For each split pick 60% of the samples
    ### and 40% of the features
    samp_features <- sample(1:dim(X)[2], size = floor(dim(X)[2]*percent_features))

    X_new = X[, samp_features]
    y_new = y

    err.mat = data.frame(Error=numeric(0),
                         Split_Value=numeric(0),
                         Split_Var = character(0),
                         stringsAsFactors = F)

    ### This can be parallelized
      err.mat = mclapply(1:length(names(X_new)), function(i)
          evaluate_var.classification(names(X_new)[i], unique(X_new[, names(X_new)[i]]), X_new, y_new, loss=loss))

      err.mat = do.call(rbind, err.mat)
      err.mat$NodeSize <- length(y_new)

    best_combo = which.min(err.mat$Error)
    return((err.mat[best_combo,]))
  }

}


##### findAncestors
#####
#####   desc: recursively finds all ancestral
#####         nodes from current node to the root
#####

findAncestors <- function(nodeNumber){

    while(nodeNumber[length(nodeNumber)] != 1){
        nodeNumber = c(nodeNumber, floor(nodeNumber[length(nodeNumber)] / 2))
    }
    return(nodeNumber)
}

##### findNodeIndex
#####
#####   desc: recursively finds all set of observations
#####         for a current node to choose its next splits
#####

is.even <- function(x) x %% 2 == 0

findEvenOddIndex <- function(split_matrix, X, ancestor){

    if(is.even(ancestor)){
        return(which(X[, split_matrix$Split_Var[ancestor]] < split_matrix$Split_Value[ancestor]))
    } else{
        return(which(X[, split_matrix$Split_Var[ancestor]] > split_matrix$Split_Value[ancestor]))
    }

}


findNodeIndex <- function(split_matrix, X, ancestors){

    index_list= list()

    for(i in 1:length(ancestors)){
        x = ancestors[i]

        if(x ==1){

        } else {

            if(is.even(x)){
               index_list[[i]] <- which(X[, split_matrix$Split_Var[ancestors[i+1]]] < split_matrix$Split_Value[ancestors[i+1]])
            } else {
                index_list[[i]] <- which(X[, split_matrix$Split_Var[ancestors[i+1]]] > split_matrix$Split_Value[ancestors[i+1]])
            }


        }
    }

    return(list(Index_list = index_list, NodeIndex =Reduce(intersect, index_list)))
}



##### build_classification_tree
#####
#####   desc: builds decisions tree of depth 3
#####         based on choose_best_split results
#####
#####   input:
#####     - node_depth        = how deep should tree be built (NA for moment)
#####     - X                 = 'subsampled (in n, p) design matrix'
#####     - y                 = subsampled target
#####     - percent_data      = how much data should be used in building tree (sampled at each node)
#####     - percent_features  = how many features should be sampled at each node
#####
#####   output:
#####     - split_matrix      = matrix specifying variable and varaible value split for each node
#####                         where each node is labelled by counting left to right and node 1
#####                         is the root node.

#####   future work:          = at the moment first 3 nodes are hard coded;
#####                         need to generalize later to m-nodes

build_classification_tree <- function(node_depth=2, X, bool.y, percent_data=.6, percent_features=.5, minNode=5){

    y = bool.y
  
    require(parallel)

    if(percent_data < 0 | percent_data > 1 | percent_data < 0 | percent_features > 1){
        stop("Error: percent_data or percent_features are outside their acceptable (0-1) range")
    }

    if(class(y) != 'logical'){
        stop('Error: Make Y vector a logical vector for classification')

    }else if(class(y) == 'logical'){
      
        #### Data for decision tree should be sampled only once
        #### but features should be subsampled at each node
        samp_rows     <- sample(1:dim(X)[1], size = floor(dim(X)[1]*percent_data))
        X = X[samp_rows,]
        y = y[samp_rows]


        ## Root node split
        root_var_split = choose_best_split.classification(X,y,percent_features=percent_features, minNode=5 )
        root_var_name  = as.character(root_var_split$Split_Var)
        root_var_value = root_var_split$Split_Value

        if(root_var_name=='Leaf Node'){
            return(root_var_split)

        } else{

            root_idx = which(X[, root_var_name] < root_var_value)


            split_matrix = root_var_split
            ## To calculate left/ride node split
            ## Find previous index
            ## Calculate left/right node split

            for(i in 2:node_depth){ #1st split is hard-coded

                startNodeNum = 2^(i-1)
                endNodeNum = 2^i -1

                for(currNode in startNodeNum:endNodeNum){

                    ancestors = findAncestors(currNode)
                    index = findNodeIndex(split_matrix, X, ancestors)$NodeIndex

                    curr_split_matrix <- choose_best_split.classification(X[index,],y[index],percent_features=percent_features, minNode=5 )
                    split_matrix <- rbind(split_matrix, curr_split_matrix)
                }



            }
            row.names(split_matrix) <- paste("node",c(1:(2^node_depth-1)))

            return(list(Tree= split_matrix, TreeError = calculateTreeError(split_matrix) ))
        }

    }
 }

findKInteractions.tree <- function(split_matrix, k){

    require(parallel)

    ## Look at K-way interactions, starting
    ## at root node
  
    if('Split_Var' %in% names(split_matrix)){
      return(split_matrix['Split_Var'])
    } else if('Tree' %in% names(split_matrix)){
      split_matrix <- split_matrix$Tree
      
      startNodeNum = 2^(k-1)
      endNodeNum = 2^k -1
      idx = startNodeNum:endNodeNum
      
      interaction_list = mclapply(idx, function(x) as.character(split_matrix$Split_Var[findAncestors(x)]))
      removeLeafNodes <- which(sapply(1:length(idx), function(x) !"Leaf Node" %in% interaction_list[[x]])==T)
      interaction_list <- interaction_list[removeLeafNodes]
      return(interaction_list)
      
    }
}

calculateBinomialP <- function(L, percent_features){
  m = floor(L * percent_features)
  prod.vector = sapply(1:m, function(x) (L-x)/(L-(x-1)))
  (1-prod(prod.vector))*(1/m)
}

findTreeWeight <- function(Tree, k){
  mat = data.frame(Var = findKInteractions.tree(Tree, k), Weight=Tree$TreeError['Weight'])
  names(mat)= c('Var','Weight')
  mat
}

findKInteractions.rf <- function(RF, k=1, plot=F, ntrees=20, percent_features, numVars){
  require(data.table)
   
    interaction_list <- mclapply(RF, function(x) as.character(x$Tree$Split_Var))
    interaction_list <- data.frame(do.call(rbind, interaction_list), stringsAsFactors = F)
    interaction_list <- data.table(t(do.call(cbind, lapply(1:ncol(interaction_list) , function(x)  t(interaction_list[,x])))))
    interaction_list <- data.frame(table(interaction_list))
    interaction_list <- interaction_list[order(-interaction_list$Freq),]
    
    p = calculateBinomialP(numVars, percent_features)
    
    interaction_list$Pvalue <- as.numeric(sapply(interaction_list$Freq, function(x) binom.test(x, n= 20, p, alternative='greater')$p.value))
    interaction_list$AdjPvalue <- p.adjust(interaction_list$Pvalue, method = fdr.method)
    
    colnames(interaction_list) <- c("Var",'Freq', 'Pval', "adjPval")
    # interaction_list = mclapply(1:length(RF), function(x) findTreeWeight(RF[[x]], k))
    # interaction_list = data.table(do.call(rbind, interaction_list))
    # interaction_list = interaction_list[, list(Weight=sum(Weight), UnweightedRank=.N), by=Var] 
    # interaction_list = interaction_list[order(-Weight)]
    # 
    # interaction_list$Weight <- interaction_list$Weight / sum(interaction_list$Weight)
    # interaction_list$UnweightedRank <- interaction_list$UnweightedRank / sum(interaction_list$UnweightedRank)
    
    if(plot){
        require(ggplot2)
        barplot(interaction_list$Weight, arg.names=interaction_list$Var, main='Interactions', xlab= 'Interaction Type', ylab='Frequency')
        ggplot(interaction_list, aes(x=Var,y=Weight)) +  geom_bar(stat='identity')
    }

    return(interaction_list)
}

findImportantKInteractionVars <- function(X, y, ntrees= 30, k){
    idx = sample(length(X[,1]), .8 * length(X[,1]))
    X = X[idx,]
    y = y[idx]

    RF = mclapply(1:ntrees, function(x) build_classification_tree(node_depth=3, X, y))
    impVars    = findKInteractions.rf(RF, k =k, plot=F)
    return(impVars)
}

####### 
### U(i)/(S(x))

calculateFSR <- function(modelVars,trueVars){
  
  Uninformative <- sum(!modelVars %in% trueVars)
  selected <- length(modelVars) +1
  fsr <- Uninformative/selected
  return(fsr)
  
}

calculateDiscoveryRate <- function(modelVars,trueVars){
  
  selected <- sum(modelVars %in% trueVars)
  NTrueVars <- length(trueVars) 
  discoveryRate <- selected/NTrueVars
  return(discoveryRate)
  
}

calculateTreeError <- function(split_matrix){
  idx = which(split_matrix$Split_Var !='Leaf Node')
  treeErr = weighted.mean(split_matrix$Error[idx], split_matrix$NodeSize[idx])
  treeErr = c(treeErr, round(1/(treeErr+1),2))
  names(treeErr) <- c('Error','Weight')
  treeErr
}

RF_Soil <- function(X, y , Xtest,ytest, ntrees=20){
  require(SOIL)
  require(parallel)
  require(data.table)
  require(randomForest)
  
  b_BIC <- SOIL(X, y, family = "binomial", weight_type = "BIC")
  
  M = length(b_BIC$candidate_models_cleaned[,1])
  
  ### Generate SOIL Var List 
  varlist <-  mclapply(1:M,  function(i) b_BIC$candidate_models_cleaned[i,] * c(1:length(b_BIC$candidate_models_cleaned[1,])))
  
  ### Generate RF Equivalent Model
  RFSoil <- mclapply(1:M, function(i) 
    if(sum(varlist[[i]]!= 0)>1){
        randomForest(as.matrix(X[, varlist[[i]]]), factor(y), ntree=ntrees)
      }else {
      return(NULL)
    }
  )
  
  SOIL.gini = mclapply(1:M, function(x) 
    if(sum(varlist[[x]]!= 0)>1){
      importance(RFSoil[[x]])
    } 
  )

  
  ## Correct for Matrices 
  
  GiniMatrix <- data.frame(VariableName= character(0),
                           ErrorWeight = numeric(0),
                           MeanDecreaseGini=numeric(0))
  
  trainErr <- testErr <-  numeric(0)
  for(i in 1:M){
    
    if(sum(varlist[[i]]!= 0)>1){
      new.mat = data.frame(SOIL.gini[[i]])
      new.mat$VariableName =row.names(new.mat)
      new.mat$ErrorWeight  = 1/last(RFSoil[[i]]$err.rate[,1])
      
      GiniMatrix = rbind(GiniMatrix, new.mat)
      trainErr <-c(trainErr, round(mean(RFSoil[[i]]$predicted != y, na.rm=T), 2))
      testErr  <-c(testErr , round(mean(predict(RFSoil[[i]] , Xtest[, varlist[[i]]]) != ytest, na.rm=T),2))
    }
  }
  
  GiniMatrix <- data.table(GiniMatrix)  
  
  GiniMatrix <- GiniMatrix[, list(MeanDecreaseGini_Average = sum(MeanDecreaseGini), 
                                  MeanDecreaseGini_WeightedAverage = sum(MeanDecreaseGini * ErrorWeight)), by=VariableName]
  
  GiniMatrix$MeanDecreaseGini_Average <- GiniMatrix$MeanDecreaseGini_Average / sum(GiniMatrix$MeanDecreaseGini_Average)
  GiniMatrix$MeanDecreaseGini_WeightedAverage <- GiniMatrix$MeanDecreaseGini_WeightedAverage / sum(GiniMatrix$MeanDecreaseGini_WeightedAverage)
  
  return(list(GiniMatrix =GiniMatrix, 
              ModelList  = RFSoil ,
              TrainErr   = trainErr,
              TestErr    = testErr))
}



