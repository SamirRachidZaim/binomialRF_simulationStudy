
#' Generate several random forests with different random number seeds using Random Jungle.
#'
#' This function generates several random forests for the same data set using different
#' random number seeds. It uses the software Random Jungle
#' \url{http://imbs-luebeck.de/imbs/de/node/227}.
#'
#' Note:
#' This function assumes that the programs rjungle and sed are in the path.
#'
#' @param data.file a character string with the name of the input file. This file might be
#' generated from binary PLINK files using \code{\link{prepare.input.file.rj}}.
#' @param out.prefix a character string with the prefix of the output files. For description
#' of the different files see manual for Random Jungle.
#' @param no.runs number of forests to be generated
#' @param type type of random forest analysis (must be either 'classification' or 'regression')
#' @param ntree number of trees in each forest
#' @param mtry number of variables selected for each split
#' @param nodesize minimum size of terminal nodes. Note that it should be set to about 10\% of the sample size.
#' @param var.pheno name of column containing phenotype information (default value 'PHENOTYPE' used by \code{\link{prepare.input.file.rj}}).
#' @param seed random number seed for first forest (additional forests will use <seed> + 1, <seed + 2>, ...)
#' @param no.threads number of threads for generating random forests in parallel mode
#' 
#' @export

random.jungle.runs <- function(data.file, out.prefix, no.runs = 1, type = "classification",
                               ntree, mtry, nodesize, var.pheno = "PHENOTYPE",
                               seed = 123456, no.threads = 1) {
  
  ## check OS
  windows = ifelse(.Platform$OS.type == "windows", TRUE, FALSE)
  
    ## check type
    if (!(type %in% c("classification", "regression"))) {
        stop("type must be either 'classification' or 'regression'!")
    }
    
    ## check if data file exists
    if (!file.exists(data.file)) {
        stop(paste("file", data.file, "does not exist!"))
    }
    
    ## check if data file contains missing values
    if (length(grep("gz$", data.file)) > 0) {
        if (windows) { 
          cmd1 = paste("gzip -d -c", data.file, ">", paste(out.prefix, "data", sep="."))
          shell(cmd1)
          cmd = paste("sed '/NA/!d;q'", paste(out.prefix, "data", sep="."))
        } else { 
          cmd = paste("gunzip -c", data.file, "| sed '/NA/!d;q'")
        }   
    } else {
        cmd = paste("sed '/NA/!d;q'", data.file)
    }
    
    no.miss = length(try(system(cmd, intern = TRUE)))
    if (no.miss > 0) {
        stop(paste(data.file, "contains missing values!"))
    }
    
    for (r in 1:no.runs) {
        print(paste("run", r))
        random.jungle(data.file = data.file,
                      out.prefix = paste(out.prefix, "_run_", r, sep = ""),
                      type = type,
                      ntree = ntree,
                      mtry = mtry,
                      nodesize = nodesize,
                      var.pheno = var.pheno,
                      seed = seed + r - 1,
                      no.threads = no.threads)
    }
}



#' Generate random jungle.
#'
#' This function generates a random forest using the software Random Jungle
#' \url{http://imbs-luebeck.de/imbs/de/node/227}.
#'
#' Note:
#' This function assumes that the program rjungle is in the path.
#' 
#' @param data.file a character string with the name of the input file. This file might be
#' generated from binary PLINK files using \code{\link{prepare.input.file.rj}}.
#' @param out.prefix a character string with the prefix of the output files. For description
#' of the different files see manual for Random Jungle.
#' @param type type of random forest analysis (must be either 'classification' or 'regression')
#' @param ntree number of trees in the forest
#' @param mtry number of variables selected for each split
#' @param nodesize minimum size of terminal nodes. Note that it should be set to about 10\% of the sample size.
#' @param var.pheno name of column containing phenotype information (default value 'PHENOTYPE' used by \code{\link{prepare.input.file.rj}}).
#' @param seed random number seed
#' @param no.threads number of threads for running Random Jungle in parallel mode
#' 
#' @export

random.jungle <- function(data.file, out.prefix, type = "classification", ntree, mtry, nodesize, var.pheno = "PHENOTYPE", seed = 123456, no.threads = 1) {

  if (!file.exists(data.file)) {
    stop(paste("file", data.file, "does not exist!"))
  }
  
  ## set parameters
  if (type == "classification") {
    tree.type = 1
  } else if (type == "regression") {
    tree.type = 3
  } else {
    stop(paste("type", type, "not known!"))
  }
  
  cmd = paste("rjungle",
    "--file", data.file,
    "--treetype", tree.type,
    "--ntree", formatC(ntree, format = "d"),
    "--mtry", formatC(mtry, format = "d"),
    "--targetpartitionsize", formatC(nodesize, format = "d"), ## terminal node size
#    "--maxtreedepth", depth,
    "--impmeasure", 4, ## permutation (unscaled)
    "--depvarname",  var.pheno,
    "--seeed", seed,
    "--nthreads", no.threads,
    "--outprefix", out.prefix)
  try(system(cmd))
}

#' Combine results of random jungle runs.
#'
#' This function combines variable importance scores and errors of several random jungle runs. It calculates relative importance scores for each variable in each run as well as minimal and median relative importance.
#' 
#' @param res.prefix a character string with the prefix of the output files of Random Jungle runs ('run_<x>' will be added automatically)
#' @param no.runs number of forests that should be combined (function will combined forests 1, 2, ..., <no.runs>
#' @param map.file a character string with name of map file (either MAP or BIM format)
#' @param out.prefix a character string with the prefix of the output files. The function generates two files: <out.prefix>_all_runs.info (errors and minimal importance scores) and <out.prefix>_all_runs_relative.importance
#' 
#'  @export

combine.results.random.jungle.runs <- function(res.prefix, no.runs, map.file = NULL,
                                               out.prefix) {

  ## check OS
  windows = ifelse(.Platform$OS.type == "windows", TRUE, FALSE)
  
  ## extract variable importance for each run
    for (r in 1:no.runs) {
        cols = "5"
        if (r == 1) {
            cols = paste("3,", cols, sep = "")
        }
        file.imp.run = paste(res.prefix, "_run_", r, ".importance2", sep = "")
        cmd = paste("sed 1d", file.imp.run, 
            "| sort -n -k3",
            "| cut -d' '", paste("-f", cols, sep = ""),
            ">", paste(file.imp.run, ".sorted", sep = "")) 
        if (windows) {
          shell(cmd)
        } else {
          try(system(cmd))
          }
        
    }

    ## combine in one file
    cmd = paste('paste -d" "',
        paste(paste(res.prefix, '_run_', 1:no.runs, '.importance2.sorted', sep = ''),
              collapse = ' '),
        '>', paste(res.prefix, '_all_runs.importance', sep = ''))
    if (windows) {
      shell(cmd)
    } else {
      try(system(cmd))
    }
    
    ## add header line
    cmd = paste("sed -i '1i ",
        paste("varname", paste(paste("raw_score_run_", 1:no.runs, sep = ""), collapse = " ")),
        "' ", paste(res.prefix, "_all_runs.importance", sep = ""), sep = "")
    try(system(cmd))

    ## load importance info
    var.imp = read.table(file = paste(res.prefix, "_all_runs.importance", sep = ""),
        header = TRUE, as.is = TRUE, row.names = 1)

    ## minimal importance for each run
    min.imp = apply(var.imp, 2, min)

    if (any(min.imp == 0)) {
        stop("Minimal importance score is not negative for at least one run!")
    }
    
    ## relative importance
    var.imp.rel = .calculate.rel.var.imp(var.imp = var.imp)
    
    info.imp = t(apply(var.imp.rel, 1, function(x) {
        return(c(rel_imp_min = min(x), rel_imp_med = median(x)))}))
    
    ## add positional info (if map file is given)
    if (!is.null(map.file)) {
        info.var = read.table(file = map.file, header = FALSE, as.is = TRUE)[, c(2, 1, 4)]
        dimnames(info.var) = list(info.var[, 1], c("varname", "chr", "pos"))
        info.var = info.var[rownames(var.imp),]
    } else {
        info.var = data.frame(varname = rownames(var.imp), stringsAsFactors = FALSE)
    }

    ## check order of variables
    if (!all.equal(info.var$varname, rownames(var.imp.rel)) |
        !all.equal(info.var$varname, rownames(info.imp))) {
        stop("Different order of variable names!")
    }
  
    ## write file
    var.imp.new = data.frame(info.var, var.imp.rel, info.imp)
    if ("chr" %in% colnames(var.imp.new)) {
        var.imp.new = var.imp.new[order(var.imp.new$chr, var.imp.new$pos),]
    }
    write.table(var.imp.new, file = paste(out.prefix, "_all_runs_relative.importance", sep = ""),
                col.names = TRUE, row.names = FALSE, quote = FALSE)
    try(system(paste("rm", paste(res.prefix, "_all_runs.importance", sep = ""))))

    ## summarize errors
    info.error = NULL
    for (r in 1:no.runs) {
        res.prefix.run = paste(res.prefix, "_run_", r, sep = "")
        
        if (file.exists(paste(res.prefix.run, "confusion2", sep = "."))) {
            ## extract sensitivity, specificity and accuracy
            temp = read.table(file = paste(res.prefix.run, "confusion2", sep = "."),
                header = TRUE, as.is = TRUE)
            res = c(1-temp[4], 1-temp[3], temp$error)
            names(res) = c("sens", "spec", "error")
            info.error = rbind(info.error, res)
        } else {
            ## extract mse and rsq
            error = readLines(paste(res.prefix.run, "confusion", sep = "."))
            line.rsq = grep("^Accuracy", error, value = TRUE)
            line.rsq = unlist(strsplit(line.rsq, " "))
            rsq = as.numeric(line.rsq[length(line.rsq)])
            
            line.mse = grep("(SS_Residual)", error, value = TRUE)
            line.mse = unlist(strsplit(line.mse, " "))
            mse = as.numeric(line.mse[length(line.mse)-1])
            info.error = rbind(info.error,
                c(rsq = rsq, mse = mse))
        }
    }
    info = cbind(run = 1:no.runs, info.error, min.imp)
    write.table(info, file = paste(out.prefix, "_all_runs.info", sep = ""),
              col.names = TRUE, row.names = FALSE, quote = FALSE)
}


#' Calculate relative variable importance.
#'
#' This internal function calculates relative variable importance.
#' 
#' @param var.imp matrix with variable importance for each variable (rows) in each run (columns)
#' 
#' @return matrix with relative variable importance for each variable (rows) in each run (columns)

.calculate.rel.var.imp <- function(var.imp) {
  ## check if global min is negative
  global.min = min(var.imp)
  if (global.min >= 0) {
    stop("No variable with negative importance score in all runs!")
  }

  ## relative variable importance
  no.neg = 0
  rel.var.imp = apply(var.imp, 2, function(x) {
    min = min(x)
    if (min >= 0) {
      no.neg = no.neg + 1
      min = global.min 
    }
    return(x / abs(min))})
  if (no.neg > 0) {
    warning(paste("No variable with negative importance score in", no.neg, "runs (using global minimum)!"))
  }              
  colnames(rel.var.imp) = paste("rel_imp_run_", 1:ncol(rel.var.imp), sep = "")
  return(rel.var.imp)
}
