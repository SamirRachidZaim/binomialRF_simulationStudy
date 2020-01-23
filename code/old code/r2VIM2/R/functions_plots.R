
#' Manhattan plot.
#'
#' This function generates a Manhattan plot for minimal relative importance scores of random forest variable selection or log transformed p-values of a statistical test.
#' 
#' @param res data.frame with results or NULL if results should be loaded from results.file
#' @param results.file a character string with the name of the results file or NULL if parameter res is used
#' @param rel.imp logical. If 'TRUE' function assumes that results from random forest variable selection should be plotted (columns 'varname', 'chr', 'pos', 'rel_imp_min'). If 'FALSE' function assumes that results of statistical test should be plotted (columns 'SNP', 'CHR', 'BP', 'P').
#' @param col vector of colors (if the length is smaller than the number of chromosomes colors will be recycled)
#' @param var.highlight vector with SNPs that will be plotted in red 
#' @param col.highlight color for SNPs given in var.highlight (default: black) 
#' @param threshold vector giving values at which a horizontal line should be drawn (Note: given as untransformed threshold for statistical test)
#' 
#'  @export

manhattanplot <- function(res = NULL, results.file = NULL, rel.imp = TRUE, col = NULL,
                          var.highlight = NULL, col.highlight = "black", threshold = NULL) {

    ## load results
    if (is.null(res)) {
        res = read.table(file = results.file, header = TRUE, as.is = TRUE)
    }
    
    ## define columns with variable names and y values
    if (rel.imp) {
        col.varname = "varname"
        col.chr = "chr"
        col.pos = "pos"
        col.y = "rel_imp_min"
        ylab = "Minimal relative importance"
    } else {
        col.varname = "SNP"
        col.chr = "CHR"
        col.pos = "BP"
        col.y = "P"
        ylab = "-log10(P-value)"
        threshold = -log10(threshold)

        ## take logarithm
        res[, col.y] = -log10(res[, col.y])
    }

    ## use only autosomal chromosomes
    if (!is.numeric(res[, col.chr])) {
        stop("Chromosome has to be numeric!")
    }
    res = res[res[, col.chr] %in% 1:22,]
    
    ## sort res by chromosome and position
    res = res[order(res[, col.chr], res[, col.pos]),]

    ## check if res contains results for all or only a single chromosome
    chr = unique(res[, col.chr])
    no.chr = length(chr)

    ## define color
    if (is.null(col)) {
        col = rainbow(length(chr))
    } 
    if (length(col) < length(chr)) {
        col = rep(col, length(chr))[1:length(chr)]
    }

    ## positions on x-axis for each SNP
    if (no.chr == 1) {
        ## use bp information
        res$pos.plot = res[, col.pos]
        ticks = seq(min(res$pos.plot), max(res$pos.plot), 50 * 10^6)
        labs = ticks / 10^6
        if (length(ticks) == 1) {
            ticks = seq(min(res$pos.plot), max(res$pos.plot), 5*10^5)
            labs = round(ticks / 10^6, digits = 1)
        }
        xlab = paste("Chromosome", chr, "(Mb)")
        col.plot = col
    } else {
        ## set position relative to first SNP on each chromosome
        pos.plot = rep(NA, nrow(res))
        start = rep(NA, length(chr) + 1)
        start[1] = 0
        for (i in chr) {
            ind = which(res[, col.chr] == i)
            min = res[1, col.pos]
            pos.plot[ind] = res[ind, col.pos] - min + 1 + start[i]
            start[i+1] = max(pos.plot[ind])
        }
        res$pos.plot = pos.plot
        
        ## get positions for chromosome labels
        ticks = start[-length(start)] + (start[-1] - start[-length(start)])/2
        xlab = "Chromosome"
        labs = chr

        col.plot = col[res[, col.chr]]
    }


    ## plot
    plot(res[, c("pos.plot", col.y)],
         main = "Manhattan plot", xlab = xlab, ylab = ylab,
         pch = 19, las = 1, cex = 0.5, xaxt = "n", col = col.plot)
    axis(side = 1, at = ticks, labels = labs)
    if (!is.null(var.highlight)) {
        ind = which(res[, col.varname] %in% var.highlight)
        points(res[ind, c("pos.plot", col.y)], pch = 19, col = col.highlight)
    }
    if (!is.null(threshold)) {
        abline(h = threshold, lty = 2)
    }
    
}


#' Plot to compare results of random forest variable selection and statistical test.
#'
#' This function generates a scatterplot of minimal relative importance of random forest variable selection and log transformed p-values of a statistical test.
#'
#' @param res.rf data.frame with random forest results (columns 'varname', 'chr', 'pos', 'rel_imp_min') or NULL if results should be loaded from rel.imp.file
#' @param res.test data.frame with test results (columns 'SNP', 'CHR', 'BP', 'P') or NULL if results should be loaded from pvalue.file
#' @param rel.imp.file a character string with name of file with results from random forest variable selection (columns 'varname', 'chr', 'pos', 'rel_imp_min') or NULL if parameter res.rf is used.
#' @param pvalue.file a character string with name of file with PLINK results from test (columns 'SNP', 'CHR', 'BP', 'P') or NULL if parameter res.test is used.
#' @param thresholds vector with thresholds for relative importance and p-value (defaults: 2, 5*10^-8)
#' @param var.highlight vector with SNPs that will be plotted in red 
#' @param main title for plot (default none)
#' 
#'  @export

scatterplot.rel.imp.pvalue <- function(res.rf = NULL, res.test = NULL,
  rel.imp.file = NULL, pvalue.file = NULL, 
  thresholds = c(2, 5*10^-8), var.highlight = NULL, main = "") {

  if (is.null(res.rf)) {
    res.rf = read.table(file = rel.imp.file, header = TRUE, as.is = TRUE)
  }
  if (is.null(res.test)) {
    res.test = read.table(file = pvalue.file, header = TRUE, as.is = TRUE)
  }
  if (nrow(res.rf) != nrow(res.test)) {
    stop("Different number of variables!")
  }

  ## sort
  res.rf = res.rf[order(res.rf$varname),]
  res.test = res.test[order(res.test$SNP),]

  if (nrow(res.rf) > 100000) {
    warning("only most interesting SNPs are plotted!")

    ## select SNPs selected with either method
    snps.lr = res.test[res.test$P < 0.05, "SNP"]
    snps.rf = res.rf[res.rf$rel_imp_min > 0, "varname"]
    snps = union(snps.lr, snps.rf)
    res.rf = res.rf[res.rf$varname %in% snps,]
    res.test = res.test[res.test$SNP %in% snps, ]
  }

  if (!all.equal(res.rf$varname, res.test$SNP)) {
    stop("Different order of variables!")
  }

  ## combine
  data = data.frame(res.rf[, c("varname", "rel_imp_min")], log.p = -log10(res.test$P))

  ## plot
  plot(data[, 2:3], pch = 19,
       xlab = "Minimal relative importance", ylab = "-log10(P-value)",
       main = main, las = 1)
  if (!is.null(thresholds)) {
      abline(v = thresholds[1], lty = 2)
      abline(h = -log10(thresholds[2]), lty = 2)
  }
  if (!is.null(var.highlight)) {
    points(data[data$varname %in% var.highlight, 2:3],
           pch = 19, col = "red")
  }
 
}


