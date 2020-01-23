
#' Prepare input files for Random Jungle software.
#'
#' This function converts binary PLINK files to input files for Random Jungle. It uses PLINK
#' (\url{http://pngu.mgh.harvard.edu/~purcell/plink/}) to
#' change genotypes to 0, 1, 2 coding giving the number of minor alleles. This file is
#' further modified (e.g. the minor alleles are cut from the SNP names and the columns
#' with pedigree and gender information are removed).
#'
#' Note:
#' This function assumes that the following programs are in the path:
#' \itemize{
#'   \item PLINK
#'   \item sed
#'   \item cut
#'   \item awk
#'   \item paste
#'   \item gzip
#' }
#' 
#' @param plink.prefix character string with the prefix of binary PLINK files or NULL if bim, bed and fam files are specified separately
#' @param bim.file a character string with with the bim file name or NULL if plink.prefix is given
#' @param fam.file a character string with with the fam file name or NULL if plink.prefix is given
#' @param bed.file a character string with with the bed file name or NULL if plink.prefix is given
#' @param out.prefix a character string with the prefix of the output file (will be named <out.prefix>.data)
#' @param no.web logical. If 'TRUE' PLINK will be called with noweb option
#' @param recode.pheno logical. If 'TRUE' phenotype coding will be changed from 1, 2 to 0, 1
#' for unaffecteds and affecteds, respectively
#' @param zip logical. If 'TRUE' ouput file will be zipped
#' 
#' @export

prepare.input.file.rj<- function(plink.prefix = NULL, bim.file = NULL, bed.file = NULL,
                                  fam.file = NULL, out.prefix, no.web = TRUE,
                                  recode.pheno = TRUE, zip = TRUE) {

    ## check OS
    windows = ifelse(.Platform$OS.type == "windows", TRUE, FALSE)
    
    if (is.null(plink.prefix) & any(is.null(c(bim.file, bed.file, fam.file)))) {
        stop("bim.file, fam.file and bed.file need to be specified!")
    }
    if (!is.null(plink.prefix)) {
        bim.file = paste(plink.prefix, ".bim", sep = "")
        bed.file = paste(plink.prefix, ".bed", sep = "")
        fam.file = paste(plink.prefix, ".fam", sep = "")
    }
    
    ## extract as genotype info
    plink.cmd = "plink"
    if (no.web) plink.cmd = paste(plink.cmd, "--noweb")
    cmd = paste(plink.cmd,
        "--bim", bim.file,
        "--bed", bed.file,
        "--fam", fam.file,
        "--recodeA",
        "--out", out.prefix)
    try(system(cmd, ignore.stdout = TRUE))
    
    ## remove alleles from SNP names
    cmd = paste("sed -i 's/_.\\>//g'",
        paste(out.prefix, "raw", sep = "."))
    try(system(cmd)) 
    
    ## remove first six columns
    cmd = paste("cut -d' ' -f7-", paste(out.prefix, "raw", sep = "."),
        ">", paste(out.prefix, "dat", sep = "."))
    if (windows) { 
      shell(cmd)
    } else { 
      try(system(cmd))
    }  
    
    
    ## extract phenotype
    cmd = paste("cut -d' ' -f6", paste(out.prefix, "raw", sep = "."),
        ">", paste(out.prefix, "pheno", sep = "."))
    if (windows) { 
      shell(cmd)
    } else { 
      try(system(cmd))
    }
    
    ## replace -9 by NA
    cmd = paste("sed -i 's/-9/NA/g'",
        paste(out.prefix, "pheno", sep = "."))
    try(system(cmd))
    
    ## change phenotype coding to 0/1
    if (recode.pheno) {
        cmd = paste("awk '{ if (NR != 1) {$1 = $1 -1}; print}'",
            paste(out.prefix, "pheno", sep = "."),
            ">", paste(out.prefix, "pheno.recoded", sep = "."))
        if (windows) { 
          shell(cmd)
        } else  { 
        try(system(cmd))
        }
    } else {
        cmd = paste("cp", paste(out.prefix, "pheno", sep = "."),
                    paste(out.prefix, "pheno.recoded", sep = "."))
        try(system(cmd))
    }
    
    ## combine
    cmd = paste('paste -d" "', paste(out.prefix, 'pheno.recoded', sep = '.'),
        paste(out.prefix, 'dat', sep = '.'),
        '>', paste(out.prefix, 'data', sep = '.'))
    if (windows) { 
      shell(cmd)
    } else {
    try(system(cmd))
    }

    ## zip
    if (zip) {
        try(system(paste("gzip", paste(out.prefix, "data", sep = "."))))
    
    }
    
    ## tidy up
    temp = file.remove(paste(out.prefix,
        c("log", "raw", "dat", "pheno", "pheno.recoded"), sep = "."))
    
}
