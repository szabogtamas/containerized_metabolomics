#!/usr/bin/env Rscript

scriptDescription <- "A script that parses metabolomics measurement results."

scriptMandatoryArgs <- list(
  inFile = list(
    abbr="-i",
    type="table",
    readoptions=list(sep="\t", stringsAsFactors=FALSE),
    help="Table of metabolite abundances in the MetaboAnalyser standard format."
  )
)

scriptOptionalArgs <- list(
  outFile = list(
    default="normalization_summary",
    help="File path without extension to normalization summary figure."
  ),
  fileType = list(
    default="png",
    help="File extension to normalization summary figure."
  ),
  commandRpath = list(
    default="commandR.r",
    help="Path to command line connectivity script (if not in cwd)."
  )
)

opt <- list()
for (rn in names(scriptOptionalArgs)){
  opt[[rn]] <- scriptOptionalArgs[[rn]][["default"]]
}

for (pk in c("tidyr", "dplyr")){
  if(!(pk %in% (.packages()))){
    library(pk, character.only=TRUE)
  }
}

#' The main function of the script, executed only if called from command line.
#' Calls subfunctions according to supplied command line arguments.
#' 
#' @param opt list. a named list of all command line options; will be passed on 
#' 
#' @return Not intended to return anything, but rather to save outputs to files.
main <- function(opt){

  cat("Parsing dataset\n")
  input <- convert_to_mSet(inFile)

  cat("Normalizing dataset\n")
  results <- normalize_mSet(input)

  cat("Saving figure\n")
  tmp_wd <- getwd()
  setwd(dirname(opt$outFile))
  #PlotNormSummary(results, basename(opt$outFile), opt$fileType)
  PlotSampleNormSummary(results, basename(opt$outFile), opt$fileType)
  setwd(tmp_wd)

  invisible(NULL)
}


#' Preprocess and normalize data inside a MetaboAnalyst mSet object
#' 
#' @param inSet mSet. a proprietary object of MetaboAnalyst after sanitycheck
#' 
#' @return normalized metabo Set.
normalize_mSet <- function(inSet){

  inSet %>%
    ReplaceMin() %>%
    FilterVariable("iqr", "F", 25) %>%
    PreparePrenormData() %>%
    Normalization("MedianNorm", "LogNorm", "NULL", ratio=FALSE, ratioNum=20)

}

# Ensuring command line connectivity by sourcing an argument parser
source(opt$commandRpath, local=TRUE)