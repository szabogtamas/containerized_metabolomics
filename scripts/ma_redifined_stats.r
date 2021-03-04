#!/usr/bin/env Rscript

scriptDescription <- "A script that runs descreptive statistics via MetaboAnlyst."

scriptMandatoryArgs <- list(
  inFile = list(
    abbr="-i",
    type="table",
    readoptions=list(sep="\t", stringsAsFactors=FALSE),
    help="Table of metobolite abundances."
  )
)

scriptOptionalArgs <- list(
  conditionOrder = list(
    default=NULL,
    type="vector",
    help="Order of conditions in the experimental design formula. Makes sense to put control as first."
  ),
  outFile = list(
    default="fold_change_summary",
    help="File path without extension to fold change summary figure."
  ),
  fileType = list(
    default="png",
    help="File extension to fold change summary figure."
  ),
  figureType = list(
    default="volcano",
    help="Type of plot to save as output."
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

for (pk in c("tidyr", "dplyr", "MetaboAnalystR")){
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
  
  outFile <- "metabolomics_results"
  opt$outFile <- NULL
  opt$help <- NULL
  opt$verbose <- NULL

  cat("Parsing dataset\n")
  input <- convert_to_mSet(inFile)

  cat("Normalizing dataset\n")
  mSet <- normalize_mSet(input)

  cat("Calculating basic descriptive statistics\n")
  mSet <- FC.Anal(mSet, 2.0, 0)
  
  cat("Saving table\n")
  tab2tsv(mSet$table, outFile)
  
  cat("Saving figure\n")
  if(opt$figureType == "volcano"){
    plotMetaboVolcano(mSet)
  } else {
    tmp_wd <- getwd()
    setwd(dirname(opt$outFile))
    PlotFC(mSet, basename(opt$outFile), opt$fileType)
    setwd(tmp_wd)
    unlink(dirname(opt$outFile))
  }

  invisible(NULL)
}

# Ensuring command line connectivity by sourcing an argument parser
source(opt$commandRpath, local=TRUE)