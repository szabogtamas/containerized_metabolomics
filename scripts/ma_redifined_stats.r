#!/usr/bin/env Rscript

scriptDescription <- "Some functions from Metaboanalyst redefined."

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
  commandRpath = list(
    default="commandR.r",
    help="Path to command line connectivity script (if not in cwd)."
  )
)

opt <- list()
for (rn in names(scriptOptionalArgs)){
  opt[[rn]] <- scriptOptionalArgs[[rn]][["default"]]
}

for (pk in c("tidyr", "dplyr", "dplyr", "MetaboAnalystR")){
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

  cat("Parsing metabolomics input\n")
  mSet <- mSetFromFile(opt$input)

  cat("Calculating basic descriptive statistics\n")
  results <- do.call(calculate_basic_stats, opt)
  
  cat("Saving table\n")
  tab2tsv(results$table, outFile)
  
  cat("Saving figure\n")
  fig2pdf(results$figure, outFile, height=8.64, width=7.2)

  invisible(NULL)
}


#' Calculates basic descriptive statistics for metabolite quantities in conditions.
#' 
#' @param in_df dataframe. Metabolomics data with abundance values and standardized compound names 
#' 
#' @return top results and overview plots.
calculate_basic_stats <- function(in_df, tmp_work_dir){
  
  in_df %>%
    make_mSet_df() %>%
    ReplaceMin() %>%
    FilterVariable("iqr", "F", 25) %>%
    PreparePrenormData() %>%
    Normalization("MedianNorm", "LogNorm", "NULL", ratio=FALSE, ratioNum=20)

  mSet <- PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
  mSet <- PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)

  #mSet <- FC.Anal.unpaired(mSet, 2.0, 0)
  mSet <- FC.Anal(mSet, 2.0, 0)
  mSet <- PlotFC(mSet, "fc_0_", "png", 72, width=NA)
  
}


#' Explores differentially regulated metabolic pathways.
#' 
#' @param in_df dataframe. Metabolomics data with abundance values and standardized compound names 
#' 
#' @return top results and overview plots.
find_de_paths <- function(in_df){
  
  in_df %>%
    msea()
  
}


# Ensuring command line connectivity by sourcing an argument parser
source(opt$commandRpath, local=TRUE)

# Ensuring command line connectivity by sourcing an argument parser
source(opt$commandRpath, local=TRUE)