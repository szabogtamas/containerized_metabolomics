#!/usr/bin/env Rscript

scriptDescription <- "A script that runs descriptive statistics via MetaboAnlyst."

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
    default="pdf",
    help="File extension to fold change summary figure."
  ),
  figureType = list(
    default="volcano",
    help="Type of plot to save as output."
  ),
  commandRpath = list(
    default="/home/rstudio/git_repo/scripts/commandR.r",
    help="Path to command line connectivity script (if not in cwd)."
  )
)

if(!exists("opt")){
  opt <- list()
}

rg <- commandArgs()
if("--commandRpath" %in% rg){
  opt$commandRpath <- rg[[which(rg == "--commandRpath") + 1]]
}

opt <- list()
for (rn in names(scriptOptionalArgs)){
  opt[[rn]] <- scriptOptionalArgs[[rn]][["default"]]
}

for (pk in c("tidyr", "dplyr", "MetaboAnalystR")){
  if(!(pk %in% (.packages()))){
    library(pk, character.only=TRUE)
  }
}

source("data_norm.r", local=TRUE)

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
  mSet <- calcMetaboStat(mSet, keep_mSet=TRUE)
  
  cat("Saving table\n")
  tab2tsv(mSet$table, outFile)
  
  cat("Saving figure\n")
  if(opt$figureType == "volcano"){
    plotMetaboVolcano(mSet)
  } else {
    tmp_wd <- getwd()
    setwd(dirname(opt$outFile))
    if(opt$figureType == "FCscatter"){
      PlotFC(mSet, basename(opt$outFile), opt$fileType)
    } else {
      PlotTT(mSet, basename(opt$outFile), opt$fileType)
    }
    setwd(tmp_wd)
  }

  invisible(NULL)
}


#' Calculate basic escriptive statistics, like pairwise t-tests and fold changes on mSet
#' 
#' @param norm_data dataframe or mSet. Metabolomics data with Fold Changes and p-values.
#' @param norm_path string. Path to the normalization set.
#' @param keep_mSet logical. If the mSet obeject should be returned or the dataframe only.
#' 
#' @return Not intended to return anything, but rather to save outputs to files.
calcMetaboStat <- function(norm_data, norm_path="row_norm.qs", keep_mSet=FALSE){
  
  if(norm_path != "row_norm.qs"){
    file.copy(norm_path, "row_norm.qs")
  }
  
  mSet <- norm %>%
    FC.Anal(2.0, 0) %>%
    Ttests.Anal(F, 0.05, FALSE, TRUE)
  
  if(!keep_mSet){
    
    fc_df <- stat_results %>%
      .$analSet %>%
      .$fc %>%
      .$sig.mat %>%
      data.frame() %>%
      tibble::rownames_to_column("Metabolite")
    
    colnames(fc_df) <- c("Metabolite", "FC", "logFC")
    
    tt_df <- stat_results %>%
      .$analSet %>%
      .$tt %>%
      .$sig.mat %>%
      data.frame() %>%
      tibble::rownames_to_column("Metabolite")
    
    colnames(tt_df) <- c("Metabolite", "t.stat", "p.value", "logP", "FDR")
    
    mSet <- full_join(fc_df, tt_df, by="Metabolite")
    
  }
  
  unlink(norm_path)
  unlink("row_norm.qs")
  
  invisible(mSet)

}


#' Create a simple Volcano plot with Fold changes and p-values
#' 
#' @param stats_data dataframe or mSet. Metabolomics data with Fold Changes and p-values.
#' 
#' @return Not intended to return anything, but rather to save outputs to files.
plotMetaboVolcano <- function(stats_data){
  
  if(is(stats_data, "mSet")){
    stats_data <- stats_data$data
  }

  stats_data %>%
    ggplot(aes(x=FC, y=pvalue)) +
    geom_point(size=2)

}

# Ensuring command line connectivity by sourcing an argument parser
source(opt$commandRpath, local=TRUE)