#!/usr/bin/env Rscript

scriptDescription <- "A script that calculates multivariate ROC for metabolites via MetaboAnlyst."

scriptMandatoryArgs <- list(
  inFile = list(
    abbr="-i",
    type="table",
    readoptions=list(sep="\t", stringsAsFactors=FALSE),
    help="Table of metobolite abundances."
  )
)

scriptOptionalArgs <- list(
  outFile = list(
    default="fold_change_summary",
    help="File path without extension to fold change summary figure."
  ),
  fileType = list(
    default="pdf",
    help="File extension for fold change summary figure."
  ),
  figureType = list(
    default="volcano",
    help="Type of plot to save as output."
  ),
  tmpLocation = list(
    default="tmp",
    help="Path to a temporary folder polluted by mSet."
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

for (pk in c("tidyr", "dplyr", "ggplot2", "MetaboAnalystR")){
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

  cat("Parsing dataset\n")
  input <- inFile %>%
    convert_cc_to_mSet(
      tmpLocation=file.path(opt$tmpLocation, "tmp.csv"), analysis_type="roc"
    ) %>%
    normalize_mSet(tmpLocation=opt$tmpLocation)

  if(opt$figureType == "volcano"){

    featureMat <- calcMultiROC(input)

    cat("Drawing volcano plot\n")
    featureMat %>%
      plotROCVolcano() %>%
      fig2pdf(opt$outFile)
    
  } else {
    mSet <- calcMultiROC(
      input, tmpLocation=opt$tmpLocation, outFile=opt$outFile, fileType=opt$fileType
    )
  
  }
  
  invisible(NULL)
}


#' Calculate multivariate ROC
#' 
#' @param norm_data dataframe or mSet. Metabolomics data with Fold Changes and p-values.
#' @param norm_path string. Path to the normalization set.
#' @param tmpLocation string. Path to temporary file for mSet init.
#' @param figureLocation string. Path where figure should be saved (if wanted).
#' @param fileType string. Figure file format.
#' @param keep_mSet logical. If the mSet obeject should be returned or the dataframe only.
#' @param cleanUp logical. If temporary files should be removed after execution.
#' 
#' @return dataframe or mSet  Contains importance of individual metabolites.
calcMultiROC <- function(norm_data, norm_path="tmp/row_norm.qs", tmpLocation="tmp", figureLocation=NULL, fileType="pdf", keep_mSet=FALSE, cleanUp=TRUE){
    
  if(!is(norm_data, "list")){
    stats_data <- norm_data %>%
      convert_cc_to_mSet(tmpLocation=file.path(tmpLocation, "tmp.csv")) %>%
      normalize_mSet(tmpLocation=tmpLocation)
  }
  
  old_wd <- getwd()
  if(!dir.exists(tmpLocation)) dir.create(tmpLocation)
  if(norm_path != file.path(tmpLocation, "row_norm.qs")){
    file.copy(norm_path, file.path(tmpLocation, "row_norm.qs"))
  }
  setwd(tmpLocation)
  
  mSet <- norm_data %>%
    SetAnalysisMode("explore") %>%
    PrepareROCData() %>%
    PerformCV.explore(cls.method = "svm", rank.method = "svm", lvNum = 2)
  
  if(!keep_mSet){
    PlotImpVars( # The side effect of the figure is a table we actually need
      mSet, imgName="multiROC_importance", format=fileType,
      mdl.inx=-1, measure="freq", feat.num=15
    )
    mSet <- read.csv("imp_features_cv.csv")
  }

  setwd(old_wd)

  if(!is.null(figureLocation)){
    "multiROC_importance" %>%
      paste0("72dpi.", fileType) %>%
      file.path(tmpLocation, .) %>%
      file.copy(figureLocation)
  }

  if(cleanUp) unlink(tmpLocation, recursive=TRUE)
  
  invisible(mSet)

}


#' Create a simple Volcano plot with Fold changes and p-values
#' 
#' @param stats_data dataframe. MultiROC data.
#' 
#' @return A ggplot with the Volcano.
plotROCVolcano <- function(stats_data){

  stats_data %>%
    mutate(
      p.value = ifelse(is.na(AUC.value), 1, 1 - AUC.value),
    )
  
}


# Ensuring command line connectivity by sourcing an argument parser
source(opt$commandRpath, local=TRUE)