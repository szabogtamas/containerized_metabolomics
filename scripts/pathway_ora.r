#!/usr/bin/env Rscript

scriptDescription <- "Find top overrepresented pathways in a metabolomics dataset via MetaboAnlyst."

scriptMandatoryArgs <- list(
  hitList= list(
    abbr="-i",
    type="vector",
    readoptions=list(sep="\t", stringsAsFactors=FALSE),
    help="List of top metabolites."
  )
)

scriptOptionalArgs <- list(
  outFile = list(
    default="path_ora",
    help="File path without extension to overrepresented pathways summary figure."
  ),
  fileType = list(
    default="pdf",
    help="File extension for figure."
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
    convert_cc_to_mSet(tmpLocation=file.path(opt$tmpLocation, "tmp.csv")) %>%
    normalize_mSet(tmpLocation=opt$tmpLocation)
  
  cat("Calculating ORA on metabolic pathways\n")
  mSet <- find_metabo_ora(mSet, tmpLocation=opt$tmpLocation, keep_mSet=TRUE)
  
  cat("Saving figure\n")
  if(opt$figureType == "volcano"){
    
    mSet %>%
      plotPathVolcano() %>%
      fig2pdf(opt$outFile)
    
  } else {
    
    tmp_wd <- getwd()
    setwd(dirname(opt$outFile))
    PlotORA(mSet, basename(opt$outFile), opt$fileType)
    setwd(tmp_wd)
    
  }
  
  invisible(NULL)
}


#' Run PATHORA to find overrepresented pathways in mSet
#' 
#' @param hitlist character. Top metabolites.
#' @param tmpLocation string. Path to temporary file for mSet init.
#' @param keep_mSet logical. If the mSet obeject should be returned or the dataframe only.
#' @param cleanUp logical. If temporary files should be removed after execution.
#' 
#' @return dataframe or mSet  Contains descriptive stats.
find_metabo_ora <- function(hitlist, tmpLocation="tmp", keep_mSet=FALSE, cleanUp=TRUE, lipid=FALSE){
  
  old_wd <- getwd()
  if(!dir.exists(tmpLocation)) dir.create(tmpLocation)
  setwd(tmpLocation)
  
  mSet <- NULL
  mSet <- InitDataObjects("conc", "msetora", FALSE) %>%
    Setup.MapData(selected_compounds) %>%
    CrossReferencing("name")
  
  mappable_compunds <- mSet$name.map$query.vec[mSet$name.map$match.state == 1]
  
  mSet <- NULL
  mSet <- InitDataObjects("conc", "msetora", FALSE) %>%
    Setup.MapData(mappable_compunds) %>%
    CreateMappingResultTable() %>%
    SetMetabolomeFilter(FALSE) %>%
    SetCurrentMsetLib("smpdb_pathway", 2) %>%
    CalculateHyperScore()
  
  setwd(old_wd)
  if(cleanUp) unlink(tmpLocation, recursive=TRUE)
  
  invisible(mSet)
  
}

#' Create a simple Volcano plot with Fold changes and p-values
#' 
#' @param stats_data dataframe or mSet. Metabolomics data with Fold Changes and p-values.
#' 
#' @return A ggplot with the Volcano.
plotPathVolcano <- function(stats_data){
  print("TODO")
}


# Ensuring command line connectivity by sourcing an argument parser
source(opt$commandRpath, local=TRUE)