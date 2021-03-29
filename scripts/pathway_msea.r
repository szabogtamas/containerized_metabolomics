#!/usr/bin/env Rscript

scriptDescription <- "Find top enriched pathways in a metabolomics dataset via MetaboAnlyst."

scriptMandatoryArgs <- list(
  changeValues = list(
    abbr="-i",
    type="table",
    readoptions=list(sep="\t", stringsAsFactors=FALSE),
    help="Metabolite names and a change measure or score."
  )
)

scriptOptionalArgs <- list(
  outFile = list(
    default="path_msea",
    help="File path without extension to enriched pathways summary figure."
  ),
  fileType = list(
    default="pdf",
    help="File extension for figure."
  ),
  figureType = list(
    default="dotplot",
    help="Type of plot to save as output."
  ),
  tmpLocation = list(
    default="tmp",
    help="Path to a temporary folder where mSet side effects are dumped."
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

for (pk in c("tidyr", "dplyr", "tibble", "ggplot2", "MetaboAnalystR")){
  if(!(pk %in% (.packages()))){
    library(pk, character.only=TRUE)
  }
}

source("pathway_ora.r", local=TRUE)

#' The main function of the script, executed only if called from command line.
#' Calls subfunctions according to supplied command line arguments.
#' 
#' @param opt list. a named list of all command line options; will be passed on 
#' 
#' @return Not intended to return anything, but rather to save outputs to files.
main <- function(opt){
  
  hitList <- opt$hitList
  
  cat("Calculatingenrichment of metabolic pathways\n")
  mSet <- find_metabo_msea(changeValues, tmpLocation=opt$tmpLocation, keep_mSet=TRUE)
  
  cat("Saving figure\n")
  if(opt$figureType == "dotplot"){
    
    mSet %>%
      plotPathHits() %>%
      fig2pdf(opt$outFile)
    
    file.rename
    
  } else {
    
    old_wd <- getwd()
    setwd(dirname(opt$outFile))
    
    if(opt$figureType == "bar"){
      PlotQEA.Overview(mSet, basename(opt$outFile), "bar", opt$fileType)
    } else {
      PlotQEA.Overview(mSet, basename(opt$outFile), "net", opt$fileType)
    }
    file.rename(
      paste(basename(opt$outFile), "dpi72.", opt$fileType, sep=""),
      paste(basename(opt$outFile), opt$fileType, sep=""),
    )
    
    setwd(old_wd)
    
  }
  
  invisible(NULL)
}


#' Run msea/QEA to find enriched pathways in mSet
#' 
#' @param hitlist character. Top metabolites.
#' @param tmpLocation string. Path to temporary file for mSet init.
#' @param keep_mSet logical. If the mSet obeject should be returned or the dataframe only.
#' @param cleanUp logical. If temporary files should be removed after execution.
#' 
#' @return dataframe or mSet  Contains descriptive stats.
find_metabo_msea <- function(metabo_change, tmpLocation="tmp", keep_mSet=FALSE, cleanUp=TRUE){
  
  old_wd <- getwd()
  if(!dir.exists(tmpLocation)) dir.create(tmpLocation)
  setwd(tmpLocation)
  
  mappable_compunds <- filter_mappable_compounds(hitlist)
  
  if(!is(metabo_change, "list")){
    mSet <- metabo_change %>%
      convert_cc_to_mSet(tmpLocation="tmp.csv", analysis_type="msetqea") %>%
      normalize_mSet(tmpLocation=tmpLocation)
  }
  
  mSet <- mSet %>%
    CrossReferencing("name") %>%
    CreateMappingResultTable() %>%
    ReplaceMin() %>%
    PreparePrenormData() %>%
    Normalization("NULL", "NULL", "NULL", "PIF_178", ratio=FALSE, ratioNum=20) %>%
    SetMetabolomeFilter(FALSE) %>%
    SetCurrentMsetLib("smpdb_pathway", 2) %>%
    CalculateGlobalTestScore()
  
  if(keep_mSet){
    mSet$summary_df <- pathway_data
  } else {
    mSet <- pathway_data
  }
  
  setwd(old_wd)
  if(cleanUp) unlink(tmpLocation, recursive=TRUE)
  
  invisible(mSet)
  
}


# Ensuring command line connectivity by sourcing an argument parser
source(opt$commandRpath, local=TRUE)