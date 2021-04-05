#!/usr/bin/env Rscript

scriptDescription <- "Find top enriched KEGG pathways in a metabolomics dataset via MetaboAnlyst."

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
    default="path_kegg",
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
    default="/home/rstudio/repo_files/scripts/commandR.r",
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
  
  cat("Calculating enrichment of metabolic KEGG pathways\n")
  mSet <- find_metabo_kegg(changeValues, tmpLocation=opt$tmpLocation, keep_mSet=TRUE)
  
  cat("Saving figure\n")
  if(opt$figureType == "dotplot"){
    
    mSet %>%
      plotPathHits() %>%
      fig2pdf(opt$outFile)
    
    file.rename
    
  } else {
    
    old_wd <- getwd()
    setwd(dirname(opt$outFile))
    
    PlotPathSummary(mSet, basename(opt$outFile), opt$fileType)
    #PlotKEGGPath(mSet, "Glycine, serine and threonine metabolism", 528, 480, "png", NULL)
    
    setwd(old_wd)
    
  }
  
  invisible(NULL)
}


#' Run msea on KEGG to find overrepresented pathways in mSet
#' 
#' @param hitlist character. Top metabolites.
#' @param tmpLocation string. Path to temporary file for mSet init.
#' @param keep_mSet logical. If the mSet obeject should be returned or the dataframe only.
#' @param cleanUp logical. If temporary files should be removed after execution.
#' 
#' @return dataframe or mSet  Contains descriptive stats.
find_metabo_kegg <- function(metabo_change, tmpLocation="tmp", keep_mSet=FALSE, cleanUp=TRUE){
  
  old_wd <- getwd()
  if(!dir.exists(tmpLocation)) dir.create(tmpLocation)
  setwd(tmpLocation)
  
  mappable_compunds <- filter_mappable_compounds(hitlist)
  
  if(!is(metabo_change, "list")){
    mSet <- metabo_change %>%
      filter(Metabolite %in% mappable_compunds) %>%
      convert_cc_to_mSet(tmpLocation="tmp.csv", analysis_type="pathora") %>%
      normalize_mSet(tmpLocation=tmpLocation)
  }
  
  pathway_data <- mSet %>%
    CrossReferencing("name") %>%
    CreateMappingResultTable() %>%
    SetKEGG.PathLib("hsa", "current") %>%
    SetMetabolomeFilter(FALSE) %>%
    CalculateOraScore("rbc", "hyperg")
  
  
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