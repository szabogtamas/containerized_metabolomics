#!/usr/bin/env Rscript

block_ora <- TRUE
source("/home/rstudio/repo_files/scripts/pathway_ora.r", local=TRUE)
if (exists("block_kegg")) eval_blocker <- TRUE else eval_blocker <- NULL

scriptDescription <- "Find top enriched KEGG pathways in a metabolomics dataset via MetaboAnlyst."

scriptMandatoryArgs <- list(
  changeValues = list(
    abbr="-i",
    type="tables",
    readoptions=list(stringsAsFactors=FALSE),
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

for (pk in c("tidyr", "dplyr", "purrr", "tibble", "ggplot2", "MetaboAnalystR")){
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
  
  cat("Calculating enrichment of metabolic KEGG pathways\n")
  mSetLs <- opt$changeValues %>%
    map(
      ~find_metabo_kegg(.x, tmpLocation=opt$tmpLocation, keep_mSet=TRUE)
    )
  
  cat("Saving table\n")
  mSummary <- mSetLs %>%
    imap(function(x, y) mutate(x$summary_df, Group = y)) %>%
    bind_rows()
  
  tab2tsv(mSummary, opt$outFile)
  
  cat("Plotting\n")
  if(opt$figureType == "dotplot"){
    
    cat("Saving custom dotplot\n")
    mSummary %>%
      plotPathHits() %>%
      fig2pdf(opt$outFile)
    
  } else {
    
    old_wd <- getwd()
    setwd(dirname(opt$outFile))
    
    for(condition in names(mSetLs)){
      
      mSet <- mSetLs[[condition]]
      mSet <- PlotPathSummary(mSet, TRUE, basename(opt$outFile), opt$fileType)
      
      file.rename(
        paste(basename(opt$outFile), "dpi72.", opt$fileType, sep=""),
        paste(basename(opt$outFile), "_", condition, ".", opt$fileType, sep="")
      )
      
    }
    
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
  
  mappable_compunds <- metabo_change %>%
    .$Metabolite %>%
    unique() %>%
    filter_mappable_compounds()
  
  if(!is(metabo_change, "list")){
    mSet <- metabo_change %>%
      filter(Metabolite %in% mappable_compunds) %>%
      convert_cc_to_mSet(tmpLocation="tmp.csv", analysis_type="pathora") %>%
      normalize_mSet(tmpLocation=".")
  }
  
  mSet <- mSet %>%
    CrossReferencing("name") %>%
    CreateMappingResultTable() %>%
    SetKEGG.PathLib("hsa", "current") %>%
    SetMetabolomeFilter(FALSE) %>%
    CalculateOraScore("rbc", "hyperg")

  mSet$analSet$type <- "pathora" # This information is lost somehow
  
  pw_hit_link <- mSet %>%
    .$analSet %>%
    .$ora.hits %>%
    enframe(name="Pathway", value="Hits") %>%
    unnest(Hits) %>%
    group_by(Pathway) %>%
    mutate(
      Hits = paste(Hits, collapse="; ")
    ) %>%
    ungroup() %>%
    distinct()
  
  pathway_data <- mSet %>%
    .$analSet %>%
    .$ora.mat %>%
    data.frame() %>%
    rownames_to_column("Pathway") %>%
    rename(nHits = Hits) %>%
    mutate(hitRatio = nHits/Total) %>%
    left_join(pw_hit_link, by="Pathway")
  
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
rg <- commandArgs()
if("--commandRpath" %in% rg){
  scriptOptionalArgs$commandRpath$default <- rg[[which(rg == "--commandRpath") + 1]]
}
source(scriptOptionalArgs$commandRpath$default, local=TRUE, chdir=FALSE)