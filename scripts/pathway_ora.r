#!/usr/bin/env Rscript

scriptDescription <- "Find top overrepresented pathways in a metabolomics dataset via MetaboAnlyst."

scriptMandatoryArgs <- list(
  hitList= list(
    abbr="-i",
    type="nested",
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

source("/home/rstudio/repo_files/scripts/data_norm.r", local=TRUE)

#' The main function of the script, executed only if called from command line.
#' Calls subfunctions according to supplied command line arguments.
#' 
#' @param opt list. a named list of all command line options; will be passed on 
#' 
#' @return Not intended to return anything, but rather to save outputs to files.
main <- function(opt){
  
  cat("Calculating ORA on metabolic pathways\n")
  mSet <- find_metabo_ora(opt$hitList, tmpLocation=opt$tmpLocation, keep_mSet=TRUE)
  
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
      PlotORA(mSet, basename(opt$outFile), "bar", opt$fileType)
    } else {
      PlotORA(mSet, basename(opt$outFile), "net", opt$fileType)
    }
    file.rename(
      paste(basename(opt$outFile), "dpi72.", opt$fileType, sep=""),
      paste(basename(opt$outFile), opt$fileType, sep=""),
    )
    
    setwd(old_wd)
    
  }
  
  invisible(NULL)
}


#' Map compound names to literature and select only those that could be mapped
#' 
#' @param hitlist character. Top metabolites.
#' @param tmpLocation string. Path to temporary file for mSet init.
#' @param cleanUp logical. If temporary files should be removed after execution.
#' 
#' @return character. Compounds that could be mapped.
filter_mappable_compounds <- function(hitlist, tmpLocation="tmp", cleanUp=TRUE){
  
  old_wd <- getwd()
  if(!dir.exists(tmpLocation)) dir.create(tmpLocation)
  setwd(tmpLocation)
  
  mSet <- NULL
  mSet <- InitDataObjects("conc", "msetora", FALSE) %>%
    Setup.MapData(selected_compounds) %>%
    CrossReferencing("name")
  
  mappable_compunds <- mSet$name.map$query.vec[mSet$name.map$match.state == 1]
  
  setwd(old_wd)
  if(cleanUp) unlink(tmpLocation, recursive=TRUE)
  
  return(mappable_compunds)
  
}

#' Run PATHORA to find overrepresented pathways in mSet
#' 
#' @param hitlist character. Top metabolites.
#' @param tmpLocation string. Path to temporary file for mSet init.
#' @param keep_mSet logical. If the mSet obeject should be returned or the dataframe only.
#' @param cleanUp logical. If temporary files should be removed after execution.
#' 
#' @return dataframe or mSet  Contains descriptive stats.
find_metabo_ora <- function(hitlist, tmpLocation="tmp", keep_mSet=FALSE, cleanUp=TRUE){
  
  old_wd <- getwd()
  if(!dir.exists(tmpLocation)) dir.create(tmpLocation)
  setwd(tmpLocation)
  
  mappable_compunds <- filter_mappable_compounds(hitlist)
  
  mSet <- NULL
  mSet <- InitDataObjects("conc", "msetora", FALSE) %>%
    Setup.MapData(mappable_compunds) %>%
    CrossReferencing("name") %>%
    CreateMappingResultTable() %>%
    SetMetabolomeFilter(FALSE) %>%
    SetCurrentMsetLib("smpdb_pathway", 2) %>%
    CalculateHyperScore()
  
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
    rename(nHits = hits) %>%
    mutate(hitRatio = nHits/total) %>%
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

#' Create a dotplot showing how significant the association between hits and top pathways is
#' 
#' @param paths_data dataframe or mSet. Metabolomics data with hits on a pathway and p-values.
#' @param numPath integer. Number of pathways to be shown on plot.
#' 
#' @return A ggplot with the dotplot showing association with top pathways.
plotPathHits <- function(paths_data, numPath=15){
  
  if(!("Group" %in% colnames(paths_data))) paths_data$Group <- ""
  
  paths_data %>%
    arrange(Raw.p) %>%
    head(numPath) %>%
    mutate(
      Pathway = factor(Pathway, levels=rev(unique(.$Pathway)))
    ) %>%
    ggplot(aes(x=Group, y=Pathway, size=hitRatio, color=Raw.p)) +
    geom_point() +
    scale_color_gradientn(
      colors=rev(c('#2b8cbe', 'grey', '#e38071', '#e34a33', '#e31e00')),
      breaks=c(0.05, 0.01, 0.001, 0.0001),
      limits=c(0.00001, 1), trans='log10', oob = scales::squish
    ) +
    theme_bw() +
    theme(
      axis.text.x=element_text(angle=30, hjust=1),
      axis.text.y=element_text(size=10)
    ) +
    labs(x="", y="")
}


# Ensuring command line connectivity by sourcing an argument parser
source(opt$commandRpath, local=TRUE)