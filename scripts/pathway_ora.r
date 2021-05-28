#!/usr/bin/env Rscript

block_norm <- TRUE
source("/home/rstudio/repo_files/scripts/data_norm.r", local=TRUE)
if (exists("block_ora")) eval_blocker <- TRUE else eval_blocker <- NULL

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
  conditionLabels = list(
    default=NULL,
    type="vector",
    help="Labels to be associated with input tables."
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

for (pk in c("pryr", "pkgcond", "tidyr", "dplyr", "tibble", "ggplot2", "MetaboAnalystR")){
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
  
  cat("Calculating ORA on metabolic pathways\n")
  ora_paths <- list()

  if(is.null(opt$conditionLabels)){
    if(length(opt$hitList) == 1) names(opt$hitList) <- "_"
  } else {
    for(nm in names(opt$hitList)){
      opt$hitList[nm] <- c(nm, opt$hitList[nm])
    }
    names(opt$hitList) <- opt$conditionLabels
  }
  
  for (condition in names(opt$hitList)){
    
    hitlist <- opt$hitList[[condition]]

    if(!is.null(hitlist)){
      
      cat("Calculating ORA for", condition, "\n")
    
      if(condition == "_"){
        condition <- ""
        outFile <- opt$outFile
      } else {
        outFile <- paste(opt$outFile, condition, sep="_")
      }
      
      mSet <- find_metabo_ora(hitlist, tmpLocation=opt$tmpLocation, keep_mSet=TRUE)
      
      ora_paths[[condition]] <- mSet %>%
        .$summary_df %>%
        mutate(Group = condition)
      
      if(opt$figureType != "dotplot"){
        
        cat("Saving built-in figure for", condition, "\n")
        
        old_wd <- getwd()
        setwd(dirname(outFile))
        
        if(opt$figureType == "bar"){
          PlotORA(mSet, basename(outFile), "bar", opt$fileType)
        } else {
          PlotORA(mSet, basename(outFile), "net", opt$fileType)
        }
        file.rename(
          paste(basename(outFile), "dpi72.", opt$fileType, sep=""),
          paste(basename(outFile), opt$fileType, sep=".")
        )
        
        setwd(old_wd)
      }
    }
  }
    
  if(length(ora_paths) > 1){
    ora_paths <- bind_rows(ora_paths)
  } else {
    ora_paths <- ora_paths[[1]]
  }
    
  if(nrow(ora_paths) > 0){
    if(opt$figureType == "dotplot"){
      
      cat("Saving figure\n")
      ora_paths %>%
        plotPathHits() %>%
        fig2pdf(opt$outFile)
      
    }
  } else {

    ora_paths <- data.frame(Pathway=c(), hitRatio=c(), Raw.p=c()) 

    cat("Empty figure\n")
      p <- ggplot() + theme_void()
      fig2pdf(p, opt$outFile)
  }
    
  cat("Saving table\n")
  tab2tsv(ora_paths, opt$outFile)
      
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
  print(hitlist)
  old_wd <- getwd()
  if(!dir.exists(tmpLocation)) dir.create(tmpLocation)
  setwd(tmpLocation)
  
  mSet <- NULL
  mSet <- InitDataObjects("conc", "msetora", FALSE) %>%
    Setup.MapData(hitlist) %>%
    CrossReferencing("name")
  
  if(!is.null(mSet$name.map)){
    mappable_compunds <- mSet$name.map$query.vec[mSet$name.map$match.state == 1]
  } else {
    mappable_compunds <- c()
  }
  
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
  
  mappable_compunds <- suppress_messages(filter_mappable_compounds(hitlist))
  
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
    
  manenv <- print(where("current.msetlib"))
  rm("current.msetlib", envir=manenv)
  
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
rg <- commandArgs()
if("--commandRpath" %in% rg){
  scriptOptionalArgs$commandRpath$default <- rg[[which(rg == "--commandRpath") + 1]]
}
source(scriptOptionalArgs$commandRpath$default, local=TRUE, chdir=FALSE)