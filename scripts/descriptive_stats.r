#!/usr/bin/env Rscript

eval_blocker <- 1
source("/home/rstudio/repo_files/scripts/data_norm.r", local=TRUE)
eval_blocker <- NULL

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
  tmpLocation = list(
    default="tmp",
    help="Path to a temporary folder polluted by mSet."
  ),
  commandRpath = list(
    default="/home/rstudio/repo_files/scripts/commandR.r",
    help="Path to command line connectivity script (if not in cwd)."
  )
)

for (pk in c("tidyr", "dplyr", "ggplot2", "MetaboAnalystR")){
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

  cat("Parsing dataset\n")
  input <- inFile %>%
    convert_cc_to_mSet(tmpLocation=file.path(opt$tmpLocation, "tmp.csv")) %>%
    normalize_mSet(tmpLocation=opt$tmpLocation)

  cat("Calculating basic descriptive statistics\n")
  mSet <- calcMetaboStat(mSet, tmpLocation=opt$tmpLocation, keep_mSet=TRUE)
  
  cat("Saving table\n")
  stats_data <- extract_stat_from_mSet(stats_data)
  tab2tsv(stats_data, opt$outFile)
  
  cat("Saving figure\n")
  if(opt$figureType == "volcano"){
    
    stats_data %>%
      plotMetaboVolcano() %>%
      fig2pdf(opt$outFile)
    
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
#' @param tmpLocation string. Path to temporary file for mSet init.
#' @param keep_mSet logical. If the mSet obeject should be returned or the dataframe only.
#' @param cleanUp logical. If temporary files should be removed after execution.
#' 
#' @return dataframe or mSet  Contains descriptive stats.
calcMetaboStat <- function(norm_data, norm_path="tmp/row_norm.qs", tmpLocation="tmp", keep_mSet=FALSE, cleanUp=TRUE){
    
  if(!is(norm_data, "list")){
    norm_data <- norm_data %>%
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
    FC.Anal(2.0, 0) %>%
    Ttests.Anal(F, 0.05, FALSE, TRUE)
  
  if(!keep_mSet) mSet <- extract_stat_from_mSet(mSet)
  
  setwd(old_wd)
  if(cleanUp) unlink(tmpLocation, recursive=TRUE)
  
  invisible(mSet)

}


#' Extract descriptive stats from an Mset object into a neat dataframe
#' 
#' @param mSet dataframe or mSet. Metabolomics data with Fold Changes and p-values.
#' 
#' @return dataframe Containing FC and p-values.
extract_stat_from_mSet <- function(mSet){
  
  fc_df <- mSet %>%
    .$analSet %>%
    .$fc %>%
    .$sig.mat %>%
    data.frame() %>%
    tibble::rownames_to_column("Metabolite")
  
  colnames(fc_df) <- c("Metabolite", "FC", "logFC")
  
  tt_df <- mSet %>%
    .$analSet %>%
    .$tt %>%
    .$sig.mat %>%
    data.frame() %>%
    tibble::rownames_to_column("Metabolite")
  
  colnames(tt_df) <- c("Metabolite", "t.stat", "p.value", "logP", "FDR")
  
  full_join(fc_df, tt_df, by="Metabolite")
  
}


#' Create a simple Volcano plot with Fold changes and p-values
#' 
#' @param stats_data dataframe or mSet. Metabolomics data with Fold Changes and p-values.
#' 
#' @return A ggplot with the Volcano.
plotMetaboVolcano <- function(stats_data){
  
  if(is(stats_data, "list")){
    stats_data <- extract_stat_from_mSet(stats_data)
  }
  
  stats_data %>%
    mutate(
      p.value = ifelse(is.na(p.value), 1, p.value),
      Effect = case_when(
        abs(log2(FC)) < log2(2) ~ "Small change",
        FDR < 0.05 ~ "Significant change",
        TRUE ~ "Insignificant change"
      ),
      Effect = factor(Effect, levels=c("Significant change", "Small change", "Insignificant change"))
    ) %>%
    ggplot(aes(x=log2(FC), y=-log10(p.value), color=Effect, label=Metabolite)) +
    geom_point(size=2) +
    scale_color_manual(values=c("#E64B35B2", "#4DBBD5B2", "#00A087B2"), drop=FALSE) +
    theme_bw() +
    labs(
      title = "Volcano plot of results",
      x = "log2(FC)", y = "-log10(p)"
    )
  
}

# Ensuring command line connectivity by sourcing an argument parser
rg <- commandArgs()
if("--commandRpath" %in% rg){
  scriptOptionalArgs$commandRpath$default <- rg[[which(rg == "--commandRpath") + 1]]
}
source(scriptOptionalArgs$commandRpath$default, local=TRUE, chdir=FALSE)