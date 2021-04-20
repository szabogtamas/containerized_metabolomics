#!/usr/bin/env Rscript

if (exists("eval_blocker")) eval_blocker <- 2 else eval_blocker <- 1
source("/home/rstudio/repo_files/scripts/data_norm.r", local=TRUE)
if (eval_blocker == 1) eval_blocker <- NULL else if (eval_blocker == 2) eval_blocker <- 1

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
    default="boxes",
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

for (pk in c("tidyr", "dplyr", "tibble", "ggplot2", "MetaboAnalystR")){
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
  input <- opt$inFile %>%
    convert_cc_to_mSet(
      tmpLocation=file.path(opt$tmpLocation, "tmp.csv"), analysis_type="roc"
    ) %>%
    normalize_mSet(tmpLocation=opt$tmpLocation)
  
  if(opt$figureType == "boxes"){
    
    featureMat <- calcMultiROC(input)
    
    cat("Drawing ggplot\n")
    featureMat %>%
      plotROCfeat() %>%
      fig2pdf(opt$outFile)
    
  } else {
    mSet <- calcMultiROC(
      input, tmpLocation=opt$tmpLocation, figureLocation=opt$outFile, fileType=opt$fileType
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
    SetAnalysisMode("explore") %>%
    PrepareROCData() %>%
    PerformCV.explore(cls.method = "svm", rank.method = "svm", lvNum = 2) %>%
    PlotImpVars( # The side effect of the figure is a table we actually need
      imgName="multiROC_importance_", format=fileType,
      mdl.inx=-1, measure="freq", feat.num=15
    )
  
  if(is.null(mSet$dataSet$orig.var.nms)){
    metab_names <- data.frame(
      Metabolite_std = mSet$dataSet$cmpd,
      Metabolite = mSet$dataSet$cmpd,
      stringsAsFactors = FALSE
    )
  } else {
    metab_names <- mSet %>%
      .$dataSet %>%
      .$orig.var.nms %>%
      enframe(name="Metabolite_std", value="Metabolite") %>%
      mutate(Metabolite_std = as.character(Metabolite_std))
  }
  
  rocStats <- "imp_features_cv.csv" %>%
    read.csv() %>%
    select(Metabolite_std=X, Importance) %>%
    mutate(Metabolite_std = as.character(Metabolite_std)) %>%
    left_join(metab_names, by="Metabolite_std")
  
  sample_class_labels <- data.frame(
    Sample = rownames(mSet$dataSet$norm),
    Condition = mSet$dataSet$cls
  )
  
  rocSummary <- mSet %>%
    .$dataSet %>%
    .$norm %>%
    t() %>%
    `colnames<-`(paste("SAMPLE", colnames(.), sep="_")) %>%
    data.frame() %>%
    rownames_to_column("Metabolite_std") %>%
    pivot_longer(-Metabolite_std) %>%
    rename(Sample=name, Normalized_value=value) %>%
    mutate(Sample = str_replace(Sample, "^SAMPLE_", "")) %>%
    left_join(sample_class_labels, by="Sample") %>%
    right_join(rocStats, by="Metabolite_std")
  
  if(!keep_mSet){
    mSet <- rocSummary
  }
  
  setwd(old_wd)
  
  if(!is.null(figureLocation)){
    "multiROC_importance_" %>%
      paste0("dpi72.", fileType) %>%
      file.path(tmpLocation, .) %>%
      file.copy(paste(figureLocation, fileType, sep="."))
  }
  
  if(cleanUp) unlink(tmpLocation, recursive=TRUE)
  
  invisible(mSet)
  
}


#' Show distribution of top metabolites in multiROC by group
#' 
#' @param stats_data dataframe. MultiROC data.
#' @param num_feat integer. Number of features to be shown on plot.
#' 
#' @return A ggplot showing distributions.
plotROCfeat <- function(stats_data, num_feat=15){
  
  top_mroc_mtb <- stats_data %>%
    distinct(Metabolite, Importance) %>%
    arrange(desc(Importance)) %>%
    head(num_feat) %>%
    .$Metabolite
  
  stats_data %>%
    filter(Metabolite %in% top_mroc_mtb) %>%
    mutate(
      Metabolite = factor(Metabolite, levels=top_mroc_mtb)
    ) %>%
    ggplot(aes(x = Metabolite, y = Normalized_value, color = Condition)) + 
    geom_boxplot(outlier.shape=NA) +
    geom_point(size=1, position=position_jitterdodge()) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle=30, hjust=1)
    ) +
    labs(x = "", y = "Normalized value")
  
}

# Ensuring command line connectivity by sourcing an argument parser
rg <- commandArgs()
if("--commandRpath" %in% rg){
  scriptOptionalArgs$commandRpath$default <- rg[[which(rg == "--commandRpath") + 1]]
}
source(scriptOptionalArgs$commandRpath$default, local=TRUE, chdir=FALSE)