#!/usr/bin/env Rscript

scriptDescription <- "Some functions from Metaboanalyst redefined."

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
  commandRpath = list(
    default="commandR.r",
    help="Path to command line connectivity script (if not in cwd)."
  )
)

opt <- list()
for (rn in names(scriptOptionalArgs)){
  opt[[rn]] <- scriptOptionalArgs[[rn]][["default"]]
}

for (pk in c("tidyr", "dplyr", "dplyr", "MetaboAnalystR")){
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
  
  outFile <- "metabolomics_results"
  opt$outFile <- NULL
  opt$help <- NULL
  opt$verbose <- NULL

  cat("Just a test for now\n")
  results <- do.call(print, opt)
  
  cat("Saving table\n")
  tab2tsv(results$table, outFile)
  
  cat("Saving figure\n")
  fig2pdf(results$figure, outFile, height=8.64, width=7.2)

  invisible(NULL)
}


#' Carries out enrichment analysis of metabolites.
#' 
#' @param in_df dataframe. Metabolomics data with abundance values and standardized compound names 
#' 
#' @return top results and overview plots.
plot_metabo_enrichment <- function(in_df){
  
mSet <- in_df %>%
  mSet()

mSet <- mSet %>%
  PCA.Anal() %>%
  PlotPCAPairSummary("pca_pair_0_", "png", 72, width=NA, 5) %>%
  PlotPCAScree("pca_scree_0_", "png", 72, width=NA, 5) %>%
  PlotPCA2DScore("pca_score2d_0_", "png", 72, width=NA, 1,2,0.95,1,0) %>%
  PlotPCALoading("pca_loading_0_", "png", 72, width=NA, 1,2) %>%
  PlotPCABiplot("pca_biplot_0_", "png", 72, width=NA, 1,2) %>%
  PlotPCA3DScoreImg("pca_score3d_0_", "png", 72, width=NA, 1,2,3, 40)

mSet <- mSet %>%
  PLSR.Anal(reg=TRUE) %>%
  PlotPLSPairSummary("pls_pair_0_", "png", 72, width=NA, 5) %>%
  PlotPLS2DScore("pls_score2d_0_", "png", 72, width=NA, 1,2,0.95,1,0) %>%
  PlotPLS3DScoreImg("pls_score3d_0_", "png", 72, width=NA, 1,2,3, 40) %>%
  PlotPLSLoading("pls_loading_0_", "png", 72, width=NA, 1, 2) %>%
  PLSDA.CV("L",5, "Q2") %>%
  PlotPLS.Classification("pls_cv_0_", "png", 72, width=NA) %>%
  PlotPLS.Imp("pls_imp_0_", "png", 72, width=NA, "vip", "Comp. 1", 15, FALSE)
  
}


#' Explores differentially regulated metabolic pathways.
#' 
#' @param in_df dataframe. Metabolomics data with abundance values and standardized compound names 
#' 
#' @return top results and overview plots.
find_de_paths <- function(in_df){
  
  in_df %>%
    msea()
  
}


# Ensuring command line connectivity by sourcing an argument parser
source(opt$commandRpath, local=TRUE)

# Ensuring command line connectivity by sourcing an argument parser
source(opt$commandRpath, local=TRUE)