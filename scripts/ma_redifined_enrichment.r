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
  
  in_df %>%
    Enrich.Anal()

    # Create mSetObj for storing objects created during your analysis
mSet<-InitDataObjects("conc", "pathora", FALSE)
# Set up mSetObj with the list of compounds
mSet<-Setup.MapData(mSet, tmp.vec);
# Cross reference list of compounds against libraries (hmdb, pubchem, chebi, kegg, metlin)
mSet<-CrossReferencing(mSet, "name");
# Creates a mapping result table; shows HMDB, KEGG, PubChem, etc. IDs
# Saved as "name_map.csv" or can be found in mSet$dataSet$map.table
# Compounds with no hits will contain NAs across the columns
mSet<-CreateMappingResultTable(mSet);
# From the mapping result table, L-Isolucine has no matches
# Now, perform potential matching with our database against this compound
mSet<-PerformDetailMatch(mSet, "L-Isolucine");
# Get list of candidates for matching
# Results are found in mSet$name.map$hits.candidate.list
mSet<-GetCandidateList(mSet);
# Replace L-Isolucine with selected compound (L-Isoleucine)
mSet<-SetCandidate(mSet, "L-Isolucine", "L-Isoleucine");

# Use "current" for the latest KEGG pathway library or "v2018" for the KEGG pathway library version prior to November 2019.
mSet<-SetKEGG.PathLib(mSet, "hsa", "current")
# Set the metabolite filter
# Default set to false
mSet<-SetMetabolomeFilter(mSet, F);
# Calculate the over representation analysis score, here we selected to use the hypergeometric test (alternative is Fisher's exact test)
# A results table "pathway_results.csv" will be created and found within your working directory
mSet<-CalculateOraScore(mSet, "rbc", "hyperg")
# Plot of the Pathway Analysis Overview
mSet<-PlotPathSummary(mSet, "path_view_0_", "png", 72, width=NA)
# Plot a specific metabolic pathway, in this case "Glycine, serine and threonine metabolism"
mSet<-PlotKEGGPath(mSet, "Glycine, serine and threonine metabolism",528, 480, "png", NULL)
  
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