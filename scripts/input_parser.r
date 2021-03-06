#!/usr/bin/env Rscript

scriptDescription <- "A script that parses metabolomics measurement results for MetaboAnalyst."

scriptMandatoryArgs <- list(
  inFile = list(
    abbr="-i",
    type="table",
    readoptions=list(sep="\t", stringsAsFactors=FALSE),
    help="Table of metabolite quantities/concentrations."
  )
)

scriptOptionalArgs <- list(
  tmpLocation = list(
    default="tmp.csv",
    help="Path to the temporary file used to push results in an mSet object."
  ),
  analysis_type = list(
    default="stat",
    help="Analysis type added when mSet object is initialized."
  ),
  input_format = list(
    default="rowu",
    help="Format code used by MetaboAnalyser to parse input."
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

for (pk in c("tidyr", "dplyr", "stringr", "MetaboAnalystR")){
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

  cat("Standardizing input\n")
  results <- standardize_metabo_data(opt$inFile)
  
  cat("Saving to temporary MetoboAnalyst file\n")
  tmp_wd <- getwd()
  setwd(dirname(opt$tmpLocation))
  
  write_metabodf_tmp(results, opt$tmpLocation)
  mSet <- convert_cc_to_mSet(opt$tmpLocation, opt$analysis_type, opt$input_format)
  
  setwd(tmp_wd)
  unlink(dirname(opt$tmpLocation))
  
  invisible(mSet)
}


#' Convert primary input to a standardized, tidy format
#' 
#' @param input_data dataframe. The primary input
#' 
#' @return Not intended to return anything, but rather to save outputs to files.
standardize_metabo_data <- function(input_data){

  input_data %>%
    rename(Metabolite = X, MetaboGroup = X.1) %>%
    pivot_longer(
      c(-Metabolite, -MetaboGroup), names_to="Condition", names_repair="universal"
    ) %>%
    mutate(
      Replicate = paste("BATCH", str_extract(Condition, "\\.\\d+")),
      Condition = str_replace(Condition, "\\.\\d+", ""),
      Subject = paste(Condition, Replicate, sep="_")
    ) %>%
    select(-name)
    
}


#' Save a standardized dataframe to a temporary file that can be read in as an mSet object
#' 
#' @param metab_data datafram. Path to the preformatted input file. Cannot be dataframe unfortunately.
#' @param subject character. Column name for subject identifiers 
#' @param condition character. Column name for condition  
#' @param metabolite character. Column name for metabolite identifiers 
#' @param value character. Column name for metabolite quantity (concentration) values
#' @param tmp_out_file character. Path to the temporary output file 
#' 
#' @return Not intended to return anything, but rather to save output to file.
write_metabodf_tmp <- function(metab_data, subject="Subject", condition="Condition", metabolite="Metabolite", value="value", tmp_out_file="tmp.csv"){
  
  metab_data %>%
    select(
      !!symbol(subject), !!symbol(condition), !!symbol(metabolite), !!symbol(value)
    )  %>%
    rename(
      Subject = !!symbol(subject), Condition = !!symbol(condition),
      Metabolite = !!symbol(metabolite), value = !!symbol(value)
    )  %>%
    pivot_wider(names_from=Metabolite) %>%
    write.csv(tmp_out_file, row.names=FALSE)

}


#' Parse a file with metabolite quantities that is in MetaboAnalyst compatible format and return an initialized mSet object
#' 
#' @param analyst_compatible_data character. Path to the preformatted input file. Cannot be dataframe unfortunately.
#' @param analysis_type character. Code for the analysis we want to codoct on the mSet object later (e.g. "stat") 
#' @param input_format character. Input format of file to be parsed (almost always "rowu") 
#' 
#' @return Not intended to return anything, but rather to save outputs to files.
convert_cc_to_mSet <- function(analyst_compatible_data, analysis_type="stat", input_format="rowu"){
  
  mSet <- NULL # Have to reset otherwise the app resurrects old instance from global when chunk is rerun
  anal.type <- "stat"
  msg.vec <- list()

  InitDataObjects("conc", analysis_type) %>%
    Read.TextData(analyst_compatible_data, input_format) %>%
    SanityCheckData()

}

# Ensuring command line connectivity by sourcing an argument parser
source(opt$commandRpath, local=TRUE)