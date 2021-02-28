#!/usr/bin/env Rscript

scriptDescription <- "A script that parses metabolomics measurement results."

scriptMandatoryArgs <- list(
  inFile = list(
    abbr="-i",
    type="table",
    readoptions=list(sep="\t", stringsAsFactors=FALSE),
    help="Table of metobolite quantities/concentrations."
  )
)

scriptOptionalArgs <- list(
  conditionOrder = list(
    default=NULL,
    type="vector",
    help="Order of conditions in the experimental design formula. Makes sense to put control as first."
  ),
  tmpLocation = list(
    default="tmp.csv",
    help="Path to the temporary file used to push results in an mSet object."
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

for (pk in c("tidyr", "dplyr", "MetaboAnalystR")){
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
  
  opt$help <- NULL
  opt$verbose <- NULL

  mSet <- NULL # Have to reset otherwise the app resurrects old instance from global when chunk is rerun
  anal.type <- "stat"
  msg.vec <- list()

  cat("Standardizing input\n")
  results <- standardize_metabo_data(opt$inFile, subject, condition, metabolite, value, tmp_out_file)
  
  cat("Saving to temporary MetoboAnalyst file\n")
  tmp_wd <- getwd()
  setwd(dirname(opt$tmpLocation))
  write_metabodf_tmp(results, subject, condition, metabolite, value, opt$tmpLocation)
  mSet <- convert_cc_to_mSet(opt$tmpLocation, analysis_type, input_format)
  unlink(tmpLocation)
  setwd(tmp_wd)
  
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
    pivot_longer(c(-Metabolite, -MetaboGroup)) %>%
    mutate(
      Condition = str_extract(name, 'Condition(A|B)'),
      Replicate = case_when(
        grepl('\\.1', name) ~ 'Batch_2',
        grepl('\\.2', name) ~ 'Batch_3',
        TRUE ~ 'Batch_1'
      ),
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
write_metabodf_tmp <- function(metab_data, subject, condition, metabolite, value, tmp_out_file){
  
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
convert_cc_to_mSet <- function(analyst_compatible_data, analysis_type, input_format){
  
  InitDataObjects("conc", analysis_type) %>%
    Read.TextData(analyst_compatible_data, input_format) %>%
    SanityCheckData()

}

# Ensuring command line connectivity by sourcing an argument parser
source(opt$commandRpath, local=TRUE)