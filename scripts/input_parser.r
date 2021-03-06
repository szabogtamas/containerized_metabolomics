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
    default="tmp/tmp.csv",
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
    default="/home/rstudio/git_repo/scripts/commandR.r",
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

for (rn in names(scriptOptionalArgs)){
  if(!(rn %in% names(opt))){
    opt[[rn]] <- scriptOptionalArgs[[rn]][["default"]]
  }
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
  mSet <- convert_cc_to_mSet(results, opt$tmpLocation, opt$analysis_type, opt$input_format)

  invisible(mSet)
}


#' Convert primary input to a standardized, tidy format
#' 
#' @param input_data dataframe. The primary input
#' 
#' @return Not intended to return anything, but rather to save outputs to files.
standardize_metabo_data <- function(input_data){

  if(is.character(input_data) & length(input_data) == 1){
    warning(paste("Parsing data in file:", input_data, sep="\n"))
    input_data <- scriptMandatoryArgs %>%
      parser4tsv(., "inFile", .["inFile"], input_data) %>%
      .$inFile
  }
  
  input_data %>%
    rename(Metabolite = X, MetaboGroup = X.1) %>%
    pivot_longer(
      c(-Metabolite, -MetaboGroup), names_to="Condition", names_repair="universal"
    ) %>%
    mutate(
      Replicate = str_extract(Condition, "\\.+\\d+") %>%
        str_replace("\\.+", "") %>%
        as.numeric() %>%
        {ifelse(is.na(.), 1, . +1)} %>%
        paste("BATCH", ., sep="_"),
      Condition = str_replace(Condition, "\\.+\\d+", ""),
      Subject = paste(Condition, Replicate, sep="_")
    )
    
}

#' Write metabo dataframe to tmp file and read as mSet; if mor control is needed, use populate_mSet()
#' 
#' @param input_data dataframe. The primary input
#' @param tmp_location character. Path to a temporary file needed for mSet creation
#' @param analysis_type character. Code for the analysis we want to codoct on the mSet object later (e.g. "stat") 
#' @param input_format character. Input format of file to be parsed (almost always "rowu") 
#' 
#' @return Populated mSet object.
convert_cc_to_mSet <- function(input_data, tmpLocation="tmp/tmp.csv", analysis_type="stat", input_format="rowu"){

  tmp_wd <- getwd()
  tmp_dir <- dirname(tmpLocation)
  preexisted_dir <- dir.exists(tmp_dir)
  tmpLocation <- basename(tmpLocation)
  
  if(!preexisted_dir) dir.create(tmp_dir)
  setwd(dirname(tmpLocation))

  write_metabodf_tmp(input_data, tmp_out_file=tmpLocation)
  mSet <- populate_mSet(tmpLocation, analysis_type, input_format)
  
  setwd(tmp_wd)
  if(preexisted_dir){
    unlink(tmpLocation)
  } else {
    unlink(tmp_dir)
  }
    
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
      !!sym(subject), !!sym(condition), !!sym(metabolite), !!sym(value)
    )  %>%
    rename(
      Subject = !!sym(subject), Condition = !!sym(condition),
      Metabolite = !!sym(metabolite), value = !!sym(value)
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
populate_mSet <- function(analyst_compatible_data, analysis_type="stat", input_format="rowu"){
  
  mSet <- NULL # Have to reset otherwise the app resurrects old instance from global when chunk is rerun
  anal.type <- "stat"
  msg.vec <- list()

  InitDataObjects("conc", analysis_type) %>%
    Read.TextData(analyst_compatible_data, input_format) %>%
    SanityCheckData()

}

# Ensuring command line connectivity by sourcing an argument parser
source(opt$commandRpath, local=TRUE, chdir=TRUE)