#!/usr/bin/env Rscript

scriptDescription <- "A script that reads a table from Google Drive."

scriptMandatoryArgs <- list(
  drive_path = list(
    abbr="-i",
    readoptions=list(sep="\t", stringsAsFactors=FALSE),
    help="Data table on Drive."
  )
)

scriptOptionalArgs <- list(
  outFile = list(
    default=NULL,
    help="Output location if we insist on saving the table locally."
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

for (pk in c("tidyr", "dplyr", "readxl", "googledrive")){
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
  
  outFile <- opt$outFile
  opt$outFile <- NULL
  opt$help <- NULL
  opt$verbose <- NULL

  cat("Reading from Drive\n")
  data <- do.call(read_drive, opt)
  
  if(!is.null(outFile)){
    cat("Saving table\n")
    tab2tsv(data, outFile)
  }

  invisible(NULL)
}


#' Downloads a data table from Drive and reads it as a dataframe.
#' 
#' @param drive_path string. Pseudopath to data table on Drive.
#' 
#' @return dataframe or list of them, if multisheet.
read_drive <- function(drive_path){
  
  drive_download(
    drive_path,
    overwrite = TRUE
  )

  data_file <- basename(drive_path)
  
  if(tail(unlist(strsplit(data_table, "\\.")), -1) %in% c("xls", "xlsx")){
    data_table <- readxl(data_file)
  } else {
    data_table <- read.csv(data_file, sep = "\t", stringsAsFactors = FALSE)
  }
  
  unlink(data_file)

  return(data_table)
  
}

# Ensuring command line connectivity by sourcing an argument parser
source(opt$commandRpath, local=TRUE)