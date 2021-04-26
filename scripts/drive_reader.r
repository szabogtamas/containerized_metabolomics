#!/usr/bin/env Rscript

if (exists("block_drive")) eval_blocker <- TRUE else eval_blocker <- NULL

scriptDescription <- "A script that reads a table from Google Drive."

scriptMandatoryArgs <- list(
  drivePath = list(
    abbr="-i",
    help="Data table on Drive."
  ),
  outFile = list(
    abbr="-o",
    help="Output location for saving the table locally."
  )
)

scriptOptionalArgs <- list(
  tokenPath = list(
    default="~/local_files/.secrets",
    help="Path to folder where oauth tokens are saved."
  ),
  commandRpath = list(
    default="/home/rstudio/repo_files/scripts/commandR.r",
    help="Path to command line connectivity script (if not in cwd)."
  )
)

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
  
  cat("Authenticating Drive token\n")
  drive_auth(cache = opt$tokenPath, email = TRUE)
  
  cat("Reading from Drive\n")
  data <- read_drive(opt$drivePath)
  
  cat("Saving table\n")
  tab2tsv(data, opt$outFile)

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
  
  if(tail(unlist(strsplit(data_file, "\\.")), -1) %in% c("xls", "xlsx")){
    data_table <- read_excel(data_file)
  } else {
    data_table <- read.csv(data_file, sep = "\t", stringsAsFactors = FALSE)
  }
  
  unlink(data_file)

  return(data_table)
  
}

# Ensuring command line connectivity by sourcing an argument parser
rg <- commandArgs()
if("--commandRpath" %in% rg){
  scriptOptionalArgs$commandRpath$default <- rg[[which(rg == "--commandRpath") + 1]]
}
source(scriptOptionalArgs$commandRpath$default, local=TRUE, chdir=FALSE)