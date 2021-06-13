#!/usr/bin/env Rscript

if (exists("block_knit")) eval_blocker <- TRUE else eval_blocker <- NULL

scriptDescription <- "A wrapper around knitr that helps joining pdf figures into a report from command line."

scriptMandatoryArgs <- list(
  reportTemplate = list(
    help="Path to the template Rmd file."
  )
)

scriptOptionalArgs <- list(
  outFile = list(
    default="metabolomics_report.pdf",
    help="File path to report output."
  ),
  commandRpath = list(
    default="/home/rstudio/repo_files/scripts/commandR.r",
    help="Path to command line connectivity script (if not in cwd)."
  )
)

#' The main function of the script, executed only if called from command line.
#' Calls subfunctions according to supplied command line arguments.
#' 
#' @param opt list. a named list of all command line options; will be passed on 
#' 
#' @return Not intended to return anything, but rather to save outputs to files.
main <- function(opt){
    
  rmarkdown::render(
    opt$reportTemplate,
    output_file = opt$outFile,
    params = list(figures=opt$positionals)
  )
  invisible(NULL)
}


# Ensuring command line connectivity by sourcing an argument parser
rg <- commandArgs()
if("--commandRpath" %in% rg){
  scriptOptionalArgs$commandRpath$default <- rg[[which(rg == "--commandRpath") + 1]]
}
source(scriptOptionalArgs$commandRpath$default, local=TRUE, chdir=FALSE)