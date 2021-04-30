#!/usr/bin/env Rscript

if (exists("block_hitlist")) eval_blocker <- TRUE else eval_blocker <- NULL

scriptDescription <- "Extract top N hits for ranked table and format hitlist."

scriptMandatoryArgs <- list(
  changeValues = list(
    abbr="-i",
    type="tables",
    readoptions=list(stringsAsFactors=FALSE),
    help="Metabolite names and a change measure or score."
  )
)

scriptOptionalArgs <- list(
  nHit = list(
    default=25,
    help="Number of top hits to be extracted."
  ),
  metabCol = list(
    default="Metabolite",
    help="Name of column to be extracted."
  ),
  scoreCol = list(
    default="p.value",
    help="Name of column with stats result (inverse score)."
  ),
  scoreDesc = list(
    default=FALSE,
    help="If descending order of score should be taken."
  ),
  tabLabels = list(
    default=NULL,
    type="vector",
    help="Labels to be associated with input tables."
  ),
  outFile = list(
    default="hitlist",
    help="File path without extension to hitlist."
  ),
  commandRpath = list(
    default="/home/rstudio/repo_files/scripts/commandR.r",
    help="Path to command line connectivity script (if not in cwd)."
  )
)

for (pk in c("tidyr", "dplyr", "purrr")){
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

  if(!is.null(opt$tabLabels)){
    names(opt$changeValues) <- opt$tabLabels
  }

  cat("Generating hitlists\n")
  opt %>%
    do.call(generate_std_hitlist, .) %>%
    {paste(names(.), map(., paste, collapse=","), sep=",")} %>%
    writeLines(paste(opt$outFile, "txt", sep="."))
  
  invisible(NULL)
}


#' Generate hitlsts from score tables that can be parsed by ORA or Venn-like plots
#' 
#' @param changeValues list. Metabolite score tables.
#' @param nHit integer. Number of top hits to be extracted.
#' @param metabCol string. Name of column to be extracted.
#' @param scoreCol string. Name of column with stats result (inverse score).
#' @param scoreDesc string. If descending order of score should be taken.
#' 
#' @return list  List of hitslists for each condition.
generate_std_hitlist <- function(changeValues, nHit=25, metabCol="Metabolite", scoreCol="p.value", scoreDesc=FALSE, ...){
  
  changeValues %>%
    map(~filter.x, !is.na(!!sym(scoreCol))) %>%
    map(
      function(x) {
        if(scoreDesc){
          arrange(x, desc(!!sym(scoreCol)))
        } else {
          arrange(x, !!sym(scoreCol))
        }
      }
    ) %>%
    map(head, nHit) %>%
    map(function(x) x[[metabCol]])
  
}

# Ensuring command line connectivity by sourcing an argument parser
rg <- commandArgs()
if("--commandRpath" %in% rg){
  scriptOptionalArgs$commandRpath$default <- rg[[which(rg == "--commandRpath") + 1]]
}
source(scriptOptionalArgs$commandRpath$default, local=TRUE, chdir=FALSE)