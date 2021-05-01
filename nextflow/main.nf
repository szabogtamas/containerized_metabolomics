#!/usr/bin/env nextflow

import java.text.SimpleDateFormat
def date = new Date()
def sdf = new SimpleDateFormat("dd/MM/yyyy")

/*
        *Calculate basic descriptive statistics for experiment
*/
Channel
    .fromPath(input_folder)
    .map{[it.baseName + '_stats', it]}
    .set{input_tabs}

process runDescStats {

    publishDir '.', saveAs: { it.contains('.tsv') || it.contains('.xlsx') ? "../tables/$it" : "../figures/$it" }, mode: 'copy'

    input:
        tuple tag, conc_tab from input_tabs

    output:
        tuple tag, "${tag}.tsv", "${tag}.pdf" into desc_stats

    """
    Rscript repo_files/scripts/descriptive_stats.r\
    --outFile $tag\
    -i $conc_tab

    """
}