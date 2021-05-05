#!/usr/bin/env nextflow

import java.text.SimpleDateFormat
def date = new Date()
def sdf = new SimpleDateFormat("dd/MM/yyyy")

/*
        *Calculate basic descriptive statistics for experiment
*/
Channel
    .fromPath(params.input_folder + '/*.*')
    .map{[params.ttest_tag, it.baseName + '_stats', it]}
    .set{input_tabs}

process runDescStats {

    publishDir '.', saveAs: { it.contains('.tsv') || it.contains('.xlsx') ? "../tables/$it" : "../figures/$it" }, mode: 'copy'
    containerOptions '--bind /data:/data'
    echo true

    input:
        tuple testtag, tag, conc_tab from input_tabs

    output:
        tuple testtag, tag, "${tag}.tsv", "${tag}.pdf" into desc_stats

    """
    Rscript /home/rstudio/repo_files/scripts/descriptive_stats.r\
    --outFile $tag\
    -i $conc_tab

    """
}


/*
        *Create hitlists from stats tables
*/

process hitlistCreator {
    
    containerOptions '--bind /data:/data'
    echo true

    input:
        tuple testtag, tag, score_tab, score_plot from desc_stats

    output:
        tuple testtag, "hits_${tag}.txt" into top_hits

    """
    Rscript /home/rstudio/repo_files/scripts/hitlist_extractor.r\
    --outFile hits_$tag\
    -i $score_tab

    """
} 


/*
        *Calculate pathway overrepresentation
*/

top_hits
    .collectFile(){ item -> [ "${item[0]}.txt", item[1] + '' ]}
    .map{[it.baseName, it]}
    .set{hits_ora}

process pwORA {
    
    containerOptions '--bind /data:/data'
    echo true
    publishDir '.', saveAs: { it.contains('.tsv') || it.contains('.xlsx') ? "../tables/$it" : "../figures/$it" }, mode: 'copy'

    input:
        tuple testtag, score_tab from hits_ora

    output:
        tuple "${testtag}_ora.tsv", "${testtag}_ora.pdf" into ora_stats

    """
    Rscript /home/rstudio/repo_files/scripts/pathway_ora.r\
    -i  $score_tab
    """
} 