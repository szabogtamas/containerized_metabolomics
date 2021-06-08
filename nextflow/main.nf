#!/usr/bin/env nextflow

import java.text.SimpleDateFormat
def date = new Date()
def sdf = new SimpleDateFormat("dd/MM/yyyy")

/*
 *    Calculate basic descriptive statistics for experiment
 */
Channel
    .fromPath(params.input_folder + '/*.*')
    .map{[params.ttest_tag, it.baseName, it]}
    .set{ input_stats_tabs }

process runDescStats {

    publishDir '.', saveAs: { it.contains('.tsv') || it.contains('.xlsx') ? params.table_folder + "/$it" : params.figure_folder + "/$it" }, mode: 'copy'
    containerOptions '--bind /data:/data'

    input:
        tuple testtag, tag, conc_tab from input_stats_tabs

    output:
        tuple testtag, tag, "${tag}.tsv", "${tag}.pdf" into desc_stats

    """
    Rscript /home/rstudio/repo_files/scripts/descriptive_stats.r\
    --outFile $tag\
    -i $conc_tab

    """
}


/*
 *    Rank metabolites based on multivariate ROC performance
 */
Channel
    .fromPath(params.input_folder + '/*.*')
    .map{[params.mroc_tag, it.baseName, it]}
    .set{ input_mroc_tabs }

process runMultiROC {

    publishDir '.', saveAs: { it.contains('.tsv') || it.contains('.xlsx') ? params.table_folder + "/$it" : params.figure_folder + "/$it" }, mode: 'copy'
    containerOptions '--bind /data:/data'

    input:
        tuple testtag, tag, conc_tab from input_mroc_tabs

    output:
        tuple testtag, tag, "${tag}.tsv", "${tag}.pdf" into mroc_stats

    """
    Rscript /home/rstudio/repo_files/scripts/multivariate_roc.r\
    --outFile $tag\
    -i $conc_tab

    """
}


/*
 *    Create hitlists from stats tables
 */
desc_stats
    .mix(mroc_stats)
    .set{ stat_tabs }

process hitlistCreator {
    
    containerOptions '--bind /data:/data'

    input:
        tuple testtag, tag, score_tab, score_plot from stat_tabs

    output:
        tuple testtag, "hits_${tag}.txt" into top_hits

    """
    Rscript /home/rstudio/repo_files/scripts/hitlist_extractor.r\
    --outFile hits_$tag\
    --tabLabels $tag\
    -i $score_tab

    """
} 


/*
 *    Calculate pathway overrepresentation
 */

top_hits
    .collectFile(){ item -> [ "${item[0]}.txt", item[1] + '' ]}
    .map{[params.ora_tag, it.baseName, '"' + it.text.replaceAll("\\n", "::") + '"']}
    .set{ hits_ora }

process pwORA {
    
    containerOptions '--bind /data:/data'
    publishDir '.', saveAs: { it.contains('.tsv') || it.contains('.xlsx') ? params.table_folder + "/$it" : params.figure_folder + "/$it" }, mode: 'copy'
    
    input:
        tuple tag, testtag, score_tab from hits_ora

    output:
        tuple "${testtag}_ora.tsv", "${testtag}_ora.pdf" into ora_stats

    """
    Rscript /home/rstudio/repo_files/scripts/pathway_ora.r\
    -i  $score_tab\
    --outFile ${testtag}_${tag}
    """
} 


/*
*    Run MSEA to detect pathway enrichment
*/

Channel
    .fromPath(params.input_folder + '/*.*')
    .map{[params.msea_tag, it.baseName, it]}
    .set{ input_for_msea }

process pwMSEA {
    
    containerOptions '--bind /data:/data'
    publishDir '.', saveAs: { it.contains('.tsv') || it.contains('.xlsx') ? params.table_folder + "/$it" : params.figure_folder + "/$it" }, mode: 'copy'
    
    input:
        tuple testtag, tag, in_tab from input_for_msea

    output:
        tuple "${testtag}_${tag}.tsv", "${testtag}_${tag}.pdf" into msea_stats

    """
    Rscript /home/rstudio/repo_files/scripts/pathway_msea.r\
    -i  $in_tab\
    --outFile ${testtag}_${tag}
    """
} 