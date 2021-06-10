#!/usr/bin/env nextflow

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Run basic metabolomics analysis based on MetaboAnalyst toolkit        //
//                                                                        //
////////////////////////////////////////////////////////////////////////////


/*
 *    Calculate basic descriptive statistics for experiment
 */

Channel
    .fromPath(params.input_folder + '/*.*')
    .map{[params.ttest_tag, it.baseName, it]}
    .set{ input_stats_tabs }

process runDescStats {

    publishDir '.', saveAs: { it.contains('.tsv') || it.contains('.xlsx') ? params.table_folder + "/$it" : params.figure_folder + "/$it" }, mode: 'copy'

    input:
        tuple testtag, tag, conc_tab from input_stats_tabs

    output:
        tuple testtag, tag, "${tag}.tsv" into desc_stats
        tuple "${testtag}_${tag}", "${tag}.pdf" into stats_figures

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

    input:
        tuple testtag, tag, conc_tab from input_mroc_tabs

    output:
        tuple testtag, tag, "${tag}.tsv" into mroc_stats
        tuple "${testtag}_${tag}", "${tag}.pdf" into mroc_figures

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
 *    Calculate pathway overrepresentation based on hitlists
 */

top_hits
    .collectFile(){ item -> [ "${item[0]}.txt", item[1] + '' ]}
    .map{[params.ora_tag, it.baseName, '"' + it.text.replaceAll("\\n", "::") + '"']}
    .set{ hits_ora }

process pwORA {
    
    publishDir '.', saveAs: { it.contains('.tsv') || it.contains('.xlsx') ? params.table_folder + "/$it" : params.figure_folder + "/$it" }, mode: 'copy'
    
    input:
        tuple tag, testtag, score_tab from hits_ora

    output:
        tuple "${testtag}_${tag}", "${testtag}_ora.tsv" into ora_stats
        tuple "${testtag}_${tag}", "${testtag}_ora.pdf" into ora_figures

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
    
    publishDir '.', saveAs: { it.contains('.tsv') || it.contains('.xlsx') ? params.table_folder + "/$it" : params.figure_folder + "/$it" }, mode: 'copy'
    
    input:
        tuple testtag, tag, in_tab from input_for_msea

    output:
        tuple "${testtag}_${tag}", "${testtag}_${tag}.tsv" into msea_stats
        tuple "${testtag}_${tag}", "${testtag}_${tag}.pdf" into msea_figures

    """
    Rscript /home/rstudio/repo_files/scripts/pathway_msea.r\
    -i  $in_tab\
    --outFile ${testtag}_${tag}
    """
} 


/*
 *    Find enriched KEGG pathways
 */

Channel
    .fromPath(params.input_folder + '/*.*')
    .map{[params.kegg_tag, it.baseName, it]}
    .set{ input_for_kegg }

process pwKEGG {
    
    publishDir '.', saveAs: { it.contains('.tsv') || it.contains('.xlsx') ? params.table_folder + "/$it" : params.figure_folder + "/$it" }, mode: 'copy'
    
    input:
        tuple testtag, tag, in_tab from input_for_kegg

    output:
        tuple "${testtag}_${tag}", "${testtag}_${tag}.tsv" into kegg_stats
        tuple "${testtag}_${tag}", "${testtag}_${tag}.pdf" into kegg_figures

    """
    Rscript /home/rstudio/repo_files/scripts/pathway_kegg.r\
    -i  $in_tab\
    --outFile ${testtag}_${tag}
    """
} 


/*
 *    Gather all figures into a single report
 */

process reportGenerator {
    
    publishDir params.report_folder, mode: 'copy'
    
    input:
        file stats_tag, stats_figure from stats_figures
        file mroc_tag, mroc_figure from mroc_figures
        file ora_tag, ora_figure from ora_figures
        file msea_tag, msea_figure from msea_figures
        file kegg_tag, kegg_figure from kegg_figures
        val report_filename from params.report_filename
        val report_title from params.report_title
        val report_author from params.report_author

    output:
        file $report_filename into final_report

    """
    R --slave -e "rmarkdown::render(
        '/home/rstudio/repo_files/streamlined_report_template.Rmd',
        output_file = ${report_filename},
        params = list(report_title = report_title, report_author = report_author)
    )"
    """
} 