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
    Rscript /data/scratch/szabo/LigandScreen/metabolomics/pipeline/bin/descriptive_stats.r\
    --outFile $tag\
    -i $conc_tab

    """
} 


/*
        *Create hitlists from stats tables
*/
desc_stats
    .transpose()
    .map{it[0].join(',')}
    .collect()
    .set{hit_tabs}
Channel
    .fromPath(input_folder)
    .map{[it.baseName + '_ora', it]}
    .set{input_tabs}

process hitlistCreator {

    label 'manycpu'

    input:
        tuple tag, score_tab from desc_stats

    output:
        tuple ora_tag, ora_list into top_hits

    """
    Rscript /data/scratch/szabo/LigandScreen/metabolomics/pipeline/bin/hitlist_extractor.r\
    --outFile $tag\
    -i $score_tab

    """
} 


        /*
                *Calculate pathway overrepresentation
        */
        top_hits
            .collectFile(){ item -> [ "${item[0]}.txt", item[1] + '
' ]}
            .set{hits_ora}

        process pwORA {

            publishDir '.', saveAs: { it.contains('.tsv') || it.contains('.xlsx') ? "../tables/$it" : "../figures/$it" }, mode: 'copy'

            input:
                file 'hitlists.txt' from hits_ora

            output:
                tuple "ora.tsv", "ora.pdf" into ora_stats

            """
            Rscript /data/scratch/szabo/LigandScreen/metabolomics/pipeline/bin/pathway_ora.r\
            -i  hitlists.txt

            """
        } 

