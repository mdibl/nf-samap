/*
 *  MODULE: nonpairwise_analysis.nf
 *
 *  Description: 
 *      Consolidates upstream pairwise analysis to draw potentially >2 species analyses.
 *     
 *
 *  Inputs:
 *      
 *
 *  Outputs:
 *
 */

process NONPAIRWISE_ANALYSIS {
    tag "${run_id} - SAMap Pairwise Analysis Consolidation"

    container 'mdiblbiocore/postanalysis:latest'


    input:
        path Grouping_Analysis

    output:
        path "${run_id}_nonpairwise.log"

    script:
    """
    LOG="${run_id}_nonpairwise.log"
    nonpairwise_analysis.py --input ${samap_obj} --id1 ${id1} --anno1 ${anno1} --id2 ${id2} --anno2 ${anno2} 2>&1 | tee -a \$LOG
    """
}
