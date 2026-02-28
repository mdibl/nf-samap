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
        tuple val(id), val(anno) 
        path Grouping_Analysis

    output:
        path "${run_id}_nonpairwise.log"

    script:
    """
    LOG="${run_id}_nonpairwise.log"
    nonpairwise_analysis.py --id ${id.join(' ')} \\
        --anno ${anno.join(' ')} \\
        --grp ${Grouping_Analysis} 2>&1 | tee -a \$LOG
    """
}
