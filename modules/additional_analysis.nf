/*
 *  MODULE: additional_analysis.nf
 *
 *  Description: 
 *      Runs some auxiliary gene pair differential expression combination functions to compare across both gene pair and DE metrics. Output is stored in same Grouping_Analysis folder as above
 *     
 *
 *  Inputs:
 *      idCompare:         Channel containing a pairwise species-species comparison to analyze
 *      genepairs:         Table of Pairiwse Gene Pairs calculated in the previous Module
 *      de_dir             Directory containing the Differentially expression output tables
 *      outdir:            Directory in which to save final output
 *
 *  Outputs:
 *      One nested directory of DE / GenePair output per Alignment Family and a logfile
 */

process ADDITIONAL_ANALYSIS {
    cache = false
    tag "${run_id} - SAMap GenePair-DE analysis"

    container 'mdiblbiocore/postanalysis:latest'

    publishDir "${params.outdir}/Analysis/GroupingAnalysis/${id1}-${id2}", mode: params.publish_dir_mode ?: 'copy'

    input:
        val run_id
        tuple val(id1), val(id2),
          path(groupinganalysis),
          path(analysis),
          path(all_de_results),
          path(genepairs),
          path(pms),
          val(anno1),
          val(anno2)

    output:
        path "Grouping_Analysis/",          emit: de_results_complete
        path "${run_id}_additionalAnalysis.log"

    script:
    """
    LOG="${run_id}_additionalAnalysis.log"
    additional_analysis.py --genepairs ${genepairs} --diff ${groupinganalysis} --pms ${pms} --id1 ${id1} --anno1 ${anno1} --id2 ${id2} --anno2 ${anno2} --analysis ${analysis} --all_de_results ${all_de_results} 2>&1 | tee -a \$LOG
    """
}
