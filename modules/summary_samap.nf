/*
 *  MODULE: summary_samap.nf
 *
 *  Description: 
 *      Produces some heatmap visualizations and performs basic pairwise cell-cell summary statistics regarding SAMap output
 *     
 *
 *  Inputs:
 *      samap_results:     Timestamp of the nextflow process
 *      idCompare:         Channel containing a pairwise species-species comparison to analyze
 *      outdir:            Directory in which to save final output
 *
 *  Outputs:
 *      Several visualizations about the SAMap results, an enhanced pairwise-mapping score file, and a logfile
 */

process SUMMARY_SAMAP {
    tag "${run_id} - SAMap Top-Level Post-Analysis"

    container 'mdiblbiocore/postanalysis:latest'


    input:
        val run_id
        path samap_obj
        tuple val(id1), val(id2), val(anno1), val(anno2)
        tuple val(id), val(anno)

    output:
        path "*.png"
        path "*.csv"
        tuple val(id1), val(id2), path("GenePairs.csv"),                        emit: genepairs
        tuple val(id1), val(id2), path("pms_cluster_alignment_scores.csv"),     emit: pms
        tuple val(id1), val(id2), path("samap_results_Cleaned.pkl"),            emit: samap_cleaned
        path "${run_id}_summary.log"

    script:
    """
    LOG="${run_id}_summary.log"
    summary_samap.py --input ${samap_obj} --id1 ${id1} --anno1 ${anno1} --id2 ${id2} --anno2 ${anno2} --allId ${id.join(' ')} --allAnno ${anno.join(' ')} 2>&1 | tee -a \$LOG
    """
}
