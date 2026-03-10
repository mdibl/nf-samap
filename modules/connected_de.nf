/*
 *  MODULE: connected_de.nf
 *
 *  Description: 
 *      Runs Comprehensive Differential Expression analysis on a Pairwise species-species comparison
 *     
 *
 *  Inputs:
 *      samap_results:     Cleaned samap object from the previous module
 *      idCompare:         Channel containing a pairwise species-species comparison to analyze
 *      genepairs:         Table of Pairiwse Gene Pairs calculated in the previous Module
 *      outdir:            Directory in which to save final output
 *
 *  Outputs:
 *      One nested directory of DE output per Alignment Family and a logfile
 */

process CONNECTED_DE {
    cache = false
    tag "${run_id} - SAMap Differential Expression Post-Analysis"

    container 'mdiblbiocore/postanalysis:latest'


    input:
        val run_id
        tuple val(id1), val(id2), val(anno1), val(anno2), path(pms), path(samap_obj)

    output:
        tuple val(id1), val(id2), path("Grouping_Analysis/"),   emit: groupinganalysis
        tuple val(id1), val(id2), path("analysis.pkl"),          emit: analysis
        tuple val(id1), val(id2), path("all_de_results.pkl"),    emit: all_de_results
        path "${run_id}_connectedDE.log"

    script:
    """
    LOG="${run_id}_connectedDE.log"
    connected_de.py --input ${samap_obj} --pms ${pms} --id1 ${id1} --anno1 ${anno1} --id2 ${id2} --anno2 ${anno2} 2>&1 | tee -a \$LOG
    """
}
