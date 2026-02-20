/*
 *  MODULE: connected_de.nf
 *
 *  Description: 
 *      Runs Comprehensive Differential Expression analysis on a Pairwise species-species comparison
 *     
 *
 *  Inputs:
 *      samap_results:     Timestamp of the nextflow process
 *      idCompare:         Channel containing a pairwise species-species comparison to analyze
 *      genepairs:         Table of Pairiwse Gene Pairs calculated in the previous Module
 *      outdir:            Directory in which to save final output
 *
 *  Outputs:
 *      One nested directory of DE output per Alignment Family and a logfile
 */

process CONNECTED_DE {
    tag "${run_id} - SAMap Differential Expression Post-Analysis"

    container 'mdiblbiocore/postanalysis:latest'

    input:
        val run_id
        path samap_obj
        path pms
        tuple val(id1), val(id2), val(anno1), val(anno2)

    output:
        path "${id1}-${id2}/Grouping_Analysis/", emit: de_results_dir
        path "${run_id}_summary.log"

    script:
    """
    LOG="${run_id}_connectedDE.log"
    connected_de.py --input ${samap_obj} --pms ${pms} --id1 ${id1} --anno1 ${anno1} --id2 ${id2} --anno2 ${anno2} 2>&1 | tee -a \$LOG


    mkdir ${id1}-${id2}
    mv Grouping_Analysis/ ${id1}-${id2}
    """
}
