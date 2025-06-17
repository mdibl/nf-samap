/*
 *  MODULE: run_samap.nf
 *
 *  Description: 
 *      Runs the SAMap algorithm on a SAMAP object
 *
 *  Inputs:
 *      run_id:         Timestamp of the nextflow process
 *      results_dir:    Directory to store the results
 *      samap_object:   Channel containing a pickled SAMAP object
 *
 *  Outputs:
 *      A pickled SAMAP object
 *      results/${run_id}/samap_objects/samap_results.pkl
 */

process RUN_SAMAP {
    tag "${run_id} - run SAMap"

    publishDir "results/${run_id}/samap_objects/", mode: 'copy'

    container 'pipeline/samap:latest'

    input:
        val run_id
        path results_dir 
        path samap_object

    output:
        path "samap_results.pkl"

    script:
    """
    run_samap.py -i ${samap_object}
    """
}
