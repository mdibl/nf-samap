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
 *      A pickled SAMAP object and a logfile
 *      results/run_id/samap_objects/samap_results.pkl
 *      results/run_dir/logs/run_id_run_samap.log
 */

process RUN_SAMAP {
    tag "${run_id} - run SAMap"

    publishDir("results/${run_id}/samap_objects/", mode: 'copy', pattern: '*.pkl')
    publishDir("results/${run_id}/logs", mode: 'copy', pattern: '*.log')

    container 'mdiblbiocore/samap:latest'

    input:
        val run_id
        path results_dir 
        path samap_object

    output:
        path "samap_results.pkl", emit: results
        path "${run_id}_run_samap.log", emit: logfile

    script:
    """
    LOG=${run_id}_run_samap.log
    run_samap.py -i ${samap_object} 2>&1 | tee -a \$LOG
    """
}
