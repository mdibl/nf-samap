/*
 *  MODULE: load_sams.nf
 *
 *  Description: 
 *      Loads SAM objects using h5ad files described 
 *      in the sample sheet.
 *
 *  Inputs:
 *      run_id:         Timestamp of the nextflow process
 *      sample_sheet:   Path to the sample sheet containing the sample metadata
 *      data_dir:       Staging the data directory so the script can access it
 *
 *  Outputs:
 *      One pickled SAM object per sample and a logfile.
 *      results/run_id/sams/id2.pkl
 *      results/run_id/logs/run_id_load_sams.log
 */

process LOAD_SAMS {
    tag "${run_id} - load and pickle SAM objects"

    publishDir("results/${run_id}/sams/", mode: 'copy', pattern: '*.pkl')
    publishDir("results/${run_id}/logs/", mode: 'copy', pattern: '*.log')

    container 'pipeline/samap:latest'

    input:
        val run_id
        path sample_sheet 
        path data_dir 

    output:
        path "*.pkl", emit: sams
        path "${run_id}_load_sams.log", emit: logfile

    script:
    """
    LOG="${run_id}_load_sams.log"
    load_sams.py \\
        --sample-sheet ${sample_sheet} 2>&1 | tee -a \$LOG
    """
}