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
 *      One pickled SAM object per sample
 *      results/${run_id}/sams/{id2}.pkl
 */

process LOAD_SAMS {
    tag "${run_id} - load and pickle SAM objects"

    publishDir "results/${run_id}/sams/", mode: 'copy'

    container 'ryansonder/samap:latest'

    input:
        val run_id
        path sample_sheet 
        path data_dir 

    output:
        path "*.pkl"

    script:
    """
    load_sams.py \\
        --sample-sheet ${sample_sheet}
    """
}