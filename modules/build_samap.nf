/*
 *  MODULE: build_samap.nf
 *
 *  Description: 
 *      Creates a SAMAP object from blast mappings and SAM objects
 *
 *  Inputs:
 *      run_id:         Timestamp of the nextflow process
 *      sample_sheet:   Path to the sample sheet containing the sample metadata
 *      data_dir:       Staging the data directory so the script can access it
 *      maps_dir:       Directory containing the BLAST mappings
 *      results_dir:    Directory containing the SAM objects
 *      sams:           Channel containing the SAM objects
 *
 *  Outputs:
 *      A pickled SAMAP object and a logfile.
 *      results/run_id/samap_objects/run_id_samap.pkl
 *      results/run_id/logs/run_id_build_samap.log
 */

process BUILD_SAMAP {
    tag "${run_id} - build SAMAP object"

    publishDir("results/${run_id}/samap_objects/", mode: 'copy', pattern: '*.pkl')
    publishDir("results/${run_id}/logs/", mode: 'copy', pattern: '*.log')

    container 'pipeline/samap:latest'

    input:
        val run_id
        path sample_sheet // Sample sheet containing metadata for the samples
        path data_dir // Load all the data necessary for building the SAMap object
        path maps_dir // Directory containing the BLAST mappings
        path results_dir // Directory to store the results
        path sams // SAM objects to be used in the SAMap

    output:
        path "samap.pkl", emit: samap
        path "${run_id}_build_samap.log", emit: logfile

    script:
    """
    LOG="${run_id}_build_samap.log"
    build_samap.py \\
        --sams-dir ${results_dir}/${run_id}/sams \\
        --sample-sheet ${sample_sheet} \\
        --maps ${maps_dir} 2>&1 | tee -a \$LOG
    """
}