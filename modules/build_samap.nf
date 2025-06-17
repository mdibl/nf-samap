/*
 *  MODULE: build_samap.nf
 *
 *  Description: 
 *      Creates a SAMAP object from 
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
 *      A pickled SAMAP object
 *      results/${run_id}/samap_objects/${run_id}_samap.pkl
 */

process BUILD_SAMAP {
    tag "${run_id} - build SAMAP object"

    publishDir "results/${run_id}/samap_objects/", mode: 'copy'

    container 'pipeline/samap:latest'

    input:
        val run_id
        path sample_sheet // Sample sheet containing metadata for the samples
        path data_dir // Load all the data necessary for building the SAMap object
        path maps_dir // Directory containing the BLAST mappings
        path results_dir // Directory to store the results
        path sams // SAM objects to be used in the SAMap

    output:
        path "samap.pkl"

    script:
    """
    build_samap.py \\
        --sams-dir ${results_dir}/${run_id}/sams \\
        --sample-sheet ${sample_sheet} \\
        --f-maps ${maps_dir}
    """
}