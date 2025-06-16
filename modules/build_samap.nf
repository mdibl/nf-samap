process BUILD_SAMAP {
    tag "Build SAMAP object"
    publishDir "results/samap_objects", mode: 'copy'

    container 'ryansonder/samap:latest'

    input:
        path data_dir // Load all the data necessary for building the SAMap object
        path sams // SAM objects to be used in the SAMap
        path maps_dir // Directory containing the BLAST mappings
        path results_dir // Directory to store the results
        path sample_sheet // Sample sheet containing metadata for the samples

    output:
        path "*.pkl"

    script:
    """
    build_samap.py --sams-dir ${results_dir}/sams --sample-sheet ${sample_sheet} --f-maps ${maps_dir}
    """
}