process RUN_SAMAP {
    tag "SAMap run"
    publishDir "results/samap_objects", mode: 'copy'

    container 'ryansonder/samap:latest'

    input:
        path results_dir // Directory to store the results
        path samap_object // SAMap object to be used in the run

    output:
        path "samap_results.pkl"

    script:
    """
    run_samap.py -i ${samap_object}
    """
}
