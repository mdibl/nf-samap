process RUN_SAMAP {
    tag "SAMap run"
    publishDir "results/samap_objects", mode: 'copy'

    container 'ryansonder/samap:latest'

    input:
        path results_dir
        path data_dir
        path sample_sheet

    output:
        path "*.pkl"

    script:
    """
    run_samap.py --sample-sheet ${sample_sheet} --maps ${results_dir}/maps
    """
}
