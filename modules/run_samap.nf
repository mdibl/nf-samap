process RUN_SAMAP {
    tag "SAMap run"
    publishDir "results/samap_objects", mode: 'copy'

    container 'ryansonder/samap:latest'

    input:
        path config
        path data_dir

    output:
        path "*.pkl"

    script:
    """
    run_samap.py --config ${config}
    """
}
