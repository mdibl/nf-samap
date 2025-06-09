

process RUN_SAMAP {
    tag "SAMap run"

    container 'ryansonder/samap:latest'

    input:
        path config

    output:
        path "samap_obj.pkl"

    script:
    """
    run_samap.py --config ${config}
    """
}
