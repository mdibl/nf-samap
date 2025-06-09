

process RUN_SAMAP {
    tag "SAMap run"

    container 'avianalter/samap:latest'

    input:
        path json_path
        path data_dir

    output:
        path "samap_obj.pkl"


    script:
    """
    run_samap.py --json-path ${json_path}
    """
}
