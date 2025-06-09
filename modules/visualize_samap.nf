process VISUALIZE_SAMAP {
    tag "SAMap visualization"

    container 'ryansonder/samap:latest'

    // input:
    //     path samap_obj

    // output:

    script:
    """
    visualize_samap.py
    """
}
