#!/usr/env/bin nextflow

process visualize_results {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(samap_results)

    output:
    path "visualization/${sample_id}_results.png"

    script:
    """
    mkdir -p visualization
    python visualize.py --input ${samap_results} --output visualization/${sample_id}_results.png
    """
}