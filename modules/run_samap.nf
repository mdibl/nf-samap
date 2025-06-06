#!/usr/env/bin nextflow

process run_samap {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(fasta_file), path(blast_results)

    output:
    path "samap_results/${sample_id}.txt"

    script:
    """
    mkdir -p samap_results
    samap -i ${fasta_file} -b ${blast_results} -o samap_results/${sample_id}.txt
    """
}