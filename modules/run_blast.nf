#!/usr/env/bin nextflow

process run_blast {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(fasta_file), path(db_file)

    output:
    path "blast_results/${sample_id}.out"

    script:
    """
    mkdir -p blast_results
    blastn -query ${fasta_file} -db ${db_file} -out blast_results/${sample_id}.out -outfmt 6
    """
}   