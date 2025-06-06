#!/usr/env/bin nextflow

process build_blast_db {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(fasta_file)

    output:
    path "blast_db/${sample_id}.*"

    script:
    """
    mkdir -p blast_db
    makeblastdb -in ${fasta_file} -dbtype nucl -out blast_db/${sample_id}
    """
}