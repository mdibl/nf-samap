
process RUN_BLAST {
    tag "BLAST run"

    container 'staphb/blast:latest'

    input:
        path sample_sheet

    script:
    """
    @echo "Running BLAST with the provided input files..."
    """
}
