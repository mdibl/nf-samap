
process RUN_BLAST {
    tag "BLAST run"

    container 'staphb/blast:latest'

    script:
    """
    @echo "Running BLAST with the provided input files..."
    """
}
