process LOAD_SAMS {
    tag "Load SAM objects"
    publishDir "results/sams", mode: 'copy'

    container 'ryansonder/samap:latest'

    input:
        path data_dir // Directory containing the h5ad files
        path sample_sheet // Sample sheet with metadata about the samples

    output:
        path "*.pkl" // Processed h5ad files

    script:
    """
    # Run the h5ad files through SAM to generate pickled SAM objects
    load_sams.py \\
        --sample-sheet ${sample_sheet}
    """
}