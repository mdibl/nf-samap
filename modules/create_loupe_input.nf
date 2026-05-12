/*
 *  MODULE: create_loupe_input.nf
 *
 *  Description: 
 *      Loads SAM objects using h5ad files described 
 *      in the sample sheet.
 *
 *  Inputs:
 *      run_id:         Timestamp of the nextflow process
 *      pairwise:       Tuple of [id1, id2, anno1, anno2, h5ad1, h5ad2]
 *      samap_output:   Path to the SAMap pickle file
 *
 *  Outputs:
 *      Combined h5ad file for Loupe input and a logfile.
 */

process CREATE_LOUPE_INPUT {
    tag "${run_id} - merge raw counts matrices for pairwise comparison ${id1}-${id2}"

    container 'mdiblbiocore/postanalysis:latest'

    input:
        val run_id
        tuple val(id1), val(id2), val(anno1), val(anno2), path(h5ad1), path(h5ad2)
        path samap_output

    output:
        tuple val(id1), val(id2),
            path("${id1}_${id2}_counts.mtx"),
            path("${id1}_${id2}_barcodes.csv"),
            path("${id1}_${id2}_features.csv"),
            path("${id1}_${id2}_umap.csv"),
            path("${id1}_${id2}_meta.csv"),
            emit: loupeinput
        path "${run_id}_${id1}_${id2}_create_loupe_input.log", emit: logfile

    script:
    """
    LOG="${run_id}_${id1}_${id2}_create_loupe_input.log"
    create_loupe_input.py \\
        --samap         ${samap_output} \\
        --h5ads         ${h5ad1} ${h5ad2} \\
        --species       ${id1} ${id2} \\
        --celltype-cols ${anno1} ${anno2} \\
        2>&1 | tee -a \$LOG
    """
}
