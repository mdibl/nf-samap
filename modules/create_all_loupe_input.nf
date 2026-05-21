/*
 *  MODULE: create_all_loupe_input.nf
 *
 *  Description:
 *      Merges raw counts from all N species and extracts the full-SAMap UMAP
 *      embedding to produce the inputs needed for all-species Loupe file creation.
 *      No alignment family layer is included (pairwise concept; to be added later).
 *
 *  Inputs:
 *      run_id:       Timestamp of the nextflow process
 *      all_info:     Tuple of [ids, annos, h5ads] — collected lists of all species IDs,
 *                    annotation column names, and raw h5ad paths
 *      samap_output: Path to the SAMap pickle file (RUN_SAMAP output, uncleaned)
 *
 *  Outputs:
 *      Tuple of [ids, counts, barcodes, features, umap, meta] for CREATE_ALL_LOUPE_FILE
 *      and a logfile.
 */

process CREATE_ALL_LOUPE_INPUT {
    tag "${run_id} - merge raw counts matrices for all-species Loupe"

    container 'mdiblbiocore/postanalysis:latest'

    input:
        val run_id
        tuple val(ids), val(annos), val(h5ads)
        path samap_output

    output:
        tuple val(ids),
            path("${ids.join('_')}_counts.mtx"),
            path("${ids.join('_')}_barcodes.csv"),
            path("${ids.join('_')}_features.csv"),
            path("${ids.join('_')}_umap.csv"),
            path("${ids.join('_')}_meta.csv"),
            emit: loupeinput
        path "${run_id}_create_all_loupe_input.log", emit: logfile

    script:
    """
    LOG="${run_id}_create_all_loupe_input.log"
    create_all_loupe_input.py \\
        --samap         ${samap_output} \\
        --h5ads         ${h5ads.join(' ')} \\
        --species       ${ids.join(' ')} \\
        --celltype-cols ${annos.join(' ')} \\
        2>&1 | tee -a \$LOG
    """
}
