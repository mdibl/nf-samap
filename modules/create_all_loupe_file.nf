/*
 *  MODULE: create_all_loupe_file.nf
 *
 *  Description:
 *      Builds a shared feature space Loupe file covering all N species
 *      from SAMap UMAP embeddings and combined raw count matrices.
 *
 *  Inputs:
 *      run_id:     Timestamp of the nextflow process
 *      loupeinput: Tuple of [ids, counts, barcodes, features, umap, meta]
 *                  produced by CREATE_ALL_LOUPE_INPUT
 *
 *  Outputs:
 *      A single .cloupe file for all species and a logfile.
 */

process CREATE_ALL_LOUPE_FILE {
    tag "${run_id} - create all-species Loupe file"

    container 'mdiblbiocore/loupe:latest'

    input:
        val run_id
        tuple val(ids),
            path(counts),
            path(barcodes),
            path(features),
            path(umap),
            path(meta)

    output:
        path "${ids.join('_')}.cloupe",             emit: loupe
        path "${run_id}_create_all_loupe.log",      emit: logfile

    script:
    """
    LOG="${run_id}_create_all_loupe.log"
    create_all_loupe_file.R \\
        --ids      ${ids.join(' ')} \\
        --counts   ${counts} \\
        --barcodes ${barcodes} \\
        --features ${features} \\
        --umap     ${umap} \\
        --meta     ${meta} \\
        2>&1 | tee -a \$LOG
    """
}
