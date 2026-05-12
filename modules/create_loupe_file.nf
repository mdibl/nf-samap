/*
 *  MODULE: create_loupe.nf
 *
 *  Description:
 *      Builds a shared feature space Loupe file from SAMap output
 *      using raw count matrices and SAMap UMAP embeddings.
 *
 *  Inputs:
 *      run_id:     Timestamp of the nextflow process
 *      loupeinput: Tuple of [id1, id2, counts, barcodes, features, umap, meta]
 *
 *  Outputs:
 *      A .cloupe file for each species pair and a logfile.
 */

process CREATE_LOUPE {
    tag "${run_id} - create loupe file for ${id1}-${id2}"

    container 'mdiblbiocore/loupe:latest'

    input:
        val run_id
        tuple val(id1), val(id2),
            path(counts),
            path(barcodes),
            path(features),
            path(umap),
            path(meta)

    output:
        path "${id1}_${id2}.cloupe",                       emit: loupe
        path "${run_id}_${id1}_${id2}_create_loupe.log",   emit: logfile

    script:
    """
    LOG="${run_id}_${id1}_${id2}_create_loupe.log"
    create_loupe.R \\
        --id1      ${id1} \\
        --id2      ${id2} \\
        --counts   ${counts} \\
        --barcodes ${barcodes} \\
        --features ${features} \\
        --umap     ${umap} \\
        --meta     ${meta} \\
        2>&1 | tee -a \$LOG
    """
}