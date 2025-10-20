/*
 *  MODULE: preprocess_anndata_object.nf
 *
 *  Description: 
 *      Takes a channel of species IDs, compressed counts matrices, .obs files, and .feats files to create an Anndata (h5ad) Object
 *
 *  Inputs:
 *      run_id:   Timestamp of the nextflow process
 *      Counts:   A sparse matrix containing 
 *      Obs:      A .csv file containing Cell IDs, original identities, and the desired annotation layer (I think we can add more than one annotation
 *                layer and just specify one)
 *      Feats:    A .csv file containing necessary feature data
 *
 *  Outputs:
 *      An AnnData object for each sample containing all necessary information for SAMap to run
 *      results/run_id/logs/run_id_preprocess_seurat_object.log
 */

process PREPROCESS_ANNDATA_OBJECT {
    tag "${run_id} - use extracted Seurat Object information to build h5ad object"

    container 'mdiblbiocore/samap:latest'

    input:
        val run_id
        tuple val(id), 
          path(counts), 
          path(obs), 
          path(feats)


    output: 
        tuple val(id), path("${meta.id}_preprocessed.h5ad"),
        emit: anndata
        path "${run_id}_preprocess_anndata_object.log", emit: logfile

    script:
    """  
    LOG="${run_id}_preprocess_anndata_object.log"
        /usr/local/bin/preprocess_anndata_object.py \
        --counts ${counts} \
        --obs ${obs} \
        --feats ${feats} \
        --id ${id} 2>&1 | tee -a \$LOG
    """
}