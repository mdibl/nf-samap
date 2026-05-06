/*
 *  MODULE: preprocess_sam_object.nf
 *
 *  Description: 
 *      Takes a channel of species IDs and initialized AnnData Objects to process into a preprocessed Anndata (h5ad) Object
 *
 *  Inputs:
 *      run_id:   Timestamp of the nextflow process
 *      Anndata:   An initialized AnnData Object
 *
 *  Outputs:
 *      An AnnData object for each sample containing all necessary information for SAMap to run
 *      results/run_id/logs/run_id_preprocess_sam_object.log
 */

process PREPROCESS_SAM_OBJECT {
    tag "${run_id} - use initialized Anndata Object to build processed h5ad object"

    container 'mdiblbiocore/preprocessing:latest'

    input:
        val run_id
        tuple val(id), path(Anndata)


    output: 
        tuple val(id), path("${id}_preprocessed.h5ad"), emit: anndata
        path "${run_id}_${id}_preprocess_sam_object.log", emit: logfile

    script:
    """  
    LOG="${run_id}_${id}_preprocess_sam_object.log"
        /usr/local/bin/preprocess_sam_object.py \
        --anndata ${Anndata} \
        --id ${id} 2>&1 | tee -a \$LOG
    """
}