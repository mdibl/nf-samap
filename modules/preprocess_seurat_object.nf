/*
 *  MODULE: preprocess_seurat_object.nf
 *
 *  Description: 
 *      Takes a channel of Seurat Objects and desired annotation layer, and extracts all necessary information for constructing the h5ad objects
 *
 *  Inputs:
 *      run_id:         Timestamp of the nextflow process
 *      SeuratObject:   Paths referring to the Seurat Objects
 *      annotation:     Values of inputted annotation layer inputted in the Sample Sheet
 *
 *  Outputs:
 *      An Obs, Var, and Counts object for every sample and a logfile.
 *      results/run_id/samap_objects/run_id_samap.pkl
 *      results/run_id/logs/run_id_preprocess_seurat_object.log
 */

process PREPROCESS_SEURAT_OBJECT {
    tag "${run_id} - extract important Seurat Object information"

    container 'mdiblbiocore/samap:latest'

    input:
        val run_id
        tuple val(meta), path(SO)


    output:
        tuple val(meta.id), 
          path("${meta.id}_Counts.mtx"), 
          path("${meta.id}_Obs.csv"), 
          path("${meta.id}_Feats.csv"), 
          emit: sample_data
        path "${run_id}_preprocess_seurat_object.log", emit: logfile

    script:
    """  
    LOG="${run_id}_preprocess_seurat_object.log"
        Rscript /usr/local/bin/preprocess_seurat_object.R \
        --so ${SO} \
        --id ${meta.id} \
        --anno ${meta.annotation} 2>&1 | tee -a \$LOG
    """
}