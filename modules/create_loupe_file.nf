/*
 *  MODULE: create_loupe_file.nf
 *
 *  Description: 
 *      Takes a combined anndata object with associated ids and generates the shared feature space Loupe file
 *
 *  Inputs:
 *      run_id:         Timestamp of the nextflow process
 *      SeuratObject:   Paths referring to the Seurat Objects
 *      annotation:     Values of inputted annotation layer inputted in the Sample Sheet
 *
 *  Outputs:
 *      An Obs, Var, and Counts object for every sample and a logfile.
 *      results/run_id/logs/run_id_preprocess_seurat_object.log
 */

process CREATE_LOUPE_FILE {
    tag "${run_id} - create shared feature space Loupe file"

    container 'mdiblbiocore/loupe:latest'

    input:
        val run_id
        tuple val(id1), val(id2), path(combined_anndata)


    output:
        path "${id1}_${id2}_LoupeCreation.log", emit: logfile

    script:
    """  
    LOG="${run_id}_${id1}_${id2}_preprocess_seurat_object.log"
        Rscript /usr/local/bin/create_loupe_file.R \
        --id1 ${id1} \
        --id2 ${id2} \
        --anno ${combined_anndata} 2>&1 | tee -a \$LOG
    """
}