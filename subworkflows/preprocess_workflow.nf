include { PREPROCESS_SEURAT_OBJECT } from '../modules/preprocess_seurat_object.nf'
include { PREPROCESS_ANNDATA_OBJECT } from '../modules/preprocess_anndata_object.nf'
include { PREPROCESS_SAM_OBJECT } from '../modules/preprocess_sam_object.nf'

workflow PREPROCESSING_WORKFLOW {
    take:
    expression // Input channel with names
    run_id

    main:
    expression
        .map { meta, data ->
            def isSeurat  = data.toString() ==~ /(?i).*\.(rds|rdata|rda)$/
            def isAnndata = data.toString() ==~ /(?i).*\.(h5ad|h5|loom)$/

            if (!isSeurat && !isAnndata) {
                error "Sample '${meta.id}': '${data}' is not a recognized Seurat or AnnData format"
            }

            return [meta, data, isSeurat ? "seurat" : "anndata"]
        }
        .branch { _meta, _data, format ->
            seurat:  format == "seurat"
            anndata: format == "anndata"
        }
        .set { ch_expr_branched }
    
     // Seurat Obj -> Component Parts
    PREPROCESS_SEURAT_OBJECT(
        run_id,
        ch_expr_branched.seurat.map { meta, data, _format -> [meta, data] }
    )

    //Component Parts -> AnnData Object
    PREPROCESS_ANNDATA_OBJECT(
        run_id,
        PREPROCESS_SEURAT_OBJECT.out.seurat_data
    )
    PREPROCESS_ANNDATA_OBJECT.out.anndata
        .mix(
            ch_expr_branched.anndata.map { meta, data, _format -> [meta.id, data] }
        )
        .toSortedList { a, b -> a[0] <=> b[0] }
        .flatMap { entries -> entries }
        .set { anndata }
 
    //Initialized AnnData -> Preprocessed AnnData
    PREPROCESS_SAM_OBJECT(
        run_id,
        anndata
    )

    emit:
    processed_AnnData = PREPROCESS_SAM_OBJECT.out.anndata // Original greetings
    raw_ad = anndata
}