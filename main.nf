#!/usr/env/bin nextflow

/*
 *  PIPELINE: main.nf
 *
 *  Description:
 *      SAMap-based cross-species transcriptome mapping pipeline.
 *      Performs preprocessing of input metadata, reciprocal BLAST between
 *      species, SAMap alignment, and post-analysis of results.
 *
 *  Inputs:
 *      - sample_sheet.csv                 Sample metadata sheet
 *
 *  Workflow Overview:
 *      1. Preprocess sample sheet to coerce input data into appropriate format (.rds -> .h5ad, if needed)
 *      2. Generate all unordered species pairs
 *      3. Run reciprocal BLAST on each species pair, if BLAST maps are not specified by user beforehand
 *      4. Build a SAMap object
 *      5. Run SAMap on each unordered unique species pair
 *      6. Analyze SAMap alignment output
 *
 *  Parameters:
 *      --run_id        Run ID provided by user. If none is provided a timestamp is used. Default: null
 *      --sample_sheet  Path to the sample sheet provided by user. Default: 'sample_sheet.csv'
 *      --maps_dir      (optional) Path to a directory containing precomputed BLAST maps if any are provided. 
 *                      Any value other than null will skip the BLAST module. Default: null
 *      --outdir   The directory all where all results will be stored. S3 paths are supported, but require a mmc.config file with appropriate key setup. Default: '.'
 *
 *  Outputs:
 *      - outdir/preprocess_seurat_object/                                              Extracted counts matrix, obs, and feature tables per species
 *      - outdir/preprocess_anndata_object/                                             AnnData (.h5ad) files per species
 *      - outdir/blast/maps/{id1}{id2}/                                                 Reciprocal BLAST result files per species pair (skipped if maps_dir set)
 *      - outdir/load_sams/                                                             Pickled SAM objects per species
 *      - outdir/build_samap/samap.pkl                                                  SAMAP object (pre-algorithm)
 *      - outdir/run_samap/samap_results.pkl                                            SAMAP object (post-algorithm)
 *      - outdir/Analysis/AnalysisResults/{id1}-{id2}/GenePairs.csv                     Top gene-pair mappings with alignment scores
 *      - outdir/Analysis/AnalysisResults/{id1}-{id2}/pms_cluster_alignment_scores.csv. Pairwise cell-cluster alignment scores
 *      - outdir/Analysis/AnalysisResults/{id1}-{id2}/Grouping_Analysis/                Per-alignment-family differential expression results
 *      - outdir/{module}/run_id_{module}.log                                           Per-module logfiles (one per module per run)
 *
 *  Author:     Markus Sujansky, Ryan Sonderman
 *  Created:    2025-06-12
 *  Last Modified: 2026-05-04
 *  Version:    2.0.2
 */

// Import the required modules 
include { PREPROCESS_SEURAT_OBJECT } from './modules/preprocess_seurat_object.nf'
include { PREPROCESS_ANNDATA_OBJECT } from './modules/preprocess_anndata_object.nf'
include { RUN_BLAST_PAIR } from './modules/run_blast_pair.nf'
include { LOAD_SAMS } from './modules/load_sams.nf'
include { BUILD_SAMAP } from './modules/build_samap.nf'
include { RUN_SAMAP } from './modules/run_samap.nf'
include { VISUALIZE_SAMAP } from './modules/visualize_samap.nf'
include { SUMMARY_SAMAP } from './modules/summary_samap.nf'
include { CONNECTED_DE } from './modules/connected_de.nf'
include { ADDITIONAL_ANALYSIS } from './modules/additional_analysis.nf'
include { validateParameters; paramsHelp; samplesheetToList } from 'plugin/nf-schema'

workflow {
    // Generate run ID unless one is provided
    run_id = params.run_id ?: "${new Date().format('yyyyMMdd_HHmmss')}"
    run_id_ch = channel.value(run_id)

    // Stage static input files
    sample_sheet = channel.fromPath(params.sample_sheet)

    // Validation of necessary files
    if (!new File(params.sample_sheet).exists()) {
        error "Missing required file: sample sheet '${params.sample_sheet}'"
    }
    
    // Reformat Sample_Sheet to remove necessity of Sample_Sheet for downstream processes


    sample_sheet
        .map { file -> 
            samplesheetToList(file.toString(), "./nf-samap/assets/schema_input.json")
        }
        .flatMap { entries -> entries }
        .map { entry ->
            def meta          = entry[0]
            def counts           = entry[1]
            def transcriptome = entry.size() > 2 ? entry[2] : null

            if ((!params.maps_dir || !new File(params.maps_dir).exists()) && transcriptome == null) {
                error "Sample '${meta.id}': Input for BLAST (prot/transcriptome) must be provided when maps_dir is not precomputed and/or not set"
            }
            if(meta.type == "prot" && (!meta.mapping_dict || !new File(params.mapping_dict).exists())) {
                log.warn "Careful! You provided a proteome as input to Sample ${meta.id} but didn't provide a valid mapping
                dictionary to convery back to Gene Ids/Symbols. Be sure your data features are in the correct format!"
            }
            return [meta, counts, transcriptome]
        }
        .set { ch_samples }
        

        

    // Grab all input expression data paths to extract relevant info
    ch_samples
        .map { tuple ->
            def (meta, exprData, fasta) = tuple
            return [meta, exprData]
        }
        .set{ expr }


    expr
        .map { meta, data ->
            def isSeurat  = data.toString() ==~ /(?i).*\.(rds|rdata|rda)$/
            def isAnndata = data.toString() ==~ /(?i).*\.(h5ad|h5|loom)$/

            if (!isSeurat && !isAnndata) {
                error "Sample '${meta.id}': '${data}' is not a recognized Seurat or AnnData format"
            }

            return [meta, data, isSeurat ? "seurat" : "anndata"]
        }
        .branch { meta, data, format ->
            seurat:  format == "seurat"
            anndata: format == "anndata"
        }
        .set { ch_expr_branched }

    
     // Seurat Obj -> Component Parts
    PREPROCESS_SEURAT_OBJECT(
        run_id_ch,
        ch_expr_branched.seurat.map { meta, data, format -> [meta, data] }
    )

    //Component Partx -> AnnData Object
    PREPROCESS_ANNDATA_OBJECT(
        run_id_ch,
        PREPROCESS_SEURAT_OBJECT.out.seurat_data
    )
    anndata = PREPROCESS_ANNDATA_OBJECT.out.anndata
        .mix(
            ch_expr_branched.anndata.map { meta, data, format -> [meta.id, data] }
        )
 
    //SAMap-sepcific AnnData preprocessing module HERE, to be applied to all h5ad objects, not just those created in previous modules!!!

    // Generate unique unordered sample pairs
    pairs_channel = ch_samples
        .combine(ch_samples)
        .filter { a, _b, _c, d, _e, _f -> a.id < d.id }

    // Run BLAST or load precomputed map files 
   if (params.maps_dir) {
        maps_dir = channel.fromPath(params.maps_dir)
    }else {
        RUN_BLAST_PAIR(
            run_id_ch,
            pairs_channel
        )
        maps_dir = RUN_BLAST_PAIR.out.maps
    }

    // Join anndata with map_dict so all three fields stay associated
    sample_info = anndata
        .join(
            ch_samples.map { meta, so, fasta -> [meta.id, meta.map_dict.exists() ? meta.map_dict : null] }
        )
        // each element: [id, h5ad, map_dict]

    // Collect from the SAME channel — ordering is guaranteed consistent
    collected_ids      = sample_info.map { it[0] }.collect()
    collected_h5ads    = sample_info.map { it[1] }.collect()
    collected_mapdicts = sample_info.map { it[2] }.collect()

    condensedSampleSheet = collected_ids
        .map { ids -> [ids, collected_h5ads.getVal()] }

    LOAD_SAMS(
        run_id_ch,
        condensedSampleSheet
    )
    sams = LOAD_SAMS.out.sams

    //Make mapping_dict optional for 
    BUILD_SAMAP(
        run_id_ch,
        condensedSampleSheet,
        maps_dir,
        sams,
        collected_mapdicts
    )
    samap = BUILD_SAMAP.out.samap

    // Run SAMap on the SAMAP object to generate mapping results
    RUN_SAMAP(
        run_id_ch,
        samap
    )
    samap_results = RUN_SAMAP.out.results

    // Building channel obj for visualization module
    anno = ch_samples
    .map { tuple ->
        def (meta, h5ad, fasta) = tuple
        return meta.annotation
    }
    .collect()

    annotations = collected_ids
        .map { ids -> [ ids, anno.getVal()] }

    // Visualize the SAMap results
      VISUALIZE_SAMAP(
        run_id_ch,
        samap_results,
        annotations
    ) 

    // Necessary Context
    id_anno = ch_samples
        .map { meta, h5ad, fasta -> [meta.id, meta.annotation] }

    //Make id + Annotation pairs for pairwise post-analysis
    idCompare = id_anno.combine(id_anno)
        .filter { a_id, a_anno, b_id, b_anno -> a_id < b_id }
        .map { a_id, a_anno, b_id, b_anno ->
            [a_id, b_id, a_anno, b_anno]
    }

    all_ids   = id_anno.map { id, _annot -> id }.collect()
    all_annos = id_anno.map { _id, annot -> annot }.collect()
    id_anno_collected = id_anno
        .collect()
        .map { pairs -> 
            [pairs.collate(2).collect { it[0] }, pairs.collate(2).collect { it[1] }]
        }
    idCompare_samap = idCompare.combine(samap_results)
    idCompare_samap.view()


    SUMMARY_SAMAP(
        run_id_ch,
        idCompare_samap,
        id_anno_collected
    )
    genepairs = SUMMARY_SAMAP.out.genepairs
    cleaned = SUMMARY_SAMAP.out.samap_cleaned  // keep as tuple val(id1), val(id2), path(pkl)
    pms     = SUMMARY_SAMAP.out.pms            // keep as tuple val(id1), val(id2), path(pms)

    // Join pms and cleaned to idCompare by id1+id2 keys
    idCompare_with_inputs = idCompare
        .join(pms,     by: [0, 1])  // match on id1, id2
        .join(cleaned, by: [0, 1])  // match on id1, id2

    CONNECTED_DE(
        run_id_ch,
        idCompare_with_inputs
    )
    groupinganalysis = CONNECTED_DE.out.groupinganalysis
    //Grouping_Analysis.view()
    analysis = CONNECTED_DE.out.analysis
    //analysis.view()
    all_de_results = CONNECTED_DE.out.all_de_results
    //all_de_results.view()


    // Since genepairs, pms come from a different module, how can I make sure that the right order of
    // pms, genepairs, GroupingAnalysis, analysis, and all_de_results are passed into the next module correctly?
    combined = CONNECTED_DE.out.groupinganalysis
        .join(CONNECTED_DE.out.analysis,        by: [0,1])
        .join(CONNECTED_DE.out.all_de_results,  by: [0,1])
        .join(SUMMARY_SAMAP.out.genepairs,      by: [0,1])
        .join(SUMMARY_SAMAP.out.pms,            by: [0,1])
        .join(idCompare,                        by: [0,1])

    combined.view()


    ADDITIONAL_ANALYSIS(
        run_id_ch,
        combined
    )



    //CONSOLIDATION MODULE, make this contingent on >2 input species, otherwise skip

     NONPAIRWISE_ANALYSIS(
        id_anno
    )  

} 
