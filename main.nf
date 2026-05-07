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
include { PREPROCESSING_WORKFLOW } from './subworkflows/preprocess_workflow.nf'
include { PAIRWISE_ANALYSIS } from './subworkflows/pairwise_analysis.nf'
include { CREATE_LOUPE } from './subworkflows/create_loupe.nf'
include { RUN_BLAST_PAIR } from './modules/run_blast_pair.nf'
include { LOAD_SAMS } from './modules/load_sams.nf'
include { BUILD_SAMAP } from './modules/build_samap.nf'
include { RUN_SAMAP } from './modules/run_samap.nf'
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
            if(meta.type == "prot" && (!meta.map_dict || !new File(meta.map_dict).exists())) {
                log.warn "Careful! You provided a proteome as input to Sample ${meta.id} but didn't provide a valid mapping
                dictionary to convery back to Gene Ids/Symbols. Be sure your data features are in the correct format!"
            }
            return [meta, counts, transcriptome]
        }
        .set { ch_samples }
        

    // Grab all input expression data paths to extract relevant info
    ch_samples
        .map { tuple ->
            def (meta, exprData, _fasta) = tuple
            return [meta, exprData]
        }
        .set{ expr }


    //Preprocessing Subworkflow
    PREPROCESSING_WORKFLOW(
        expr,
        run_id_ch
    )    

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
     PREPROCESSING_WORKFLOW.out.processed_AnnData
        .join(
            ch_samples.map { meta, _so, _fasta -> 
                def mapDict = (meta.map_dict && new File(meta.map_dict.toString()).exists()) ? meta.map_dict : null
                [meta.id, mapDict]
            }        
        )
        .set{ sample_info }
        // each element: [id, h5ad, map_dict]

    // Collect from the SAME channel — ordering is guaranteed consistent
    collected_ids      = sample_info.map {id -> id[0] }.collect()
    collected_h5ads    = sample_info.map {ad -> ad[1] }.collect()
    collected_mapdicts = sample_info.map {md -> md[2] }.collect()

    condensedSampleSheet = collected_ids
        .map { ids -> [ids, collected_h5ads.getVal()] }

    LOAD_SAMS(
        run_id_ch,
        condensedSampleSheet
    )

    //Make mapping_dict optional for 
    BUILD_SAMAP(
        run_id_ch,
        condensedSampleSheet,
        maps_dir,
        LOAD_SAMS.out.sams,
        collected_mapdicts
    )

    // Run SAMap on the SAMAP object to generate mapping results
    RUN_SAMAP(
        run_id_ch,
        BUILD_SAMAP.out.samap
    )


    //Subworkflow for Pairiwise Post-Analysis
    PAIRWISE_ANALYSIS(
        run_id_ch,
        ch_samples,
        RUN_SAMAP.out.results
    )

    PAIRWISE_ANALYSIS.out.pairCompare
        .join(PREPROCESSING_WORKFLOW.out.raw_ad)
        .map { id1, id2, anno1, anno2, h5ad1 -> [id2, id1, anno1, anno2, h5ad1] }
        .join(PREPROCESSING_WORKFLOW.out.raw_ad)
        .map { id2, id1, anno1, anno2, h5ad1, h5ad2 -> [id1, id2, anno1, anno2, h5ad1, h5ad2] }
        .join(PAIRWISE_ANALYSIS.out.samap_cleaned, by: [0, 1])
        .set { idCompare_with_h5ads }

    if (params.create_loupe == "true") {
        CREATE_LOUPE(
            run_id_ch,
            idCompare_with_h5ads
        )
    }

    //CONSOLIDATION MODULE, make this contingent on >2 input species, otherwise skip

     /* NONPAIRWISE_ANALYSIS(
        id_anno
    )   */

} 
