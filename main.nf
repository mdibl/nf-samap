#!/usr/env/bin nextflow

/*
 *  PIPELINE: main.nf
 *
 *  Description:
 *      SAMap-based cross-species transcriptome mapping pipeline.
 *      Performs preprocessing of input metadata, reciprocal BLAST between
 *      species, SAMap alignment, and visualization of results.
 *
 *  Workflow Overview:
 *      1. Preprocess sample sheet to classify inputs and assign IDs
 *      2. Generate all unordered species pairs
 *      3. Run reciprocal BLAST on each species pair
 *      4. Build a SAMap object
 *      5. Run SAMap on each BLAST result
 *      6. Visualize SAMap alignment and write results
 *
 *  Inputs:
 *      - data/transcriptomes/*.fasta      Input transcriptome FASTA files
 *      - data/*.h5ad                      Precomputed AnnData files
 *      - sample_sheet.csv                 Sample metadata sheet
 *
 *  Parameters:
 *      --use_precomputed_blast=BOOL   If true, skips BLAST and uses precomputed maps. Default: false
 *
 *  Outputs:
 *      - results/${run_id}/sample_sheet.csv       Updated metadata with type and ID
 *      - results/${run_id}/maps/                  Reciprocal BLAST outputs per sample pair
 *      - results/${run_id}/samap_objects/         Pickled SAMap object (Python)
 *      - results/${run_id}/plots/sankey.html      Mapping Sankey diagram
 *      - results/${run_id}/plots/scatter.png      Gene mapping scatter plot
 *      - results/${run_id}/csv/hms.csv            Highest mapping scores (HMS) matrix
 *      - results/${run_id}/csv/pms.csv            Pairwise mapping scores (PMS) matrix
 *      - results/${run_id}/logs/vis.log           Log file for visualization step
 *
 *  Author:     Ryan Sonderman
 *  Created:    2025-06-12
 *  Version:    1.0.0
 */

// Import the required modules 
include { PREPROCESS } from './modules/preprocess.nf'
include { RUN_BLAST_PAIR } from './modules/run_blast_pair.nf'
include { LOAD_SAMS } from './modules/load_sams.nf'
include { BUILD_SAMAP } from './modules/build_samap.nf'
include { RUN_SAMAP } from './modules/run_samap.nf'
include { VISUALIZE_SAMAP } from './modules/visualize_samap.nf'

workflow {
    // Generate a workflow id from a timestamp
    run_id = "${new Date().format('yyyyMMdd_HHmmss')}"
    run_id_ch = Channel.value(run_id)


    // Stage static input files
    data_dir    = Channel.fromPath('data')
    results_dir = Channel.fromPath('results')
    sample_sheet = Channel.fromPath('sample_sheet.csv')


    // Preprocess sample sheet to add type and ID
    sample_sheet_pr = PREPROCESS(
        run_id_ch,
        sample_sheet,
        data_dir,
    )


    // Generate unique unordered sample pairs
    samples_channel = sample_sheet_pr.splitCsv(header: true, sep: ',')
    pairs_channel = samples_channel
        .combine(samples_channel)
        .filter { a, b -> a.id2 < b.id2 }


    // Run BLAST or load precomputed map files 
    if (params.use_precomputed_blast) {
        // Set path to maps from provided data
        maps_dir = data_dir.map { it -> it.resolve('maps/') } 
    } else {
        RUN_BLAST_PAIR(
            run_id_ch,
            pairs_channel,
            data_dir.first(),
        )
        // Set path to maps from BLAST results
        maps_dir = results_dir.map { it -> it.resolve('maps/') } 
    }


    // Load SAM objects from the AnnData h5ad files
    sams = LOAD_SAMS(
        run_id_ch,
        sample_sheet_pr,
        data_dir,
    )


    // Build the SAMap object from the SAM objects and the BLAST maps
    samap = BUILD_SAMAP(
        run_id_ch,
        sample_sheet_pr,
        data_dir,
        maps_dir,
        results_dir,
        sams,
    )


    // Run SAMap on the SAMAP object to generate mapping results
    samap_results = RUN_SAMAP(
        run_id_ch,
        results_dir,
        samap
    )


    // Visualize the SAMap results
    VISUALIZE_SAMAP(
        run_id_ch,
        samap_results,
        sample_sheet_pr,
    )
}
