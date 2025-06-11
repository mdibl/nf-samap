#!/usr/env/bin nextflow

// Import the required modules 
include { RUN_BLAST } from './modules/run_blast.nf'
include { RUN_SAMAP } from './modules/run_samap.nf'
include { VISUALIZE_SAMAP } from './modules/visualize_samap.nf'

// Define the main workflow
workflow {

    // Stage the data files and config JSON
    data_dir = Channel.fromPath('data')
    config_file = Channel.fromPath('config.json')
    keys_json = Channel.fromPath('keys.json')
    sample_sheet = Channel.fromPath('sample_sheet.csv')

    // Blast the data files
    RUN_BLAST(
        sample_sheet
    )

    // Testing Exit
    exit 0

    // Run the SAMAP process with the provided data and config
    samap_obj = RUN_SAMAP(
        config_file,
        data_dir
    )

    // Visualize the results of the SAMAP process
    VISUALIZE_SAMAP(
        samap_obj,
        keys_json
    )
}