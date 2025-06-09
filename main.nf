#!/usr/env/bin nextflow

// Import the required modules 
include { RUN_SAMAP } from './modules/run_samap.nf'
include { VISUALIZE_SAMAP } from './modules/visualize_samap.nf'

// Define the main workflow
workflow {

    // Stage the data files and config JSON
    data_dir = Channel.fromPath('data')
    config_file = Channel.fromPath('config.json')

    samap_obj = RUN_SAMAP(
        config_file,
        data_dir
    )

    VISUALIZE_SAMAP(
        samap_obj
    )
}