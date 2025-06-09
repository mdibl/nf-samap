#!/usr/env/bin nextflow

// Import the required modules 
include { RUN_SAMAP } from './modules/run_samap.nf'

// Define the main workflow
workflow {

    // Stage the data files and config JSON
    data_dir = Channel.fromPath('data')
    config_file = Channel.fromPath('config.json')
    
    RUN_SAMAP(
        config_file,
        data_dir
    )
}