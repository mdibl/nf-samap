#!/usr/env/bin nextflow

// Import the required modules 
include { RUN_SAMAP } from './modules/run_samap.nf'

// Define the main workflow
workflow {

    RUN_SAMAP(
        file('config.json')
    )

}