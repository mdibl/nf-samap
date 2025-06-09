#!/usr/env/bin nextflow

// Import the required modules 
include { build_blast_db } from './modules/build_blast_db.nf'
include { run_blast } from './modules/run_blast.nf'
include { RUN_SAMAP } from './modules/run_samap.nf'
include { visualize_results } from './modules/visualize_results.nf'

// Define the main workflow
workflow {

    Channel
        .fromPath('data')
        .set { data_dir }


    RUN_SAMAP(
        file('samap_input.json'),
        data_dir
    )

}