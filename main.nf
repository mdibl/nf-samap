#!/usr/env/bin nextflow

// Import the required modules 
include { PREPROCESS } from './modules/preprocess.nf'
include { RUN_BLAST_PAIR } from './modules/run_blast_pair.nf'
include { LOAD_SAMS } from './modules/load_sams.nf'
include { RUN_SAMAP } from './modules/run_samap.nf'
include { VISUALIZE_SAMAP } from './modules/visualize_samap.nf'

// Define the main workflow
workflow {

    // Stage the data files and config JSON
    data_dir    = Channel.fromPath('data')
    results_dir = Channel.fromPath('results')
    precomputed_dir = Channel.fromPath('precomputed')

    /*
    Preprocess the sample sheet by 
    1. Classifying the fasta files as `prot` or `nucl`
    2. Appending a unique 2 character ID to each sample
    
    Output is a CSV file with the following columns:
    - id: Unique ID for the sample
    - h5ad: Path to the h5ad file
    - fasta: Path to the fasta file
    - annotation: How to interpret the output of SAMap
    - type: Type of the fasta file (prot or nucl)
    - id2: Unique 2 character ID for the sample
    */
    sample_sheet = Channel.fromPath('sample_sheet.csv')
    sample_sheet_pr = PREPROCESS(
        sample_sheet,
        data_dir,
    )

    /*
    Create inordered pairs of samples from the sample sheet.
    These inordered pairs will be used to run BLAST.
    
    Output is a channel of tuples where each tuple contains two unique samples.
    */
    samples_channel = sample_sheet_pr.splitCsv(header: true, sep: ',')
    pairs_channel = samples_channel
        .combine(samples_channel)
        .filter { a, b -> a.id2 < b.id2 }

    
    /*
    Use the script provided from SAMap to run BLAST on the pairs of samples.
    
    Output is a blast mapping for each pair of samples.
    Mappings are stored in the `results/maps` directory with the format:
    maps/<sample1_id><sample2_id/<sample2_id>_to_<sample1_id>.txt
    where <sample1_id> and <sample2_id> are the unique 2 character IDs of the samples.
    */
    if (params.use_precomputed_blast) {
        samap_map_parent = precomputed_dir
        } else {
        RUN_BLAST_PAIR(
            pairs_channel,
            data_dir.first()
        )
        samap_map_parent = results_dir
    }

    /*
    Load the SAM objects from the h5ad files.

    Output is a channel of SAM objects.
    Each SAM object is a pickled file stored in the `results/sams` directory.
    */
    sams = LOAD_SAMS(
        data_dir,
        sample_sheet_pr
    )
}
