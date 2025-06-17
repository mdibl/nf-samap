/*
 *  MODULE: preprocess.nf
 *
 *  Description: 
 *      Updates the sample sheet to add:
 *          - [type] Classification of fasta (protein or nucleotide)
 *          - [id2] Unique 2-character ID
 *
 *  Inputs:
 *      run_id:         Timestamp of the nextflow process
 *      sample_sheet:   Path to the sample sheet being updated
 *      data_dir:       Staging the data directory so the script can access it
 *
 *  Outputs:
 *      results/run_id/run_id_sample_sheet.csv
 */

process PREPROCESS {
    tag "${run_id} - sample sheet preprocessing"

    publishDir("results/${run_id}/", mode: 'copy', pattern: '*.csv')

    container 'ryansonder/samap-blast:latest'

    input:
        val run_id
        path sample_sheet
        path data_dir

    output:       
        path "${run_id}_sample_sheet.csv"

    script:
    """
    update_sample_sheet.sh ${sample_sheet}
    if [[ ! -f out.csv ]]; then
        echo "ERROR: out.csv not found!" >&2
        exit 1
    fi
    mv out.csv ${run_id}_sample_sheet.csv
    """
}
