/*
 *  MODULE: visualize_samap.nf
 *
 *  Description: 
 *      Produces several visualizations from analysis.py from 
 *      https://github.com/atarashansky/SAMap
 *
 *  Inputs:
 *      run_id:         Timestamp of the nextflow process
 *      samap_object:   Channel containing a pickled SAMAP object
 *      sample_sheet:   Path to the sample sheet CSV with sample metadata
 *
 *  Outputs:
 *      Several visualizations about the SAMap results
 *      results/${run_di}/plots/chord.html
 *      results/${run_id}/plots/sankey.html
 *      results/${run_id}/plots/scatter.png
 *      results/${run_id}/csv/hms.csv 
 *      results/${run_id}/csv/pms.csv 
 */

process VISUALIZE_SAMAP {
    tag "${run_id} - SAMap visualization"

    publishDir("results/${run_id}/plots/", mode: 'copy', pattern: '*.html')
    publishDir("results/${run_id}/plots/", mode: 'copy', pattern: '*.png')
    publishDir("results/${run_id}/csv/", mode: 'copy', pattern: '*.csv')

    container 'pipeline/samap:latest'

    input:
        val run_id
        path samap_obj
        path sample_sheet

    output:
        path 'chord.html'
        path 'sankey.html'
        path 'scatter.png'
        path 'hms.csv'
        path 'pms.csv'

    script:
    """
    visualize_samap.py --input ${samap_obj} --sample-sheet ${sample_sheet}
    """
}
