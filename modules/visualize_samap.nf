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
 *      outdir:         Directory in which to save final output
 *
 *  Outputs:
 *      Several visualizations about the SAMap results and a logfile
 *      results/${run_di}/plots/chord.html
 *      results/${run_id}/plots/sankey.html
 *      results/${run_id}/plots/scatter.png
 *      results/${run_id}/csv/hms.csv 
 *      results/${run_id}/csv/pms.csv 
 */

process VISUALIZE_SAMAP {
    tag "${run_id} - SAMap visualization"

    container 'mdiblbiocore/samap:latest'

    input:
        val run_id
        path samap_obj
        tuple val(id2), val(anno)

    output:
        path "chord.html"
        path "sankey.html"
        path "scatter.png"
        path "hms.csv"
        path "pms.csv"
        path "${run_id}_viz.log"

    script:
    """
    LOG="${run_id}_viz.log"
    visualize_samap.py --input ${samap_obj} --id2 ${id2} --annotation ${anno} 2>&1 | tee -a \$LOG
    """
}
