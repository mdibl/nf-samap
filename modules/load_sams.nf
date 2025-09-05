/*
 *  MODULE: load_sams.nf
 *
 *  Description: 
 *      Loads SAM objects using h5ad files described 
 *      in the sample sheet.
 *
 *  Inputs:
 *      run_id:         Timestamp of the nextflow process
 *      sample_sheet:   Path to the sample sheet containing the sample metadata
 *      data_dir:       Staging the data directory so the script can access it
 *
 *  Outputs:
 *      One pickled SAM object per sample and a logfile.
 *      results/run_id/sams/id2.pkl
 *      results/run_id/logs/run_id_load_sams.log
 */

process LOAD_SAMS {
    tag "${run_id} - load and pickle SAM objects"

    container 'mdiblbiocore/samap:latest'

    input:
        val run_id
        tuple val(meta), val(h5ad) 

    // this might be the only module where I actually need to stage the data_dir, since the list of paths won't be recognized by memVerge and automatically staged
    // since I need them all to be strings inside the list, bc can't really make LOAD_SAMS work with the channel one-by-one standard

    output:
        path "*.pkl", emit: sams
        path "${run_id}_load_sams.log", emit: logfile

    script:
    """
    LOG="${run_id}_load_sams.log"
    load_sams.py \\
        --id2 ${meta[0].join(' ')} \\
        --h5ad ${h5ad.join(' ')} \\
        2>&1 | tee -a \$LOG
    """
}
