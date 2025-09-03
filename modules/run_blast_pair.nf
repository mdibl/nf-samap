/*
 *  MODULE: run_blast_pair.nf
 *
 *  Description: 
 *      Uses bash script provided from 
 *      https://github.com/atarashansky/SAMap
 *      on pairs of samples.
 *
 *  Inputs:
 *      run_id:         Timestamp of the nextflow process
 *      val(a):         Map with keys: fasta, type, id2 (Sample A metadata)
 *      val(b):         Map with keys: fasta, type, id2 (Sample B metadata)
 *      data_dir:       Staging the data directory so the script can access it
 *
 *  Outputs:
 *      Two BLAST result text files for each direction and a logfile
 *      results/${run_id}/maps/{pair_id}/[A_to_B.txt, B_to_A.txt]
 */

process RUN_BLAST_PAIR {
    stageInMode 'copy'
    tag "${run_id} - ${meta1.id2}_vs_${meta2.id2}"

    container 'mdiblbiocore/samap-blast:latest'

    input:
        val run_id
        tuple val(meta1), val(meta2), path(fasta1), path(fasta2)


    output:
        path "maps/*/*_to_*.txt", emit: maps
        path "${run_id}_${meta1.id2}${meta2.id2}_blast.log", emit: logfile

    script:
    """
    LOG="${run_id}_${meta1.id2}${meta2.id2}_blast.log"
    map_genes.sh \\
        --threads ${task.cpus} \\
        --tr1 ${fasta1} --t1 ${meta1.type} --n1 ${meta1.id2} \\
        --tr2 ${fasta2} --t2 ${meta2.type} --n2 ${meta2.id2} | \\
        while IFS= read -r line; do
            echo "[\$(date +'%Y-%m-%d %H:%M:%S.%3N')] \$line"
        done 2>&1 | tee -a \$LOG
    """
}
