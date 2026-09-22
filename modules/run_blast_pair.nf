/*
 *  MODULE: run_blast_pair.nf
 *
 *  Description:
 *      Runs reciprocal BLAST between two species using the map_genes.sh
 *      script from https://github.com/atarashansky/SAMap. Builds a BLAST
 *      database for each species and aligns in both directions
 *      (blastp for proteomes, tblastx/blastn for transcriptomes).
 *
 *  Inputs:
 *      val run_id:     Run ID (user-provided or timestamp)
 *      tuple:
 *          val(id1):       Species A ID (2-letter SAMap species code)
 *          val(type1):     Species A sequence type ('prot' or 'nucl')
 *          path(fasta1):   Species A proteome/transcriptome FASTA
 *          val(id2):       Species B ID (2-letter SAMap species code)
 *          val(type2):     Species B sequence type ('prot' or 'nucl')
 *          path(fasta2):   Species B proteome/transcriptome FASTA
 *
 *      Note: main.nf reshapes each pair from
 *      [meta1, so1, fasta1, meta2, so2, fasta2] into this tuple
 *      before calling the module. Pairs are unordered and unique (id1 < id2).
 *
 *  Outputs:
 *      maps:       Reciprocal BLAST results for the pair
 *                  maps/{id1}{id2}/{id1}_to_{id2}.txt
 *                  maps/{id1}{id2}/{id2}_to_{id1}.txt
 *      logfile:    Timestamped map_genes.sh output
 *                  {run_id}_{id1}{id2}_blast.log
 */

process RUN_BLAST_PAIR {
    tag "${run_id} - ${id1}_vs_${id2}"

    container 'mdiblbiocore/samap-blast:latest'

    input:
        val run_id
        tuple val(id1), val(type1), path(fasta1), val(id2), val(type2), path(fasta2)


    output:
        path "maps/*/*_to_*.txt", emit: maps
        path "${run_id}_${id1}${id2}_blast.log", emit: logfile

    script:
    """
    LOG="${run_id}_${id1}${id2}_blast.log"
    map_genes.sh \\
        --threads ${task.cpus} \\
        --tr1 ${fasta1} --t1 ${type1} --n1 ${id1} \\
        --tr2 ${fasta2} --t2 ${type2} --n2 ${id2} | \\
        while IFS= read -r line; do
            echo "[\$(date +'%Y-%m-%d %H:%M:%S.%3N')] \$line"
        done 2>&1 | tee -a \$LOG
    """
}