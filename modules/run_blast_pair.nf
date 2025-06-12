process RUN_BLAST_PAIR {
    tag "${a.id2}_vs_${b.id2}"

    publishDir('results', mode: 'copy')

    container 'ryansonder/samap-blast:latest'

    input:
        tuple val(a), val(b)
        path data_dir

    output:
        path "maps/*/*_to_*.txt"

    script:
    """
    head -n 10 ${a.fasta} > a_trunc.fasta
    head -n 10 ${b.fasta} > b_trunc.fasta
    tr1=a_trunc.fasta
    tr2=b_trunc.fasta

    echo Running BLAST for ${a.id2} vs ${b.id2}
    map_genes.sh \\
        --tr1 \$tr1 --t1 ${a.type} --n1 ${a.id2} \\
        --tr2 \$tr2 --t2 ${b.type} --n2 ${b.id2}
    """
}
