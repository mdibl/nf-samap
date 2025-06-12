process PREPROCESS {
    tag "Preprocess sample sheet"

    publishDir('results', mode: 'copy', pattern: '*.csv')

    container 'ryansonder/samap-blast:latest'

    input:
        path sample_sheet
        path data_dir

    output:       
        path 'sample_sheet_*.csv'

    script:
    """
    update_sample_sheet.sh ${sample_sheet}
    """
}
