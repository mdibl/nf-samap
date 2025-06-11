process VISUALIZE_SAMAP {
    tag "SAMap visualization"
    publishDir('results/plots', mode: 'copy', pattern: '*.html')
    publishDir('results/plots', mode: 'copy', pattern: '*.png')
    publishDir('results/csv', mode: 'copy', pattern: '*.csv')
    publishDir('results/logs', mode: 'copy', pattern: '*.log')

    container 'ryansonder/samap:latest'

    input:
        path samap_obj

    output:
        path 'sankey_*.html'
        path 'scatter_*.png'
        path 'hms_*.csv'
        path 'pms_*.csv'
        path 'visualization_*.log'

    script:
    """
    visualize_samap.py --input ${samap_obj}
    """
}
