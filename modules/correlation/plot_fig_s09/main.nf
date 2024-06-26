process PLOT_FIG_S09 {
    container 'leonhafner/plotting'

    publishDir "${params.output}/Fig_S09", mode: 'copy'

    input:
    path correlations, stageAs: 'correlations/correlations_?.tsv'

    output:
    path 'Fig_S09.png'

    script:
    correlation_string = correlations.collect{ it.toString() }.join(';')
    template 'plot_fig_s09.R'

    stub:
    """
    touch Fig_S09.png
    """
}