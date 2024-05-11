process PLOT_FIG_03 {
    container 'leonhafner/plotting'

    publishDir "${params.output}", mode: 'copy'

    input:
    val correlations

    output:
    path 'Fig_03.png'

    script:
    correlation_string = correlations.collect{ it.toString() }.join(';')
    template 'plot_fig_03.R'

    stub:
    """
    touch Fig_03.png
    """
}