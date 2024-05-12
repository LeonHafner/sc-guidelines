process PLOT_FIG_07 {
    container 'leonhafner/plotting'

    publishDir "${params.output}", mode: 'copy'

    input:
    tuple val(meta), path(plotting_data)

    output:
    path "Fig_07.png"

    script:
    template 'plot_fig_07.R'

    stub:
    """
    touch Fig_07.png
    """
}