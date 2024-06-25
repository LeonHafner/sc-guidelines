process PLOT_FIG_05 {
    container 'leonhafner/plotting'

    publishDir "${params.output}/Fig_05", mode: 'copy'

    input:
    tuple val(meta), path(plotting_data)

    output:
    path "Fig_05.png"

    script:
    template 'plot_fig_05.R'

    stub:
    """
    touch Fig_05.png
    """
}