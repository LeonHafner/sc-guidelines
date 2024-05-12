process PLOT_FIG_S05 {
    container 'leonhafner/plotting'

    publishDir "${params.output}", mode: 'copy'

    input:
    tuple val(meta), path(plotting_data_kang)
    tuple val(meta2), path(plotting_data_atlas_negative)

    output:
    path 'Fig_S05.png'

    script:
    template 'plot_fig_s05.R'

    stub:
    """
    touch Fig_S05.png
    """
}