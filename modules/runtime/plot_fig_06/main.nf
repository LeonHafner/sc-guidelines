process PLOT_FIG_06 {
    container 'leonhafner/plotting'

    publishDir "${params.output}/Fig_06", mode: 'copy'

    input:
    path fixed_cells
    path fixed_genes

    output:
    path 'Fig_06.png'

    script:
    template 'plot_fig_06.R'

    stub:
    """
    touch Fig_06.png
    """
}