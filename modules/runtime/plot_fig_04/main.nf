process PLOT_FIG_04 {
    container 'leonhafner/plotting'

    publishDir "${params.output}", mode: 'copy'

    input:
    path fixed_cells
    path fixed_genes

    output:
    path 'Fig_04.png'

    script:
    template 'plot_fig_04.R'

    stub:
    """
    touch Fig_04.png
    """
}