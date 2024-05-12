process PLOT_FIG_02 {
    container 'leonhafner/plotting'

    publishDir "${params.output}/Fig_02", mode: 'copy'

    input:
    tuple val(meta), path(cell_counts)

    output:
    path 'Fig_02.png', emit: 'png1'
    path 'Fig_02A.png', emit: 'png2'
    path 'Fig_02A.drawio', emit: 'drawio'

    script:
    template 'plot_fig_02.R'

    stub:
    """
    touch Fig_02.png
    touch Fig_02A.png
    touch Fig_02A.drawio
    """
}