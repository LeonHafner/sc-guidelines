process PLOT_FIG_S09 {
    tag "${meta.scenario}_${meta.run}"
    
    container 'leonhafner/plotting'

    publishDir "${params.output}/Fig_S09", mode: 'copy'

    input:
    tuple val(meta), path(input_anndata)

    output:
    path "Fig_S09_run${meta.run}.png"

    script:
    template 'plot_fig_s09.R'

    stub:
    """
    touch Fig_S09_run${meta.run}.png
    """
}