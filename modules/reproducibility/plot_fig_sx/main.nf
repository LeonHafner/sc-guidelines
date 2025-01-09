process PLOT_FIG_SX {
    
    container 'leonhafner/plotting_rev'

    publishDir "${params.output}/Fig_SX", mode: 'copy'
    
    input:
    path pvalues

    output:
    path "Fig_SX.png"

    script:
    template 'plot_fig_sx.R'

    stub:
    """
    touch Fig_SX.png
    """
}