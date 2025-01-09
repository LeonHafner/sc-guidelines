process PLOT_FIG_SY {
    
    container 'leonhafner/plotting_rev'

    publishDir "${params.output}/Fig_SY", mode: 'copy'
    
    input:
    path pvalues

    output:
    path "Fig_SY.png"

    script:
    template 'plot_fig_sy.R'

    stub:
    """
    touch Fig_SY.png
    """
}