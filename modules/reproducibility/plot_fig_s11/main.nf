process PLOT_FIG_S11 {
    
    container 'leonhafner/plotting_rev'

    publishDir "${params.output}/Fig_S11", mode: 'copy'
    
    input:
    path pvalues

    output:
    path "Fig_S11.png"

    script:
    template 'plot_fig_s11.R'

    stub:
    """
    touch Fig_S11.png
    """
}