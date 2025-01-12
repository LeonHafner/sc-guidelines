process PLOT_FIG_S10 {
    
    container 'leonhafner/plotting_rev'

    publishDir "${params.output}/Fig_S10", mode: 'copy'
    
    input:
    path pvalues

    output:
    path "Fig_S10.png"

    script:
    template 'plot_fig_s10.R'

    stub:
    """
    touch Fig_S10.png
    """
}