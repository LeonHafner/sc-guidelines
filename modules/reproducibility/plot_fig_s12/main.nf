process PLOT_FIG_S12 {
    
    container 'leonhafner/python'

    publishDir "${params.output}/Fig_S12", mode: 'copy'
    
    input:
    path pvalues

    output:
    path "Fig_S12.png"

    script:
    template 'plot_fig_s12.py'

    stub:
    """
    touch Fig_S12.png
    """
}