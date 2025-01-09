process PLOT_FIG_SZ {
    
    container 'leonhafner/python'

    publishDir "${params.output}/Fig_SZ", mode: 'copy'
    
    input:
    path pvalues

    output:
    path "Fig_SZ.png"

    script:
    template 'plot_fig_sz.py'

    stub:
    """
    touch Fig_SZ.png
    """
}