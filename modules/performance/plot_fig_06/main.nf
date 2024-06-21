process PLOT_FIG_06 {
    container 'leonhafner/plotting'

    publishDir "${params.output}/Fig_06", mode: 'copy'

    input:
    tuple val(meta), path(prc)

    output:
    path "Fig_06.png"

    script:
    prc_path_string = prc.join(';')
    template 'plot_fig_06.R'

    stub:
    """
    touch Fig_06.png
    """
}