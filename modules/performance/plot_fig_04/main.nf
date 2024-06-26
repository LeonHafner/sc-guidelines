process PLOT_FIG_04 {
    container 'leonhafner/plotting'

    publishDir "${params.output}/Fig_04", mode: 'copy'

    input:
    tuple val(meta), path(prc)

    output:
    path "Fig_04.png"

    script:
    prc_path_string = prc.join(';')
    template 'plot_fig_04.R'

    stub:
    """
    touch Fig_04.png
    """
}