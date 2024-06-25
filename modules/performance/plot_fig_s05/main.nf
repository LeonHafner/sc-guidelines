process PLOT_FIG_S05 {
    container 'leonhafner/plotting'

    publishDir "${params.output}/Fig_S05", mode: 'copy'

    input:
    tuple val(meta), path(prc_files)

    output:
    path 'Fig_S05.png'

    script:
    meta_string = meta.join(';')
    path_string = prc_files.join(';')
    template 'plot_fig_s05.R'

    stub:
    """
    touch Fig_S05.png
    """
}