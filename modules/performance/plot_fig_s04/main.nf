process PLOT_FIG_S04 {
    container 'leonhafner/plotting'

    publishDir "${params.output}/Fig_S04", mode: 'copy'

    input:
    tuple val(meta), path(auc_files)

    output:
    path 'Fig_S04.png'

    script:
    meta_string = meta.join(';')
    path_string = auc_files.join(';')
    template 'plot_fig_s04.R'

    stub:
    """
    touch Fig_S04.png
    """
}