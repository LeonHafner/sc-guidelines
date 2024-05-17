process PLOT_FIG_05 {
    container 'leonhafner/plotting'

    publishDir "${params.output}/Fig_05", mode: 'copy'

    input:
    tuple val(meta), path(auc_files)

    output:
    path 'Fig_05.png'

    script:
    meta_string = meta.join(';')
    path_string = auc_files.join(';')
    template 'plot_fig_05.R'

    stub:
    """
    touch Fig_05.png
    """
}