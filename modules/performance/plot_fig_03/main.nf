process PLOT_FIG_03 {
    container 'leonhafner/plotting'

    publishDir "${params.output}/Fig_03", mode: 'copy'

    input:
    tuple val(meta), path(auc_files)

    output:
    path 'Fig_03.png'

    script:
    meta_string = meta.join(';')
    path_string = auc_files.join(';')
    template 'plot_fig_03.R'

    stub:
    """
    touch Fig_03.png
    """
}