process PLOT_FIG_S02 {
    container 'leonhafner/plotting'

    publishDir "${params.output}/Fig_S02", mode: 'copy'

    input:
    tuple val(meta), path(auc_files)

    output:
    path "Fig_S02.png"

    script:
    meta_string = meta.join(';')
    path_string = auc_files.join(';')
    template 'plot_fig_s02.R'

    stub:
    """
    touch Fig_S02.png
    """
}