process PLOT_FIG_S08 {
    container 'leonhafner/plotting'

    publishDir "${params.output}", mode: 'copy'

    input:
    tuple val(meta), path(auc_files)

    output:
    path "Fig_S08.png"

    script:
    meta_string = meta.join(';')
    path_string = auc_files.join(';')
    template 'plot_fig_s08.R'

    stub:
    """
    touch Fig_S08.png
    """
}