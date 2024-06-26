process PLOT_FIG_S06 {
    tag "run_${meta_prc.run}"

    container 'leonhafner/plotting'

    publishDir "${params.output}/Fig_S06", mode: 'copy'

    input:
    tuple val(meta_prc), path(prc_files)
    tuple val(meta_auc), path(auc_files)

    output:
    path "Fig_S06_run${meta_prc.run}.png"
    
    script:
    prc_path_string = prc_files.join(';')
    auc_meta_string = meta_auc.join(';')
    auc_path_string = auc_files.join(';')
    template 'plot_fig_s06.R'

    stub:
    """
    touch Fig_S06_run${meta_prc.run}.png
    """
}