process PLOT_FIG_S07 {
    tag "run_${meta.run}"

    container 'leonhafner/plotting'

    publishDir "${params.output}/Fig_S07", mode: 'copy'

    input:
    tuple val(meta), val(prc)
    val auc

    output:
    path "Fig_S07_run${meta.run}.png"
    
    script:
    // Preprocess items for prc channel
    prc_path_string = prc.join(';')

    // Preprocess items for auc channel
    auc_meta_list = []
    auc_path_list = []
    for (element in auc) {
        auc_meta_list.add(element[0])
        auc_path_list.add(element[1])
    }
    auc_meta_string = auc_meta_list.join(';')
    auc_path_string = auc_path_list.join(';')
    template 'plot_fig_s07.R'

    stub:
    """
    touch Fig_S07_run${meta.run}.png
    """
}