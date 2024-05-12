process PLOT_FIG_S06 {
    container 'leonhafner/plotting'

    publishDir "${params.output}", mode: 'copy'

    input:
    val tuple_list

    output:
    path 'Fig_S06.png'

    script:
    meta_list = []
    path_list = []
    for (element in tuple_list) {
        meta_list.add(element[0])
        path_list.add(element[1])
    }
    meta_string = meta_list.join(';')
    path_string = path_list.join(';')
    template 'plot_fig_s06.R'

    stub:
    """
    touch Fig_S06.png
    """
}