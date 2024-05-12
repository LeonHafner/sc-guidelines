process PLOT_FIG_05 {
    container 'leonhafner/plotting'

    publishDir "${params.output}/Fig_05", mode: 'copy'

    input:
    val tuple_list

    output:
    path 'Fig_05.png'

    script:
    // Process channel input which is a list of tuples consisting of a meta map and a path object
    meta_list = []
    path_list = []
    for (element in tuple_list) {
        meta_list.add(element[0])
        path_list.add(element[1])
    }
    meta_string = meta_list.join(';')
    path_string = path_list.join(';')
    template 'plot_fig_05.R'

    stub:
    """
    touch Fig_05.png
    """
}