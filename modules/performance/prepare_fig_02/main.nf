process PREPARE_FIG_02 {
    container 'leonhafner/python'

    input:
    tuple val(meta), path(path_list)

    output:
    tuple val(meta), path('cell_counts.txt')

    script:
    template 'prepare_fig_02.py'

    stub:
    """
    touch cell_counts.txt
    """
}