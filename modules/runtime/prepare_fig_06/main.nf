process PREPARE_FIG_06 {
    container 'leonhafner/python'

    input:
    val meta

    output:
    path 'fixed_cells.tsv', emit: fixed_cells
    path 'fixed_genes.tsv', emit: fixed_genes

    script:
    meta_string = meta.collect{ it.toString() }.join(';')
    template 'prepare_fig_06.py'

    stub:
    """
    touch fixed_cells.tsv
    touch fixed_genes.tsv
    """
}