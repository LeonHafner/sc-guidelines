process PVALUES {
    tag "${tuple_list[0][0].scenario}_${tuple_list[0][0].run}"

    container 'leonhafner/python'

    input:
    val tuple_list

    output:
    tuple val(meta), path("pval_${meta.scenario}_${meta.run}.tsv")

    script:
    // Process channel input which is a list of tuples consisting of a meta map and a path object
    def meta_list = []
    def path_list = []
    for (element in tuple_list) {
        meta_list.add(element[0])
        path_list.add(element[1])
    }
    def meta_string = meta_list.join(';')
    def path_string = path_list.join(';')
    meta = [scenario: tuple_list[0][0].scenario, run: tuple_list[0][0].run]
    template 'pvalues.py'

    stub:
    meta = [scenario: tuple_list[0][0].scenario, run: tuple_list[0][0].run]
    """
    touch pval_${meta.scenario}_${meta.run}.tsv
    """
}