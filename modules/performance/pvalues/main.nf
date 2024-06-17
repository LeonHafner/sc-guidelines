process PVALUES {
    tag "${meta.scenario}_${meta.run}"

    container 'leonhafner/python'

    input:
    tuple val(meta), path(results)

    output:
    tuple val(meta), path("pval_${meta.scenario}_${meta.run}.tsv")

    script:
    meta_string = meta.join(';')
    path_string = results.join(';')
    meta = [scenario: meta[0].scenario, run: meta[0].run]
    template 'pvalues.py'

    stub:
    meta = [scenario: meta[0].scenario, run: meta[0].run]
    """
    touch pval_${meta.scenario}_${meta.run}.tsv
    """
}