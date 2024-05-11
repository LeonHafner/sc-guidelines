process CORRELATIONS {
    tag "${meta.scenario}_${meta.run}"

    container 'leonhafner/python'

    input:
    tuple val(meta), path(anndata)

    output:
    tuple val(meta), path('correlations.tsv')

    script:
    template 'correlations.py'

    stub:
    """
    touch correlations.tsv
    """
}