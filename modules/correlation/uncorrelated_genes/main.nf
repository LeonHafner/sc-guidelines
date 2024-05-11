process UNCORRELATED_GENES {
    tag "${meta.scenario}_${meta.run}"

    container 'leonhafner/python'

    input:
    tuple val(meta), path(anndata)

    output:
    tuple val(meta), path('uncorrelated_genes.h5ad')

    script:
    template 'uncorrelated_genes.py'

    stub:
    """
    touch uncorrelated_genes.h5ad
    """
}