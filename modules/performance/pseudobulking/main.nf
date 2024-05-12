process PSEUDOBULKING {
    tag "${meta.scenario}_${meta.run}"

    container 'leonhafner/python'

    input:
    tuple val(meta), path(input_anndata, stageAs: 'input_anndata')

    output:
    tuple val(meta), path("${meta.scenario}_${meta.run}.h5ad")

    script:
    """
    pseudobulking.py \
        --input ${input_anndata} \
        --output ${meta.scenario}_${meta.run}.h5ad \
        --scenario ${meta.scenario}
    """

    stub:
    """
    touch ${meta.scenario}_${meta.run}.h5ad
    """
}