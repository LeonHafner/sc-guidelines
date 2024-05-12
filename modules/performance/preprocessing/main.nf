process PREPROCESSING {
    tag "${meta.scenario}_${meta.run}"

    container 'leonhafner/python'

    input:
    tuple val(meta), path(input_anndata, stageAs: 'input_anndata.h5ad')
    val threshold

    output:
    tuple val(meta), path("${meta.scenario}_${meta.run}.h5ad")

    script:
    """
    preprocessing.py \
        --input ${input_anndata} \
        --output ${meta.scenario}_${meta.run}.h5ad \
        --threshold ${threshold}
    """

    stub:
    """
    touch ${meta.scenario}_${meta.run}.h5ad
    """
}