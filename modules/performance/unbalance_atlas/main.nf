process UNBALANCE_ATLAS {
    tag "${meta.scenario}_${meta.run}"

    container 'leonhafner/python'

    input:
    tuple val(meta), path(input_anndata, stageAs: 'input_anndata.h5ad')

    output:
    tuple val(meta), path("${meta.scenario}_${meta.run}.h5ad")

    script:
    template 'unbalance_atlas.py'

    stub:
    """
    touch ${meta.scenario}_${meta.run}.h5ad
    """
}