process FILTER_HVG {
    tag "${meta.scenario}_${meta.run}"
    
    container 'leonhafner/python'

    input:
    tuple val(meta), path(input_anndata, stageAs: 'input_anndata.h5ad')
    val hvg_ratio

    output:
    tuple val(meta), path("${meta.scenario}_${meta.run}.h5ad")

    script:
    template 'filter_hvg.py'

    stub:
    """
    touch ${meta.scenario}_${meta.run}.h5ad
    """
}