process HIERARCHICAL_BOOTSTRAPPING {
    tag "${meta.scenario}_${meta.run}"

    container 'leonhafner/python'

    input:
    tuple val(meta), path(input_anndata, stageAs: 'input_anndata.h5ad')

    output:
    tuple val(meta), path("${meta.scenario}_hierarchical-bootstrapping_${meta.run}.tsv")

    script:
    meta = meta + [method: 'hierarchical-bootstrapping']
    """
    hierarchical_bootstrapping.py \
        --input ${input_anndata} \
        --output1 ${meta.scenario}_hierarchical-bootstrapping_${meta.run}_pre.tsv \
        --output2 ${meta.scenario}_hierarchical-bootstrapping_${meta.run}.tsv \
        --scenario ${meta.scenario} 
    """

    stub:
    meta = meta + [method: 'hierarchical-bootstrapping']
    """
    touch ${meta.scenario}_hierarchical-bootstrapping_${meta.run}.tsv
    """
}