process DISTINCT {
    tag "${meta.scenario}_${meta.run}"

    container 'leonhafner/distinct'
    
    input:
    tuple val(meta), path(input_anndata, stageAs: 'input_anndata.h5ad')

    output:
    tuple val(meta), path("${meta.scenario}_distinct_${meta.run}.tsv")

    script:
    meta = meta + [method: 'distinct']
    """
    distinct.R \
        --input ${input_anndata} \
        --output ${meta.scenario}_distinct_${meta.run}.tsv \
        --scenario ${meta.scenario}
    """

    stub:
    meta = meta + [method: 'distinct']
    """
    touch ${meta.scenario}_distinct_${meta.run}.tsv
    """
}