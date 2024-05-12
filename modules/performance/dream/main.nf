process DREAM {
    tag "${meta.scenario}_${meta.run}"

    container 'leonhafner/dream'

    input:
    tuple val(meta), path(input_anndata, stageAs: 'input_anndata.h5ad')

    output:
    tuple val(meta), path("${meta.scenario}_dream_${meta.run}.tsv")

    script:
    meta = meta + [method: 'dream']
    """
    dream.R \
        --input ${input_anndata} \
        --output ${meta.scenario}_dream_${meta.run}.tsv \
        --scenario ${meta.scenario}
    """

    stub:
    meta = meta + [method: 'dream']
    """
    touch ${meta.scenario}_dream_${meta.run}.tsv
    """
}