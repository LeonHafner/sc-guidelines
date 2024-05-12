process MAST {
    tag "${meta.scenario}_${meta.run}"

    container 'leonhafner/mast'

    input:
    tuple val(meta), path(input_anndata, stageAs: 'input_anndata.h5ad')

    output:
    tuple val(meta), path("${meta.scenario}_mast_${meta.run}.tsv")

    script:
    meta = meta + [method: 'mast']
    """
    mast.R \
        --input ${input_anndata} \
        --output ${meta.scenario}_mast_${meta.run}.tsv \
        --scenario ${meta.scenario}
    """

    stub:
    meta = meta + [method: 'mast']
    """
    touch ${meta.scenario}_mast_${meta.run}.tsv
    """
}