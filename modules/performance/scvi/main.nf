process SCVI {
    tag "${meta.scenario}_${meta.run}"

    container 'leonhafner/scvi'

    input:
    tuple val(meta), path(input_anndata, stageAs: 'input_anndata.h5ad')

    output:
    tuple val(meta), path("${meta.scenario}_scvi_${meta.run}.tsv")

    script:
    meta = meta + [method: 'scvi']
    """
    scvi-de.py \
        --input ${input_anndata} \
        --output ${meta.scenario}_scvi_${meta.run}.tsv \
        --scenario ${meta.scenario}
    """

    stub:
    meta = meta + [method: 'scvi']
    """
    touch ${meta.scenario}_scvi_${meta.run}.tsv
    """
}