process SCDD {
    tag "${meta.scenario}_${meta.run}"

    container 'leonhafner/scdd'

    input:
    tuple val(meta), path(input_anndata, stageAs: 'input_anndata.h5ad')

    output:
    tuple val(meta), path("${meta.scenario}_scdd_${meta.run}.tsv")

    script:
    meta = meta + [method: 'scdd']
    """
    scdd.R \
        --input ${input_anndata} \
        --output ${meta.scenario}_scdd_${meta.run}.tsv \
        --n_cores ${task.cpus}
    """

    stub:
    meta = meta + [method: 'scdd']
    """
    touch ${meta.scenario}_scdd_${meta.run}.tsv
    """
}