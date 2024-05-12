process DESEQ2 {
    tag "${meta.scenario}_${meta.run}"

    container 'leonhafner/deseq2'

    input:
    tuple val(meta), path(input_anndata, stageAs: 'input_anndata.h5ad')

    output:
    tuple val(meta), path("${meta.scenario}_deseq2_${meta.run}.tsv")

    script:
    meta = meta + [method: 'deseq2']
    """
    deseq2.R \
        --input ${input_anndata} \
        --output ${meta.scenario}_deseq2_${meta.run}.tsv \
        --scenario ${meta.scenario}
    """

    stub:
    meta = meta + [method: 'deseq2']
    """
    touch ${meta.scenario}_deseq2_${meta.run}.tsv
    """
}