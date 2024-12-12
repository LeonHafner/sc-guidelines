process TTEST {
    tag "${meta.scenario}_${meta.run}"

    container 'leonhafner/python'

    input:
    tuple val(meta), path(input_anndata, stageAs: 'input_anndata.h5ad')

    output:
    tuple val(meta), path("${meta.scenario}_ttest_${meta.run}.tsv")

    script:
    meta = meta + [method: 'ttest']
    """
    ttest.py \
        --input ${input_anndata} \
        --output ${meta.scenario}_ttest_${meta.run}.tsv
    """

    stub:
    meta = meta + [method: 'ttest']
    """
    touch ${meta.scenario}_ttest_${meta.run}.tsv
    """
}