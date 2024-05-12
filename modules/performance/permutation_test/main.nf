process PERMUTATION_TEST {
    tag "${meta.scenario}_${meta.run}"

    container 'leonhafner/python'

    input:
    tuple val(meta), path(input_anndata, stageAs: 'input_anndata.h5ad')

    output:
    tuple val(meta), path("${meta.scenario}_permutation-test_${meta.run}.tsv")

    script:
    meta = meta + [method: 'permutation-test']
    """
    permutation.py \
        --input ${input_anndata} \
        --output ${meta.scenario}_permutation-test_${meta.run}.tsv \
        -n 10000 \
        -n_max 100000 \
        --scenario ${meta.scenario}
    """

    stub:
    meta = meta + [method: 'permutation-test']
    """
    touch ${meta.scenario}_permutation-test_${meta.run}.tsv
    """
}