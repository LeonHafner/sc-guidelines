process SCVI {
    tag "${meta.scenario}_${meta.run}"

    container 'leonhafner/scvi'

    input:
    tuple val(meta), path(input_anndata, stageAs: 'input_anndata.h5ad')

    output:
    tuple val(meta), path("result-*.tsv")

    script:
    meta = meta + [method: 'scvi']
    template 'scvi-de.py'

    stub:
    meta = meta + [method: 'scvi']
    """
    touch result-0_01.tsv
    touch result-0_02.tsv
    touch result-0_03.tsv
    touch result-0_04.tsv
    touch result-0_05.tsv
    touch result-0_06.tsv
    touch result-0_07.tsv
    """
}