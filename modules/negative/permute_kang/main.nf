process PERMUTE_KANG {
    container 'leonhafner/python'

    input:
    tuple val(meta), path(input_anndata, stageAs: 'input_anndata.h5ad')

    output:
    tuple val(meta), path('kang2018.h5ad')

    script:
    template 'permute_kang.py'

    stub:
    """
    touch kang2018.h5ad
    """
}