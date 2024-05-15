process PSEUDOBULKING {
    tag "${meta.scenario}_${meta.n_cells}_${meta.n_genes}_${meta.run}"
    
    container 'leonhafner/python'

    input:
    tuple val(meta), path(anndata)

    output:
    tuple val(meta), path('pseudobulked.h5ad'), path('time.txt')

    script:
    """
    START=\$( date +%s )
    pseudobulking.py \
        --input ${anndata} \
        --output pseudobulked.h5ad \
        --scenario dataset
    END=\$( date +%s )
    echo \$((\$END-\$START)) > time.txt
    """

    stub:
    """
    touch pseudobulked.h5ad
    echo 44 > time.txt
    """
}