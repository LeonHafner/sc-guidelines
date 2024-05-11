process PREPROCESSING {
    tag "${meta.scenario}_${meta.n_cells}_${meta.n_genes}_${meta.run}"
    
    container 'leonhafner/python'

    input:
    tuple val(meta), path(anndata)
    val threshold

    output:
    tuple val(meta), path('preprocessed.h5ad'), path('time.txt')

    script:
    """
    START=\$( date +%s )
    preprocessing.py \
        --input ${anndata} \
        --output preprocessed.h5ad \
        --threshold ${threshold}
    END=\$( date +%s )
    echo \$((\$END-\$START)) > time.txt
    """

    stub:
    """
    touch preprocessed.h5ad
    echo 43 > time.txt
    """
}